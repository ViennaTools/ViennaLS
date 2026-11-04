"""Read OpenVDB files (.vdb) and turn grids of material labels into level sets.

Reading needs NumPy. Most .vdb files are compressed with Blosc, OpenVDB's
default, which additionally needs the blosc package::

    pip install numpy blosc

Example, for a grid of material IDs whose x axis points down into the
substrate::

    import viennals as vls
    import viennals.vdb

    grid = vls.vdb.readGrids("structure.vdb")[0]
    for label, levelSet in vls.vdb.levelSetsFromLabels(grid, 3, ("z", "y", "-x")):
        ...  # e.g. insert into a ViennaPS domain, from the bottom up
"""

import struct
import warnings
import zlib
from dataclasses import dataclass, field

import numpy as np

from . import BoundaryConditionEnum, Mesh, d2, d3

__all__ = ["Grid", "readGrids", "levelSetsFromLabels"]

_MAGIC = 0x56444220
# the oldest file format storing per-node compression metadata
_MIN_FILE_VERSION = 222
_COMPRESS_ZIP, _COMPRESS_ACTIVE_MASK, _COMPRESS_BLOSC = 0x1, 0x2, 0x4
_VALUE_TYPES = {"int32": "<i4", "int64": "<i8", "float": "<f4", "double": "<f8"}
_METADATA_TYPES = {
    "bool": "<?",
    "int32": "<i",
    "int64": "<q",
    "float": "<f",
    "double": "<d",
    "vec3i": "<3i",
    "vec3s": "<3f",
    "vec3d": "<3d",
}
# which inactive values a node stores in front of its values
_NO_MASK_AND_ONE_INACTIVE_VAL = 2
_MASK_AND_NO_INACTIVE_VALS = 3
_MASK_AND_ONE_INACTIVE_VAL = 4
_MASK_AND_TWO_INACTIVE_VALS = 5
_NO_MASK_AND_ALL_VALS = 6


@dataclass
class Grid:
    """A grid read from a .vdb file.

    indices holds the index coordinates of the active voxels, values their
    values; all other voxels have the background value. A voxel's centre lies
    at translation + voxelSize * index in world coordinates.
    """

    name: str
    valueType: str
    voxelSize: tuple
    translation: tuple
    background: object
    metadata: dict = field(repr=False)
    indices: np.ndarray = field(repr=False)
    values: np.ndarray = field(repr=False)


class _Stream:
    def __init__(self, data):
        self._data = data
        self.pos = 0

    def read(self, size):
        if size < 0 or self.pos + size > len(self._data):
            raise ValueError("unexpected end of the .vdb file")
        chunk = self._data[self.pos : self.pos + size]
        self.pos += size
        return chunk

    def unpack(self, fmt):
        values = struct.unpack(fmt, self.read(struct.calcsize(fmt)))
        return values[0] if len(values) == 1 else values

    def string(self):
        return self.read(self.unpack("<I")).decode("utf-8", "replace")

    def array(self, dtype, count):
        dtype = np.dtype(dtype)
        return np.frombuffer(self.read(dtype.itemsize * count), dtype)

    def mask(self, size):
        bits = np.unpackbits(self.array(np.uint8, size // 8), bitorder="little")
        return bits.astype(bool)


@dataclass
class _Tree:
    compression: int
    valueType: np.dtype
    storedType: np.dtype
    leaves: list = field(default_factory=list)
    tiles: list = field(default_factory=list)


def _readMetadata(stream):
    metadata = {}
    for _ in range(stream.unpack("<I")):
        name, typeName = stream.string(), stream.string()
        raw = stream.read(stream.unpack("<I"))
        if typeName == "string":
            metadata[name] = raw.decode("utf-8", "replace")
        elif typeName in _METADATA_TYPES:
            value = struct.unpack(_METADATA_TYPES[typeName], raw)
            metadata[name] = value[0] if len(value) == 1 else value
        # other types, e.g. OpenVDB's internal delayed-load data, are skipped
    return metadata


def _readTransform(stream):
    mapType = stream.string()
    if mapType == "TranslationMap":
        return (1.0, 1.0, 1.0), stream.unpack("<3d")
    if mapType in ("UniformScaleMap", "ScaleMap"):
        translation, values = (0.0, 0.0, 0.0), stream.unpack("<15d")
    elif mapType in ("UniformScaleTranslateMap", "ScaleTranslateMap"):
        values = stream.unpack("<18d")
        translation, values = values[:3], values[3:]
    else:
        raise ValueError(f"unsupported transform {mapType}")
    # the scale is followed by the voxel size and derived values
    scale, voxelSize = values[0:3], values[3:6]
    if not np.allclose(scale, voxelSize):
        raise ValueError(f"unexpected layout of the transform {mapType}")
    return scale, translation


def _bloscDecompress(data):
    try:
        import blosc
    except ImportError as error:
        raise ImportError(
            "This .vdb file is compressed with Blosc; install it with "
            "'pip install blosc'."
        ) from error
    return blosc.decompress(data)


def _readData(stream, count, tree):
    size = count * tree.storedType.itemsize
    if tree.compression & (_COMPRESS_BLOSC | _COMPRESS_ZIP):
        # a size of zero or below marks uncompressed data
        stored = stream.unpack("<q")
        if stored <= 0:
            raw = stream.read(-stored)
        elif tree.compression & _COMPRESS_BLOSC:
            raw = _bloscDecompress(stream.read(stored))
        else:
            raw = zlib.decompress(stream.read(stored))
    else:
        raw = stream.read(size)
    if len(raw) != size:
        raise ValueError("corrupt .vdb file: a block of values has the wrong size")
    return np.frombuffer(raw, tree.storedType).astype(tree.valueType)


def _readActiveValues(stream, size, valueMask, tree):
    """Values of a node's active entries, in order."""
    flag = stream.unpack("<b")
    if flag in (
        _NO_MASK_AND_ONE_INACTIVE_VAL,
        _MASK_AND_ONE_INACTIVE_VAL,
        _MASK_AND_TWO_INACTIVE_VALS,
    ):
        stream.read(tree.valueType.itemsize)
        if flag == _MASK_AND_TWO_INACTIVE_VALS:
            stream.read(tree.valueType.itemsize)
    if flag in (
        _MASK_AND_NO_INACTIVE_VALS,
        _MASK_AND_ONE_INACTIVE_VAL,
        _MASK_AND_TWO_INACTIVE_VALS,
    ):
        # selects between the two inactive values
        stream.read(size // 8)
    if tree.compression & _COMPRESS_ACTIVE_MASK and flag != _NO_MASK_AND_ALL_VALS:
        return _readData(stream, int(valueMask.sum()), tree)
    return _readData(stream, size, tree)[valueMask]


def _localCoordinates(offsets, log2Dim):
    offsets = np.asarray(offsets)
    low = (1 << log2Dim) - 1
    return np.stack(
        [offsets >> (2 * log2Dim), (offsets >> log2Dim) & low, offsets & low], -1
    )


def _readInternalNode(stream, origin, log2Dim, childLog2Dim, tree):
    size = 1 << (3 * log2Dim)
    childMask = stream.mask(size)
    valueMask = stream.mask(size)
    values = _readActiveValues(stream, size, valueMask, tree)
    tiles = valueMask & ~childMask
    if tiles.any():
        offsets = np.flatnonzero(valueMask)
        isTile = ~childMask[offsets]
        tree.tiles.append(
            (
                origin + (_localCoordinates(offsets[isTile], log2Dim) << childLog2Dim),
                1 << childLog2Dim,
                values[isTile],
            )
        )
    for offset in np.flatnonzero(childMask):
        childOrigin = origin + (_localCoordinates(offset, log2Dim) << childLog2Dim)
        if childLog2Dim == 3:
            # the leaf's value mask; it is stored again with its values
            stream.read(64)
            tree.leaves.append(childOrigin)
        else:
            _readInternalNode(stream, childOrigin, 4, 3, tree)


def _readTree(stream, blockPos, endPos, tree, name):
    if stream.unpack("<i") != 1:
        raise ValueError(f"grid {name}: trees with several buffers are not supported")
    background = stream.array(tree.valueType, 1)[0]
    numTiles, numChildren = stream.unpack("<2I")
    for _ in range(numTiles):
        stream.unpack("<3i")
        stream.array(tree.valueType, 1)
        if stream.unpack("<?"):
            raise ValueError(f"grid {name}: active root tiles are not supported")
    for _ in range(numChildren):
        _readInternalNode(stream, np.array(stream.unpack("<3i")), 5, 7, tree)
    if stream.pos != blockPos:
        raise ValueError(f"grid {name}: unexpected layout of the tree")

    indices, values = [], []
    for origin, size, tileValues in tree.tiles:
        voxels = np.indices((size,) * 3).reshape(3, -1).T
        for tileOrigin, value in zip(origin, tileValues):
            indices.append(tileOrigin + voxels)
            values.append(np.full(len(voxels), value, tree.valueType))
    for origin in tree.leaves:
        mask = stream.mask(512)
        values.append(_readActiveValues(stream, 512, mask, tree))
        indices.append(origin + _localCoordinates(np.flatnonzero(mask), 3))
    if stream.pos != endPos:
        raise ValueError(f"grid {name}: unexpected layout of the values")
    if not indices:
        return background, np.zeros((0, 3), np.int64), np.zeros(0, tree.valueType)
    return background, np.concatenate(indices), np.concatenate(values)


def readGrids(fileName):
    """Read the grids of scalar values (int32, int64, float, double) of a .vdb
    file. Grids of other types are skipped with a warning."""
    with open(fileName, "rb") as file:
        stream = _Stream(file.read())
    if stream.unpack("<q") != _MAGIC:
        raise ValueError(f"{fileName} is not a .vdb file")
    version = stream.unpack("<I")
    if version < _MIN_FILE_VERSION:
        raise ValueError(f"{fileName}: .vdb file format {version} is too old")
    stream.unpack("<2I")  # library version
    if not stream.unpack("<?"):
        raise ValueError(f"{fileName}: files without grid offsets are not supported")
    stream.read(36)  # UUID
    _readMetadata(stream)

    grids, trees = [], {}
    for _ in range(stream.unpack("<i")):
        uniqueName, gridType = stream.string(), stream.string()
        parent = stream.string()
        gridPos, blockPos, endPos = stream.unpack("<3q")
        name = uniqueName.split("\x1e")[0]

        half = gridType.endswith("_HalfFloat")
        parts = (gridType[: -len("_HalfFloat")] if half else gridType).split("_")
        if len(parts) != 5 or parts[0] != "Tree" or parts[2:] != ["5", "4", "3"]:
            warnings.warn(f"{fileName}: skipping grid {name} of type {gridType}")
            stream.pos = endPos
            continue
        if parts[1] not in _VALUE_TYPES:
            warnings.warn(f"{fileName}: skipping grid {name} of type {gridType}")
            stream.pos = endPos
            continue
        valueType = np.dtype(_VALUE_TYPES[parts[1]])

        stream.pos = gridPos
        compression = stream.unpack("<I")
        metadata = _readMetadata(stream)
        voxelSize, translation = _readTransform(stream)
        if parent:
            if parent not in trees:
                raise ValueError(f"{fileName}: grid {name} shares a missing tree")
            background, indices, values = trees[parent]
        else:
            tree = _Tree(compression, valueType, np.dtype("<f2") if half else valueType)
            background, indices, values = _readTree(
                stream, blockPos, endPos, tree, name
            )
            trees[uniqueName] = (background, indices, values)
        grids.append(
            Grid(
                name,
                parts[1],
                tuple(voxelSize),
                tuple(translation),
                background.item(),
                metadata,
                indices,
                values,
            )
        )
        stream.pos = endPos
    return grids


def _levelSet(indices, values, low, high, gridDelta, dimension):
    module = d3 if dimension == 3 else d2
    # FromMesh needs the points sorted with the last coordinate slowest
    order = np.lexsort(indices.T)
    coordinates = np.zeros((len(indices), 3))
    coordinates[:, :dimension] = indices[order] * gridDelta
    mesh = Mesh()
    for node in coordinates.tolist():
        mesh.insertNextNode(node)
    mesh.getPointData().insertNextScalarData(values[order].tolist(), "LSValues")

    # a quarter of a grid delta inside the outermost grid points, so that the
    # rounding of the bounds to grid points cannot add one
    bounds = []
    for d in range(dimension):
        bounds += [(low[d] + 0.25) * gridDelta, (high[d] - 0.25) * gridDelta]
    boundaryConditions = [BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * (dimension - 1)
    boundaryConditions.append(BoundaryConditionEnum.INFINITE_BOUNDARY)
    levelSet = module.Domain(bounds, boundaryConditions, gridDelta)
    module.FromMesh(levelSet, mesh).apply()
    # two layers of points, as the other importers create
    module.Expand(levelSet, 2).apply()
    return levelSet


def levelSetsFromLabels(grid, dimension=3, axes=None, ambient=None, order=None):
    """Level sets of the materials in a grid of integer material labels.

    The interfaces are the faces between voxels of different labels. The voxel
    centres become the grid points of the level sets, whose grid delta is the
    voxel size. Voxels with the background value or the label ambient are
    empty space.

    axes names the grid axis behind each level set axis, from the first to the
    vertical one, e.g. ("z", "y", "-x") for a grid whose x axis points down into
    the substrate; "-" flips an axis. In 2D, the grid axis left out is the one
    the structure is extruded along, and its middle slice is used. The default
    is ("x", "y", "z") or ("x", "y").

    order lists the labels from the bottom up; by default, they are sorted by
    their mean height.

    Returns (label, level set) pairs from the bottom up. Each level set
    encloses its own material and all materials below it, so they can be
    inserted into a ViennaPS domain in this order. The lateral boundaries are
    reflective and the vertical one is infinite: the materials continue past
    the lateral edges and below the bottom of the grid.
    """
    if dimension not in (2, 3):
        raise ValueError("dimension must be 2 or 3")
    if not np.issubdtype(grid.values.dtype, np.integer):
        raise TypeError(f"grid {grid.name} does not hold integer labels")
    voxelSize = np.asarray(grid.voxelSize, float)
    if not np.allclose(voxelSize, voxelSize[0]):
        raise ValueError(f"grid {grid.name} does not have cubic voxels")
    gridDelta = float(voxelSize[0])
    if len(grid.indices) == 0:
        raise ValueError(f"grid {grid.name} has no active voxels")

    low = grid.indices.min(0)
    labels = np.full(grid.indices.max(0) - low + 1, grid.background, grid.values.dtype)
    labels[tuple((grid.indices - low).T)] = grid.values
    shift = np.asarray(grid.translation) / gridDelta
    if np.allclose(shift, np.round(shift)):
        low = low + np.round(shift).astype(low.dtype)
    else:
        warnings.warn(
            f"grid {grid.name}: the translation is not a multiple of the voxel "
            "size and is ignored"
        )

    if axes is None:
        axes = ("x", "y", "z")[:dimension]
    if len(axes) != dimension:
        raise ValueError(f"axes must name {dimension} grid axes")
    source = ["xyz".index(axis.lstrip("+-")) for axis in axes]
    flipped = [axis.startswith("-") for axis in axes]
    if len(set(source)) != dimension:
        raise ValueError("axes must name different grid axes")

    kept = sorted(source)
    if dimension == 2:
        extruded = ({0, 1, 2} - set(source)).pop()
        middle = np.take(labels, labels.shape[extruded] // 2, axis=extruded)
        if not (labels == np.expand_dims(middle, extruded)).all():
            warnings.warn(
                f"grid {grid.name} is not uniform along its {'xyz'[extruded]} axis; "
                "using its middle slice"
            )
        labels = middle
    labels = np.transpose(labels, [kept.index(s) for s in source])
    shape = np.array(labels.shape)
    newLow = np.empty(dimension, np.int64)
    for d, (s, flip) in enumerate(zip(source, flipped)):
        if flip:
            labels = np.flip(labels, d)
            newLow[d] = -(low[s] + shape[d] - 1)
        else:
            newLow[d] = low[s]
    high = newLow + shape - 1

    empty = {grid.background} if ambient is None else {grid.background, ambient}
    present = [label for label in np.unique(labels).tolist() if label not in empty]
    if order is None:
        order = sorted(
            present, key=lambda label: np.nonzero(labels == label)[-1].mean()
        )
    elif sorted(order) != sorted(present):
        raise ValueError(f"order must list the labels {sorted(present)}")

    # replicating the outermost voxels keeps faces off the grid's edges
    padded = np.pad(labels, 1, mode="edge")
    core = (slice(1, -1),) * dimension
    inside = np.zeros(padded.shape, bool)
    levelSets = []
    for label in order:
        inside |= padded == label
        boundary = np.zeros(padded.shape, bool)
        for d in range(dimension):
            for step in (1, -1):
                boundary |= inside != np.roll(inside, step, d)
        points = np.argwhere(boundary[core])
        # the voxel centres next to a face lie half a voxel from it
        values = np.where(inside[core][tuple(points.T)], -0.5, 0.5)
        levelSets.append(
            (
                label,
                _levelSet(points + newLow, values, newLow, high, gridDelta, dimension),
            )
        )
    return levelSets
