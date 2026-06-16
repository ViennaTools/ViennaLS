#!/usr/bin/env python3
"""
Build and install ViennaLS with GPU-accelerated BiCGSTAB support.

The GPU solver requires CUDA 12+ and a compatible GCC version on Linux.
A virtual environment is created (or reused) and ViennaLS is installed into it.
The script auto-detects the ViennaLS source directory from its own location
(ViennaLS/python/scripts/install_ViennaLS.py → repo root is ../../).

Usage:
    python python/scripts/install_ViennaLS.py          # GPU build (default)
    python python/scripts/install_ViennaLS.py --no-gpu # CPU-only build
"""

import argparse
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

REQUIRED_NVCC_MAJOR = 12

IS_WINDOWS = sys.platform == "win32" or os.name == "nt"
IS_LINUX = sys.platform.startswith("linux")
OS_NAME = platform.system()

# Filled in by ensure_cuda() on Linux.
REQUIRED_GCC: list[str] | None = None


# ── Helpers ───────────────────────────────────────────────────────────────────

def run(cmd, **kwargs):
    print("+", " ".join(str(c) for c in cmd))
    return subprocess.run(cmd, check=True, **kwargs)


def run_capture(cmd, **kwargs):
    print("+", " ".join(str(c) for c in cmd))
    return subprocess.run(
        cmd, check=True, stdout=subprocess.PIPE, text=True, **kwargs
    ).stdout.strip()


def which_or_fail(name: str) -> str:
    p = shutil.which(name)
    if not p:
        sys.exit(f"ERROR: '{name}' is required but was not found in PATH.")
    return p


# ── CUDA / compiler checks ────────────────────────────────────────────────────

def parse_nvcc_version():
    out = run_capture(["nvcc", "--version"])
    for line in out.splitlines():
        if "release" in line:
            part = line.split("release", 1)[1].strip().split(",")[0].strip()
            try:
                parts = part.split(".")
                return int(parts[0]), int(parts[1]) if len(parts) > 1 else 0, part
            except Exception:
                break
    sys.exit("Could not parse nvcc version.")


def determine_required_gcc(nvcc_major: int, nvcc_minor: int) -> list[str] | None:
    v = nvcc_major * 10 + nvcc_minor
    if v < 124:
        return ["11", "12"]
    if v < 128:
        return ["11", "12", "13"]
    return None  # any GCC is fine


def get_default_gcc_version() -> tuple[int, int] | None:
    if not shutil.which("gcc"):
        return None
    try:
        out = run_capture(["gcc", "--version"])
        m = re.search(r"(\d+)\.(\d+)", out.splitlines()[0])
        if m:
            return int(m.group(1)), int(m.group(2))
    except Exception:
        pass
    return None


def ensure_cuda():
    global REQUIRED_GCC
    which_or_fail("nvcc")
    major, minor, full = parse_nvcc_version()
    if major < REQUIRED_NVCC_MAJOR:
        sys.exit(f"CUDA {REQUIRED_NVCC_MAJOR}.0+ required (found {full}).")
    print(f"CUDA toolkit: {full}")
    if IS_LINUX:
        REQUIRED_GCC = determine_required_gcc(major, minor)
        if REQUIRED_GCC:
            print(f"CUDA {full} requires GCC: {', '.join(REQUIRED_GCC)}")
        else:
            print(f"CUDA {full} works with any GCC version")


def ensure_compilers():
    if IS_WINDOWS:
        cl = shutil.which("cl")
        if not cl:
            print(
                "WARNING: 'cl.exe' not found. "
                "Run from a Visual Studio Developer Command Prompt."
            )
        else:
            print(f"Found MSVC: {cl}")
        return

    global REQUIRED_GCC
    if REQUIRED_GCC is None:
        which_or_fail("gcc")
        which_or_fail("g++")
        ver = get_default_gcc_version()
        print(f"Using default GCC {ver[0]}.{ver[1]}" if ver else "Using default GCC")
    else:
        found = None
        for v in REQUIRED_GCC:
            if shutil.which(f"gcc-{v}") and shutil.which(f"g++-{v}"):
                found = v
                break
        if found is None:
            default = get_default_gcc_version()
            if default and str(default[0]) in REQUIRED_GCC:
                found = str(default[0])
                print(f"Using default GCC-{found} (compatible with CUDA)")
            else:
                sys.exit(
                    f"ERROR: None of the required GCC versions {REQUIRED_GCC} found.\n"
                    "Please install one of: "
                    + ", ".join(f"gcc-{v}/g++-{v}" for v in REQUIRED_GCC)
                )
        else:
            print(f"Using GCC-{found}")
        REQUIRED_GCC = found


# ── venv ─────────────────────────────────────────────────────────────────────

def venv_paths(venv_dir: Path):
    bindir = "Scripts" if IS_WINDOWS else "bin"
    python = venv_dir / bindir / ("python.exe" if IS_WINDOWS else "python")
    pip    = venv_dir / bindir / ("pip.exe"    if IS_WINDOWS else "pip")
    return python, pip


def create_or_reuse_venv(venv_dir: Path):
    _, pip = venv_paths(venv_dir)
    if venv_dir.exists() and pip.exists():
        print(f"Reusing existing venv at {venv_dir}")
    else:
        print(f"Creating venv at {venv_dir}")
        run([sys.executable, "-m", "venv", str(venv_dir)])


# ── ViennaLS install ──────────────────────────────────────────────────────────

def get_viennals_dir() -> Path:
    # Script lives at ViennaLS/python/scripts/install_ViennaLS.py,
    # so ../.. from here is always the repo root.
    candidate = Path(__file__).resolve().parent.parent.parent
    if (candidate / "CMakeLists.txt").exists():
        return candidate
    # Fallback: current working directory if it looks like the repo root.
    cwd = Path.cwd()
    if (cwd / "CMakeLists.txt").exists():
        return cwd
    sys.exit(
        "Could not locate the ViennaLS source directory.\n"
        "Run this script from inside the ViennaLS repository."
    )


def install_viennals(pip_path: Path, viennals_dir: Path, gpu: bool,
                     debug: bool, verbose: bool):
    if not viennals_dir.exists() or not (viennals_dir / "CMakeLists.txt").exists():
        sys.exit(f"Not a valid ViennaLS directory: {viennals_dir}")

    env = os.environ.copy()
    if IS_LINUX and REQUIRED_GCC is not None:
        env["CC"]  = f"gcc-{REQUIRED_GCC}"
        env["CXX"] = f"g++-{REQUIRED_GCC}"

    cmake_args = ["-DVIENNALS_USE_GPU=ON" if gpu else "-DVIENNALS_USE_GPU=OFF"]
    env["CMAKE_ARGS"] = " ".join(cmake_args)

    cmd = [str(pip_path), "install", "--no-deps", "."]
    if debug:
        cmd += ["--config-settings=cmake.build-type=Debug"]
    if verbose:
        cmd += ["-v"]

    run(cmd, cwd=viennals_dir, env=env)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=f"Install ViennaLS with GPU BiCGSTAB support ({OS_NAME})."
    )
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Verbose pip build output.")
    parser.add_argument("--venv",
                        default=os.environ.get("VIRTUAL_ENV", ".venv"),
                        help="Path to virtual environment (default: .venv).")
    parser.add_argument("--debug-build", action="store_true",
                        help="Build in Debug mode.")
    parser.add_argument("--no-gpu", action="store_true",
                        help="Disable GPU support (CPU-only build).")
    parser.add_argument("--skip-toolchain-check", action="store_true",
                        help="Skip CUDA/compiler checks.")
    args = parser.parse_args()

    gpu = not args.no_gpu
    if gpu and not args.skip_toolchain_check:
        print("Checking toolchain...")
        ensure_cuda()
        ensure_compilers()

    venv_dir = Path(args.venv).expanduser().resolve()
    create_or_reuse_venv(venv_dir)
    _, venv_pip = venv_paths(venv_dir)

    install_viennals(venv_pip, get_viennals_dir(), gpu, args.debug_build, args.verbose)

    bindir = "Scripts" if IS_WINDOWS else "bin"
    activate = venv_dir / bindir / ("activate.bat" if IS_WINDOWS else "activate")
    print("\nInstallation complete.")
    if IS_WINDOWS:
        print(f"Activate venv:\n  {activate}")
    else:
        print(f"Activate venv:\n  source {activate}")
    print("Deactivate:\n  deactivate")


if __name__ == "__main__":
    main()
