#pragma once

#include <lsPreCompileMacros.hpp>

#include <map>

#include <hrleSparseCellIterator.hpp>
#include <hrleSparseStarIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFiniteDifferences.hpp>
#include <lsMarchingCubes.hpp>
#include <lsMesh.hpp>

namespace viennals {

using namespace viennacore;

/// Extract an explicit Mesh<> instance from an lsDomain.
/// The interface is then described by explicit surface elements:
/// Lines in 2D, Triangles in 3D.
template <class T, int D> class ToSurfaceMesh {
  typedef typename Domain<T, D>::DomainType hrleDomainType;

  SmartPointer<Domain<T, D>> levelSet = nullptr;
  SmartPointer<Mesh<T>> mesh = nullptr;
  const T epsilon;
  bool updatePointData = true;
  bool sharpCorners = true;

public:
  explicit ToSurfaceMesh(double eps = 1e-12) : epsilon(eps) {}

  ToSurfaceMesh(const SmartPointer<Domain<T, D>> passedLevelSet,
                SmartPointer<Mesh<T>> passedMesh, double eps = 1e-12)
      : levelSet(passedLevelSet), mesh(passedMesh), epsilon(eps) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void setUpdatePointData(bool update) { updatePointData = update; }

  void setSharpCorners(bool check) { sharpCorners = check; }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addError("No level set was passed to ToSurfaceMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addError("No mesh was passed to ToSurfaceMesh.")
          .print();
      return;
    }

    if (levelSet->getNumberOfPoints() == 0) {
      VIENNACORE_LOG_WARNING(
          "ToSurfaceMesh: Level set is empty. No mesh will be created.");
      return;
    }

    mesh->clear();
    const unsigned int corner0[12] = {0, 1, 2, 0, 4, 5, 6, 4, 0, 1, 3, 2};
    const unsigned int corner1[12] = {1, 3, 3, 2, 5, 7, 7, 6, 4, 5, 7, 6};
    const unsigned int direction[12] = {0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2};

    // test if level set function consists of at least 2 layers of
    // defined grid points
    if (levelSet->getLevelSetWidth() < 2) {
      VIENNACORE_LOG_WARNING("Levelset is less than 2 layers wide. Expanding "
                             "levelset to 2 layers.");
      Expand<T, D>(levelSet, 2).apply();
    }

    viennahrle::ConstSparseIterator<hrleDomainType> valueIt(levelSet->getDomain());

    typedef std::map<viennahrle::Index<D>, unsigned> nodeContainerType;

    nodeContainerType nodes[D];
    nodeContainerType faceNodes[D];
    typename nodeContainerType::iterator nodeIt;
    const bool updateData = updatePointData;

    // save how data should be transferred to new level set
    // list of indices into the old pointData vector
    std::vector<std::vector<unsigned>> newDataSourceIds;
    // there is no multithreading here, so just use 1
    if (updateData)
      newDataSourceIds.resize(1);

    // Cache for gradients to avoid re-calculation
    std::vector<Vec3D<T>> nodeNormals;
    if (sharpCorners) {
      nodeNormals.resize(levelSet->getNumberOfPoints());

      auto &domain = levelSet->getDomain();
      auto &grid = levelSet->getGrid();

      // Parallel pre-calculation of gradients
#pragma omp parallel num_threads(domain.getNumberOfSegments())
      {
        int p = 0;
#ifdef _OPENMP
        p = omp_get_thread_num();
#endif

        viennahrle::Index<D> startVector =
            (p == 0) ? grid.getMinGridPoint() : domain.getSegmentation()[p - 1];
        viennahrle::Index<D> endVector =
            (p != static_cast<int>(domain.getNumberOfSegments() - 1))
                ? domain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint());

        viennahrle::SparseStarIterator<hrleDomainType, 1> neighborIt(
            domain, startVector);

        for (; neighborIt.getIndices() < endVector; neighborIt.next()) {
          if (neighborIt.getCenter().isDefined()) {
            Vec3D<T> grad{};
            for (int i = 0; i < D; ++i) {
              bool negDefined = neighborIt.getNeighbor(i + D).isDefined();
              bool posDefined = neighborIt.getNeighbor(i).isDefined();

              if (negDefined && posDefined) {
                T valNeg = neighborIt.getNeighbor(i + D).getValue();
                T valCenter = neighborIt.getCenter().getValue();
                T valPos = neighborIt.getNeighbor(i).getValue();

                bool centerSign = valCenter > 0;
                bool negSign = valNeg > 0;
                bool posSign = valPos > 0;

                if (centerSign != negSign && centerSign == posSign) {
                  grad[i] = (valPos - valCenter);
                } else if (centerSign != posSign && centerSign == negSign) {
                  grad[i] = (valCenter - valNeg);
                } else {
                  grad[i] = 0;
                }
              } else if (negDefined) {
                grad[i] = (neighborIt.getCenter().getValue() - neighborIt.getNeighbor(i + D).getValue());
              } else if (posDefined) {
                grad[i] = (neighborIt.getNeighbor(i).getValue() - neighborIt.getCenter().getValue());
              } else {
                grad[i] = 0;
              }
            }
            Normalize(grad);
            nodeNormals[neighborIt.getCenter().getPointId()] = grad;
          }
        }
      }
    }

    // iterate over all active points
    for (viennahrle::ConstSparseCellIterator<hrleDomainType> cellIt(
             levelSet->getDomain());
         !cellIt.isFinished(); cellIt.next()) {

      for (int u = 0; u < D; u++) {
        while (!nodes[u].empty() &&
               nodes[u].begin()->first <
                   viennahrle::Index<D>(cellIt.getIndices()))
          nodes[u].erase(nodes[u].begin());

        while (!faceNodes[u].empty() &&
               faceNodes[u].begin()->first <
                   viennahrle::Index<D>(cellIt.getIndices()))
          faceNodes[u].erase(faceNodes[u].begin());
      }

      unsigned signs = 0;
      for (int i = 0; i < (1 << D); i++) {
        if (cellIt.getCorner(i).getValue() >= T(0))
          signs |= (1 << i);
      }

      // all corners have the same sign, so no surface here
      if (signs == 0)
        continue;
      if (signs == (1 << (1 << D)) - 1)
        continue;

      bool perfectCornerFound = false;
      if (sharpCorners) {
        auto getGradient = [&](int cornerID) {
          auto corner = cellIt.getCorner(cornerID);

          if (corner.isDefined()) {
            return nodeNormals[corner.getPointId()];
          }
          return Vec3D<T>{};
        };

        auto handleFaceNode = [&](int axis, int dir, Vec3D<T> pos) {
          viennahrle::Index<D> fIdx(cellIt.getIndices());
          if (dir == 1) fIdx[axis]++;

          auto it = faceNodes[axis].find(fIdx);
          if (it != faceNodes[axis].end()) {
            return it->second;
          }
          unsigned nid = mesh->insertNextNode(pos);
          if (updateData)
            newDataSourceIds[0].push_back(newDataSourceIds[0].back());
          faceNodes[axis][fIdx] = nid;

          // Stitching logic
          viennahrle::Index<D> neighborIdx = cellIt.getIndices();
          if (dir == 1) {
            neighborIdx[axis]++;
          } else {
            neighborIdx[axis]--;
          }

          if (neighborIdx < cellIt.getIndices()) {
            // Neighbor has already been processed. Assume it was a linear cell.
            // We need to find the edge it should have created on the shared face.
            unsigned int neighborSigns = 0;
            for (int i = 0; i < (1 << D); i++) {
              valueIt.goToIndices(neighborIdx + viennahrle::BitMaskToIndex<D>(i));
              if (valueIt.isDefined() && valueIt.getValue() >= T(0))
                neighborSigns |= (1 << i);
            }

            if (neighborSigns > 0 && neighborSigns < (1 << (1 << D)) - 1) {
              const int *Triangles = lsInternal::MarchingCubes::polygonize3d(neighborSigns);
              for (; Triangles[0] != -1; Triangles += 3) {
                  std::vector<unsigned> face_edge_nodes;
                  int num_on_face = 0;
                  for (int i = 0; i < 3; ++i) {
                      int edge = Triangles[i];
                      int c0 = corner0[edge];
                      int c1 = corner1[edge];
                      
                      // An edge is on the shared face if both its points are on that face plane.
                      // The shared face is at the boundary of the neighbor cell.
                      // If dir is 1 (current cell is at i), neighbor is at i+1. Shared face is at i+1.
                      // If dir is 0 (current cell is at i), neighbor is at i-1. Shared face is at i.
                      int face_bit = (dir == 1) ? 0 : 1;
                      
                      if ((((c0 >> axis) & 1) == face_bit) && (((c1 >> axis) & 1) == face_bit)) {
                        num_on_face++;
                      }
                  }

                  if (num_on_face == 2) {
                    // This triangle has an edge on the shared face.
                    // Find the two nodes of that edge.
                    for (int i=0; i<3; ++i) {
                      int edge = Triangles[i];
                      int c0 = corner0[edge];
                      int c1 = corner1[edge];
                      int face_bit = (dir == 1) ? 0 : 1;
                      if ((((c0 >> axis) & 1) == face_bit) && (((c1 >> axis) & 1) == face_bit)) {
                        
                        auto n_dir = direction[edge];
                        viennahrle::Index<D> d(neighborIdx);
                        d += viennahrle::BitMaskToIndex<D>(c0);

                        auto nodeIt = nodes[n_dir].find(d);
                        if (nodeIt != nodes[n_dir].end()) {
                          face_edge_nodes.push_back(nodeIt->second);
                        } else {
                          // Node should have been created by neighbor, but let's be safe.
                          // This part is identical to getNode, but for neighborIdx
                          Vec3D<T> cc{};
                          std::size_t currentPointId = 0;
                          for (int z = 0; z < D; z++) {
                            if (z != n_dir) {
                              cc[z] = static_cast<double>(neighborIdx[z] + viennahrle::BitMaskToIndex<D>(c0)[z]);
                            } else {
                              valueIt.goToIndices(neighborIdx + viennahrle::BitMaskToIndex<D>(c0));
                              T d0 = valueIt.getValue();
                              valueIt.goToIndices(neighborIdx + viennahrle::BitMaskToIndex<D>(c1));
                              T d1 = valueIt.getValue();
                              if (d0 == -d1) {
                                cc[z] = static_cast<T>(neighborIdx[z]) + 0.5;
                              } else {
                                if (std::abs(d0) <= std::abs(d1)) {
                                  cc[z] = static_cast<T>(neighborIdx[z]) + (d0 / (d0 - d1));
                                } else {
                                  cc[z] = static_cast<T>(neighborIdx[z] + 1) - (d1 / (d1 - d0));
                                }
                              }
                              cc[z] = std::max(cc[z], neighborIdx[z] + epsilon);
                              cc[z] = std::min(cc[z], (neighborIdx[z] + 1) - epsilon);
                            }
                            cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
                          }
                          unsigned newNodeId = mesh->insertNextNode(cc);
                          nodes[n_dir][d] = newNodeId;
                          if (updateData) newDataSourceIds[0].push_back(0); // No good source info here
                          face_edge_nodes.push_back(newNodeId);
                        }
                      }
                    }

                    if(face_edge_nodes.size() == 2) {
                      mesh->insertNextElement(std::array<unsigned, 3>{face_edge_nodes[0], face_edge_nodes[1], nid});
                    }
                  }
              }
            }
          }
          return nid;
        };

        // Check for perfect corner (2D only)
        if constexpr (D == 2) {
          int countNeg = 0;
          int countPos = 0;
          int negIdx = -1;
          int posIdx = -1;

          for (int i = 0; i < 4; ++i) {
            T val = cellIt.getCorner(i).getValue();
            if (val < -epsilon) {
              countNeg++;
              negIdx = i;
            } else if (val > epsilon) {
              countPos++;
              posIdx = i;
            }
          }

          int cornerIdx = -1;
          int insideCount = 0;

          if (countNeg == 1) {
            cornerIdx = negIdx;
            insideCount = 1;
          } else if (countPos == 1) {
            cornerIdx = posIdx;
            insideCount = 3;
          }

          if (cornerIdx != -1) {
            // Check if this corner is also a corner in any neighboring cell
            bool isSharedCorner = false;
            viennahrle::Index<D> pIdx = cellIt.getIndices() + viennahrle::BitMaskToIndex<D>(cornerIdx);
            T pVal = cellIt.getCorner(cornerIdx).getValue();
            auto &grid = levelSet->getGrid();

            for(int i=0; i<(1<<D); ++i) {
                if (i == cornerIdx) continue;

                viennahrle::Index<D> neighborIndices;
                for(int k=0; k<D; ++k) neighborIndices[k] = pIdx[k] - ((i >> k) & 1);

                bool neighborIsCorner = true;

                // Check edge-connected neighbors in the neighbor cell
                for(int k=0; k<D; ++k) {
                    int neighborLocal = i ^ (1 << k);
                    viennahrle::Index<D> checkIdx = neighborIndices + viennahrle::BitMaskToIndex<D>(neighborLocal);
                    if (grid.isOutsideOfDomain(checkIdx)) {
                        checkIdx = grid.globalIndices2LocalIndices(checkIdx);
                    }
                    valueIt.goToIndices(checkIdx);
                    T nVal = valueIt.getValue();
                    if (std::abs(nVal) > epsilon && ((pVal >= 0) == (nVal > 0))) {
                        neighborIsCorner = false;
                        break;
                    }
                }

                if (neighborIsCorner) {
                    isSharedCorner = true;
                    break;
                }
            }

            if (!isSharedCorner) {
            int n1 = cornerIdx ^ 1;
            int n2 = cornerIdx ^ 2;

            Vec3D<T> norm1 = getGradient(n1);
            Vec3D<T> norm2 = getGradient(n2);

            if (std::abs(DotProduct(norm1, norm2)) < 0.9) {
              int edgeX = -1, edgeY = -1;
              if (cornerIdx == 0) {
                edgeX = 0;
                edgeY = 3;
              } else if (cornerIdx == 1) {
                edgeX = 0;
                edgeY = 1;
              } else if (cornerIdx == 2) {
                edgeX = 2;
                edgeY = 3;
              } else if (cornerIdx == 3) {
                edgeX = 2;
                edgeY = 1;
              }

              // Helper to get/create node
              auto getNode = [&](int edge) {
                unsigned p0 = corner0[edge];
                unsigned p1 = corner1[edge];
                auto dir = direction[edge];
                viennahrle::Index<D> d(cellIt.getIndices());
                d += viennahrle::BitMaskToIndex<D>(p0);

                if (nodes[dir].find(d) != nodes[dir].end()) {
                  return nodes[dir][d];
                }

                // Create node
                Vec3D<T> cc{};
                std::size_t currentPointId = 0;
                for (int z = 0; z < D; z++) {
                  if (z != dir) {
                    cc[z] = static_cast<double>(
                        cellIt.getIndices(z) +
                        viennahrle::BitMaskToIndex<D>(p0)[z]);
                  } else {
                    T d0 = cellIt.getCorner(p0).getValue();
                    T d1 = cellIt.getCorner(p1).getValue();
                    if (d0 == -d1) {
                      currentPointId = cellIt.getCorner(p0).getPointId();
                      cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                    } else {
                      if (std::abs(d0) <= std::abs(d1)) {
                        currentPointId = cellIt.getCorner(p0).getPointId();
                        cc[z] = static_cast<T>(cellIt.getIndices(z)) +
                                (d0 / (d0 - d1));
                      } else {
                        currentPointId = cellIt.getCorner(p1).getPointId();
                        cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) -
                                (d1 / (d1 - d0));
                      }
                    }
                    cc[z] = std::max(cc[z], cellIt.getIndices(z) + epsilon);
                    cc[z] =
                        std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
                  }
                  cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
                }
                unsigned nodeId = mesh->insertNextNode(cc);
                nodes[dir][d] = nodeId;
                if (updateData)
                  newDataSourceIds[0].push_back(currentPointId);
                return nodeId;
              };

              unsigned nX = getNode(edgeX);
              unsigned nY = getNode(edgeY);

              auto pX = mesh->getNodes()[nX];
              auto pY = mesh->getNodes()[nY];

              double d1 = DotProduct(norm1, pX);
              double d2 = DotProduct(norm2, pY);
              double det = norm1[0] * norm2[1] - norm1[1] * norm2[0];

              Vec3D<T> cornerPos{};
              if (std::abs(det) > 1e-6) {
                cornerPos[0] = (d1 * norm2[1] - d2 * norm1[1]) / det;
                cornerPos[1] = (d2 * norm1[0] - d1 * norm2[0]) / det;
              } else {
                for (int i = 0; i < D; ++i)
                  cornerPos[i] = (pX[i] + pY[i]) * 0.5;
              }

              perfectCornerFound = true;
              unsigned nCorner = mesh->insertNextNode(cornerPos);
              if (updateData)
                newDataSourceIds[0].push_back(newDataSourceIds[0].back());

              // Add lines
              mesh->insertNextElement(std::array<unsigned, 2>{nX, nCorner});
              mesh->insertNextElement(std::array<unsigned, 2>{nCorner, nY});
            }
          }
          }
        } else if constexpr (D == 3) {
          // Helper to get/create node on edge
          auto getNode = [&](int edge) {
            unsigned p0 = corner0[edge];
            unsigned p1 = corner1[edge];
            auto dir = direction[edge];
            viennahrle::Index<D> d(cellIt.getIndices());
            d += viennahrle::BitMaskToIndex<D>(p0);

            if (nodes[dir].find(d) != nodes[dir].end()) {
              return nodes[dir][d];
            }

            // Create node
            Vec3D<T> cc{};
            std::size_t currentPointId = 0;
            for (int z = 0; z < D; z++) {
              if (z != dir) {
                cc[z] = static_cast<double>(
                    cellIt.getIndices(z) +
                    viennahrle::BitMaskToIndex<D>(p0)[z]);
              } else {
                T d0 = cellIt.getCorner(p0).getValue();
                T d1 = cellIt.getCorner(p1).getValue();
                if (d0 == -d1) {
                  currentPointId = cellIt.getCorner(p0).getPointId();
                  cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                } else {
                  if (std::abs(d0) <= std::abs(d1)) {
                    currentPointId = cellIt.getCorner(p0).getPointId();
                    cc[z] = static_cast<T>(cellIt.getIndices(z)) +
                            (d0 / (d0 - d1));
                  } else {
                    currentPointId = cellIt.getCorner(p1).getPointId();
                    cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) -
                            (d1 / (d1 - d0));
                  }
                }
                cc[z] = std::max(cc[z], cellIt.getIndices(z) + epsilon);
                cc[z] =
                    std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
              }
              cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
            }
            unsigned nodeId = mesh->insertNextNode(cc);
            nodes[dir][d] = nodeId;
            if (updateData)
              newDataSourceIds[0].push_back(currentPointId);
            return nodeId;
          };

          int insideCount = 0;
          for (int i = 0; i < 8; ++i) {
            if (signs & (1 << i))
              insideCount++;
          }

          // Check for a corner
          int cornerIdx = -1;
          if (insideCount == 1) {
            for (int i = 0; i < 8; ++i)
              if (signs & (1 << i)) {
                cornerIdx = i;
                break;
              }
          } else if (insideCount == 7) {
            for (int i = 0; i < 8; ++i)
              if (!(signs & (1 << i))) {
                cornerIdx = i;
                break;
              }
          }

          if (cornerIdx != -1) {
            // Check if this corner is also a corner in any neighboring cell
            bool isSharedCorner = false;
            viennahrle::Index<D> pIdx = cellIt.getIndices() + viennahrle::BitMaskToIndex<D>(cornerIdx);
            T pVal = cellIt.getCorner(cornerIdx).getValue();
            auto &grid = levelSet->getGrid();

            for(int i=0; i<(1<<D); ++i) {
                if (i == cornerIdx) continue;

                viennahrle::Index<D> neighborIndices;
                for(int k=0; k<D; ++k) neighborIndices[k] = pIdx[k] - ((i >> k) & 1);

                bool neighborIsCorner = true;

                // Check edge-connected neighbors in the neighbor cell
                for(int k=0; k<D; ++k) {
                    int neighborLocal = i ^ (1 << k);
                    viennahrle::Index<D> checkIdx = neighborIndices + viennahrle::BitMaskToIndex<D>(neighborLocal);
                    if (grid.isOutsideOfDomain(checkIdx)) {
                        checkIdx = grid.globalIndices2LocalIndices(checkIdx);
                    }
                    valueIt.goToIndices(checkIdx);
                    T nVal = valueIt.getValue();
                    if (std::abs(nVal) > epsilon && ((pVal >= 0) == (nVal > 0))) {
                        neighborIsCorner = false;
                        break;
                    }
                }

                if (neighborIsCorner) {
                    isSharedCorner = true;
                    break;
                }
            }

            if (!isSharedCorner) {
            // cornerIdx is the single point in- or outside the geometry
            // The three edges connected to this point are the basis for our calculation
            int edgeIdx[3];
            int c = 0;
            for (int i = 0; i < 12; i++) {
              if (corner0[i] == cornerIdx || corner1[i] == cornerIdx) {
                edgeIdx[c++] = i;
              }
            }

            // The three neighbor points to the corner point
            int n[3];
            c = 0;
            for (int i = 0; i < 3; i++) {
              if (corner0[edgeIdx[i]] == cornerIdx) {
                n[i] = corner1[edgeIdx[i]];
              } else {
                n[i] = corner0[edgeIdx[i]];
              }
            }
            
            Vec3D<T> norm1 = getGradient(n[0]);
            Vec3D<T> norm2 = getGradient(n[1]);
            Vec3D<T> norm3 = getGradient(n[2]);

            // The three points on the edges generated by marching cubes
            unsigned nId1 = getNode(edgeIdx[0]);
            unsigned nId2 = getNode(edgeIdx[1]);
            unsigned nId3 = getNode(edgeIdx[2]);

            auto p1 = mesh->getNodes()[nId1];
            auto p2 = mesh->getNodes()[nId2];
            auto p3 = mesh->getNodes()[nId3];

            // Now we define 3 planes from these points and the (one-sided) gradients
            T d1 = DotProduct(norm1, p1);
            T d2 = DotProduct(norm2, p2);
            T d3 = DotProduct(norm3, p3);

            // And find the intersection of these three planes
            Vec3D<T> cross23 = CrossProduct(norm2, norm3);
            T det = DotProduct(norm1, cross23);

            if (std::abs(det) > 1e-6) {
              Vec3D<T> cornerPos = (d1 * cross23 + d2 * CrossProduct(norm3, norm1) + d3 * CrossProduct(norm1, norm2)) / det;
              
              // Fallback check
              Vec3D<T> centroid = (p1 + p2 + p3) / static_cast<T>(3.0);
              Vec3D<T> avgNorm = (norm1 + norm2 + norm3);
              Normalize(avgNorm);
              T dot = DotProduct(cornerPos - centroid, avgNorm);

              if ((insideCount == 1 && dot < 0) || (insideCount == 7 && dot > 0)) {
                perfectCornerFound = true;

                // Also find face intersection points
                unsigned int axis[3] = {direction[edgeIdx[0]], direction[edgeIdx[1]], direction[edgeIdx[2]]};
                
                // Face between norm1 and norm2 is perpendicular to axis3
                Vec3D<T> faceNorm3 = {0,0,0};
                faceNorm3[axis[2]] = 1.0;
                T faceD3 = (cellIt.getIndices(axis[2]) + ((cornerIdx >> axis[2]) & 1)) * levelSet->getGrid().getGridDelta();
                Vec3D<T> cross12 = CrossProduct(norm1, norm2);
                T det12 = DotProduct(faceNorm3, cross12);
                Vec3D<T> facePos3;
                if (std::abs(det12) > 1e-6)
                  facePos3 = (d1*CrossProduct(norm2, faceNorm3) + d2*CrossProduct(faceNorm3, norm1) + faceD3*cross12) / det12;
                else
                  facePos3 = (p1+p2)/static_cast<T>(2.0);

                // Face between norm2 and norm3 is perpendicular to axis1
                Vec3D<T> faceNorm1 = {0,0,0};
                faceNorm1[axis[0]] = 1.0;
                T faceD1 = (cellIt.getIndices(axis[0]) + ((cornerIdx >> axis[0]) & 1)) * levelSet->getGrid().getGridDelta();
                T det23 = DotProduct(faceNorm1, cross23);
                Vec3D<T> facePos1;
                if (std::abs(det23) > 1e-6)
                  facePos1 = (d2*CrossProduct(norm3, faceNorm1) + d3*CrossProduct(faceNorm1, norm2) + faceD1*cross23) / det23;
                else
                  facePos1 = (p2+p3)/static_cast<T>(2.0);

                // Face between norm3 and norm1 is perpendicular to axis2
                Vec3D<T> faceNorm2 = {0,0,0};
                faceNorm2[axis[1]] = 1.0;
                T faceD2 = (cellIt.getIndices(axis[1]) + ((cornerIdx >> axis[1]) & 1)) * levelSet->getGrid().getGridDelta();
                Vec3D<T> cross31 = CrossProduct(norm3, norm1);
                T det31 = DotProduct(faceNorm2, cross31);
                Vec3D<T> facePos2;
                if (std::abs(det31) > 1e-6)
                  facePos2 = (d3*CrossProduct(norm1, faceNorm2) + d1*CrossProduct(faceNorm2, norm3) + faceD2*cross31) / det31;
                else
                  facePos2 = (p3+p1)/static_cast<T>(2.0);
                

                unsigned nCorner = mesh->insertNextNode(cornerPos);
                if (updateData) newDataSourceIds[0].push_back(newDataSourceIds[0].back());

                unsigned nFace1 = handleFaceNode(axis[0], (cornerIdx >> axis[0]) & 1, facePos1);
                unsigned nFace2 = handleFaceNode(axis[1], (cornerIdx >> axis[1]) & 1, facePos2);
                unsigned nFace3 = handleFaceNode(axis[2], (cornerIdx >> axis[2]) & 1, facePos3);

                auto addTri = [&](unsigned a, unsigned b, unsigned c) {
                  if (insideCount == 7) mesh->insertNextElement(std::array<unsigned,3>{a,c,b});
                  else mesh->insertNextElement(std::array<unsigned,3>{a,b,c});
                };

                // Triangulate the 3 quads meeting at nCorner
                addTri(nCorner, nId1, nFace2);
                addTri(nCorner, nFace2, nId3);

                addTri(nCorner, nId2, nFace3);
                addTri(nCorner, nFace3, nId1);

                addTri(nCorner, nId3, nFace1);
                addTri(nCorner, nFace1, nId2);
              }
            }
          }

          // if no corner was found, check for an edge
          if (!perfectCornerFound) {
            int edgeAxis = -1;
            int c1 = -1, c2 = -1;
            if (insideCount == 2) {
              for (int i = 0; i < 8; ++i) {
                if (signs & (1 << i)) {
                  if (c1 == -1) c1 = i;
                  else c2 = i;
                }
              }
            } else if (insideCount == 6) {
              for (int i = 0; i < 8; ++i) {
                if (!(signs & (1 << i))) {
                  if (c1 == -1) c1 = i;
                  else c2 = i;
                }
              }
            }

            int diff = c1 ^ c2;
            if (diff == 1) edgeAxis = 0;
            else if (diff == 2) edgeAxis = 1;
            else if (diff == 4) edgeAxis = 2;

            if (edgeAxis != -1) {
              // We have an edge along edgeAxis
              // It can be treated as two 2D corners on the slices perpendicular to the edge
              
              auto checkSlice = [&](int sliceIdx) -> std::pair<bool, std::array<unsigned, 3>> {
                // Returns {success, {sharpNode, mcNode1, mcNode2}}
                
                // find corner on this slice
                int corner = -1;
                int numInside = 0;
                for (int i=0; i < 8; ++i) {
                  if (((i >> edgeAxis) & 1) == sliceIdx) {
                    bool isInside = (signs & (1 << i)) > 0;
                    if (isInside) numInside++;
                  }
                }

                int insideOnSlice = 0;
                if (insideCount == 2) insideOnSlice = 1;
                else if (insideCount == 6) insideOnSlice = 3;
                
                if (numInside == insideOnSlice) {
                  for (int i = 0; i < 8; ++i) {
                    if (((i >> edgeAxis) & 1) == sliceIdx) {
                      bool isInside = (signs & (1 << i)) > 0;
                      if ((insideCount == 2 && isInside) || (insideCount == 6 && !isInside)) {
                        corner = i;
                        break;
                      }
                    }
                  }
                }
                
                if (corner == -1) return {false, {}};

                int axes[2];
                int idx=0;
                for(int i=0; i<3; ++i) if (i != edgeAxis) axes[idx++] = i;

                int n1_idx = corner ^ (1 << axes[0]);
                int n2_idx = corner ^ (1 << axes[1]);
                
                Vec3D<T> g1 = getGradient(n1_idx);
                Vec3D<T> g2 = getGradient(n2_idx);

                // find edge indices
                int edge1_idx = -1, edge2_idx = -1;
                for (int i=0; i<12; ++i) {
                  if ((corner0[i] == corner && corner1[i] == n1_idx) || (corner1[i] == corner && corner0[i] == n1_idx)) edge1_idx = i;
                  if ((corner0[i] == corner && corner1[i] == n2_idx) || (corner1[i] == corner && corner0[i] == n2_idx)) edge2_idx = i;
                }

                unsigned mc_n1 = getNode(edge1_idx);
                unsigned mc_n2 = getNode(edge2_idx);

                auto p1 = mesh->getNodes()[mc_n1];
                auto p2 = mesh->getNodes()[mc_n2];

                T d1 = DotProduct(g1, p1);
                T d2 = DotProduct(g2, p2);
                
                T det = g1[axes[0]]*g2[axes[1]] - g1[axes[1]]*g2[axes[0]];

                if (std::abs(det) > 1e-6) {
                  Vec3D<T> sharpPos;
                  sharpPos[edgeAxis] = levelSet->getGrid().getGridDelta() * (cellIt.getIndices(edgeAxis) + sliceIdx);
                  sharpPos[axes[0]] = (d1 * g2[axes[1]] - d2*g1[axes[1]])/det;
                  sharpPos[axes[1]] = (d2 * g1[axes[0]] - d1*g2[axes[0]])/det;
                  
                  // Fallback check
                  Vec3D<T> midPoint = (p1+p2)/static_cast<T>(2.0);
                  Vec3D<T> avgNorm = g1+g2;
                  T dot = DotProduct(sharpPos-midPoint, avgNorm);

                  bool sliceIsConvex = (insideOnSlice == 1);

                  if ((sliceIsConvex && dot < 0) || (!sliceIsConvex && dot > 0)) {
                    unsigned sharpNode = handleFaceNode(edgeAxis, sliceIdx, sharpPos);
                    if (updateData) newDataSourceIds[0].push_back(newDataSourceIds[0].back());
                    return {true, {sharpNode, mc_n1, mc_n2}};
                  }
                }
                return {false, {0, mc_n1, mc_n2}};
              };

              auto res0 = checkSlice(0);
              auto res1 = checkSlice(1);
              
              if (res0.first && res1.first) {
                perfectCornerFound = true;
                unsigned s0 = res0.second[0];
                unsigned n00 = res0.second[1];
                unsigned n01 = res0.second[2];

                unsigned s1 = res1.second[0];
                unsigned n10 = res1.second[1];
                unsigned n11 = res1.second[2];

                // Check winding
                Vec3D<T> pA = mesh->getNodes()[s0];
                Vec3D<T> pB = mesh->getNodes()[n00];
                Vec3D<T> pC = mesh->getNodes()[n10];
                Vec3D<T> center;
                for (int k=0; k<3; ++k) center[k] = (cellIt.getIndices(k) + 0.5) * levelSet->getGrid().getGridDelta();
                Vec3D<T> n = CrossProduct(pB-pA, pC-pA);
                Vec3D<T> dir = pA - center;
                bool reversed = DotProduct(n, dir) > 0;
                
                if ( (insideCount == 2 && !reversed) || (insideCount == 6 && reversed)) {
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n00, s1});
                  mesh->insertNextElement(std::array<unsigned, 3>{n00, n10, s1});
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, s1, n01});
                  mesh->insertNextElement(std::array<unsigned, 3>{s1, n11, n01});
                } else {
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, s1, n00});
                  mesh->insertNextElement(std::array<unsigned, 3>{n00, s1, n10});
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n01, s1});
                  mesh->insertNextElement(std::array<unsigned, 3>{s1, n01, n11});
                }
              } else if (res0.first && !res1.first) {
                perfectCornerFound = true;
                unsigned s0 = res0.second[0];
                unsigned n00 = res0.second[1];
                unsigned n01 = res0.second[2];
                unsigned n10 = res1.second[1];
                unsigned n11 = res1.second[2];

                // Triangulate {s0, n00, n01} with {n10, n11}
                // Fallback to just using the marching cubes points on slice 1
                // This is a quad {n00, n01, n11, n10} and a point s0
                // We make a fan triangulation from s0
                if (insideCount == 2) {
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n00, n10});
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n10, n11});
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n11, n01});
                } else {
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n10, n00});
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n11, n10});
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n01, n11});
                }
              } else if (!res0.first && res1.first) {
                perfectCornerFound = true;
                unsigned s1 = res1.second[0];
                unsigned n10 = res1.second[1];
                unsigned n11 = res1.second[2];
                unsigned n00 = res0.second[1];
                unsigned n01 = res0.second[2];
                
                if (insideCount == 2) {
                  mesh->insertNextElement(std::array<unsigned, 3>{s1, n10, n00});
                  mesh->insertNextElement(std::array<unsigned, 3>{s1, n00, n01});
                  mesh->insertNextElement(std::array<unsigned, 3>{s1, n01, n11});
                } else {
                  mesh->insertNextElement(std::array<unsigned, 3>{s1, n00, n10});
                  mesh->insertNextElement(std::array<unsigned, 3>{s1, n01, n00});
                  mesh->insertNextElement(std::array<unsigned, 3>{s1, n11, n01});
                }
              }
            }
          }
          }
        }
      }

      if (perfectCornerFound)
        continue;


      if constexpr (D == 3) {
        // Stitch to perfect corners/edges
        for (int axis = 0; axis < 3; ++axis) {
          for (int d = 0; d < 2; ++d) {
            viennahrle::Index<D> faceKey(cellIt.getIndices());
            if (d == 1)
              faceKey[axis]++;

            auto it = faceNodes[axis].find(faceKey);
            if (it != faceNodes[axis].end()) {
              const unsigned faceNodeId = it->second;
              const int *Triangles = lsInternal::MarchingCubes::polygonize3d(signs);

              auto getNodeOnEdge = [&](int edge) -> unsigned {
                unsigned p0 = corner0[edge];
                unsigned p1 = corner1[edge];
                auto dir = direction[edge];
                viennahrle::Index<D> nodeKey(cellIt.getIndices());
                nodeKey += viennahrle::BitMaskToIndex<D>(p0);

                auto nodeIt = nodes[dir].find(nodeKey);
                if (nodeIt != nodes[dir].end()) {
                  return nodeIt->second;
                }

                // if node does not exist yet, create it
                Vec3D<T> cc{}; // initialise with zeros
                std::size_t currentPointId = 0;
                for (int z = 0; z < D; z++) {
                  if (z != dir) {
                    cc[z] = static_cast<double>(
                        cellIt.getIndices(z) +
                        viennahrle::BitMaskToIndex<D>(p0)[z]);
                  } else {
                    T d0 = cellIt.getCorner(p0).getValue();
                    T d1 = cellIt.getCorner(p1).getValue();
                    if (d0 == -d1) {
                      currentPointId = cellIt.getCorner(p0).getPointId();
                      cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                    } else {
                      if (std::abs(d0) <= std::abs(d1)) {
                        currentPointId = cellIt.getCorner(p0).getPointId();
                        cc[z] = static_cast<T>(cellIt.getIndices(z)) +
                                (d0 / (d0 - d1));
                      } else {
                        currentPointId = cellIt.getCorner(p1).getPointId();
                        cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) -
                                (d1 / (d1 - d0));
                      }
                    }
                    cc[z] = std::max(cc[z], cellIt.getIndices(z) + epsilon);
                    cc[z] =
                        std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
                  }
                  cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
                }

                unsigned nodeId = mesh->insertNextNode(cc);
                nodes[dir][nodeKey] = nodeId;
                if (updateData)
                  newDataSourceIds[0].push_back(currentPointId);

                return nodeId;
              };

              for (; Triangles[0] != -1; Triangles += 3) {
                std::vector<unsigned> face_edge_nodes;
                for (int i = 0; i < 3; ++i) {
                  int edge = Triangles[i];
                  int c0 = corner0[edge];
                  int c1 = corner1[edge];
                  bool onFace = (((c0 >> axis) & 1) == d) && (((c1 >> axis) & 1) == d);
                  if (onFace) {
                    face_edge_nodes.push_back(getNodeOnEdge(edge));
                  }
                }
                if (face_edge_nodes.size() == 2) {
                  mesh->insertNextElement(std::array<unsigned, 3>{
                      face_edge_nodes[0], face_edge_nodes[1], faceNodeId});
                }
              }
            }
          }
        }
      }

      // for each element
      const int *Triangles =
          (D == 2) ? lsInternal::MarchingCubes::polygonize2d(signs)
                   : lsInternal::MarchingCubes::polygonize3d(signs);

      for (; Triangles[0] != -1; Triangles += D) {
        std::array<unsigned, D> nod_numbers;

        // for each node
        for (int n = 0; n < D; n++) {
          const int edge = Triangles[n];

          unsigned p0 = corner0[edge];
          unsigned p1 = corner1[edge];

          // determine direction of edge
          auto dir = direction[edge];

          // look for existing surface node
          viennahrle::Index<D> d(cellIt.getIndices());
          d += viennahrle::BitMaskToIndex<D>(p0);

          nodeIt = nodes[dir].find(d);
          if (nodeIt != nodes[dir].end()) {
            nod_numbers[n] = nodeIt->second;
          } else { // if node does not exist yet

            // calculate coordinate of new node
            Vec3D<T> cc{}; // initialise with zeros
            std::size_t currentPointId = 0;
            for (int z = 0; z < D; z++) {
              if (z != dir) {
                // TODO might not need BitMaskToVector here, just check if z bit
                // is set
                cc[z] =
                    static_cast<double>(cellIt.getIndices(z) +
                                        viennahrle::BitMaskToIndex<D>(p0)[z]);
              } else {
                T d0, d1;

                d0 = cellIt.getCorner(p0).getValue();
                d1 = cellIt.getCorner(p1).getValue();

                // calculate the surface-grid intersection point
                if (d0 == -d1) { // includes case where d0=d1=0
                  currentPointId = cellIt.getCorner(p0).getPointId();
                  cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                } else {
                  if (std::abs(d0) <= std::abs(d1)) {
                    currentPointId = cellIt.getCorner(p0).getPointId();
                    cc[z] =
                        static_cast<T>(cellIt.getIndices(z)) + (d0 / (d0 - d1));
                  } else {
                    currentPointId = cellIt.getCorner(p1).getPointId();
                    cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) -
                            (d1 / (d1 - d0));
                  }
                }
                cc[z] = std::max(cc[z], cellIt.getIndices(z) + epsilon);
                cc[z] = std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
              }
              cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
            }

            // insert new node
            nod_numbers[n] =
                mesh->insertNextNode(cc); // insert new surface node
            nodes[dir][d] = nod_numbers[n];

            if (updateData)
              newDataSourceIds[0].push_back(currentPointId);
          }
        }

        if (!triangleMisformed(nod_numbers))
          mesh->insertNextElement(nod_numbers); // insert new surface element
      }
    }

    // now copy old data into new level set
    if (updateData) {
      mesh->getPointData().translateFromMultiData(levelSet->getPointData(),
                                                  newDataSourceIds);
    }
  }

private:
  static bool triangleMisformed(const std::array<unsigned, D> &nodeNumbers) {
    if constexpr (D == 3) {
      return nodeNumbers[0] == nodeNumbers[1] ||
             nodeNumbers[0] == nodeNumbers[2] ||
             nodeNumbers[1] == nodeNumbers[2];
    } else {
      return nodeNumbers[0] == nodeNumbers[1];
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(ToSurfaceMesh)

} // namespace viennals
