#ifndef LS_VTK_READER_HPP
#define LS_VTK_READER_HPP

#include <fstream>
#include <string>

#include <lsMesh.hpp>
#include <lsMessage.hpp>

class lsVTKReader {
  lsMesh &mesh;

  unsigned vtk_nodes_for_cell_type[15] = {0, 1, 0, 2, 0, 3, 0, 0,
                                          4, 4, 4, 8, 8, 6, 5};

public:
  lsVTKReader(lsMesh &passedMesh) : mesh(passedMesh) {}

  void readVTKLegacy(std::string filename) {
    // open geometry file
    std::ifstream f(filename.c_str());
    if (!f)
      lsMessage::getInstance().addError("Could not open geometry file!");
    std::string temp;

    // Check if geometry is an unstructured grid as required
    while (std::getline(f, temp)) {
      if (temp.find("DATASET") != std::string::npos)
        break;
    }
    if (temp.find("UNSTRUCTURED_GRID") == std::string::npos) {
      lsMessage::getInstance().addError("DATASET is not an UNSTRUCTURED_GRID!");
    }

    // Find POINTS in file to know number of nodes to read in
    while (std::getline(f, temp)) {
      if (temp.find("POINTS") != std::string::npos)
        break;
    }
    int num_nodes = atoi(&temp[temp.find(" ") + 1]);

    mesh.nodes.resize(num_nodes);

    for (int i = 0; i < num_nodes; i++) {
      double coords[3];

      for (int j = 0; j < 3; j++)
        f >> coords[j];
      for (int j = 0; j < 3; j++) {
        mesh.nodes[i][j] = coords[j];
        // int shift_size = shift.size();
        // if (shift_size > j)
        //   mesh.nodes[i][j] += shift[j]; // Assign desired shift
        // if (InputTransformationSigns[j])
        //   mesh.nodes[i][j] = -mesh.nodes[i][j]; // Assign sign
        //   transformation, if needed
        // mesh.nodes[i][j] *= scale; // Scale the geometry according to
        // parameters file
      }
    }

    while (std::getline(f, temp)) {
      if (temp.find("CELLS") == 0)
        break;
    }

    int num_elems = atoi(&temp[temp.find(" ") + 1]);

    std::ifstream f_ct(filename.c_str()); // stream to read cell CELL_TYPES
    std::ifstream f_m(
        filename.c_str()); // stream for material numbers if they exist

    // advance to cell types and check if there are the right number
    while (std::getline(f_ct, temp)) {
      if (temp.find("CELL_TYPES") == 0)
        break;
    }
    int num_cell_types = atoi(&temp[temp.find(" ") + 1]);
    // need a cell_type for each cell
    if (num_elems != num_cell_types) {
      lsMessage::getInstance().addError(
          "Corrupt input geometry! Number of CELLS and CELL_TYPES "
          "is different!");
    }

    bool is_material = true;
    // advance to material if it is specified
    while (std::getline(f_m, temp)) {
      if (temp.find("CELL_DATA") != std::string::npos) {
        std::getline(f_m, temp);
        if ((temp.find("SCALARS material") != std::string::npos) ||
            (temp.find("SCALARS Material") != std::string::npos)) {
          std::getline(f_m, temp);
          break;
        }
      }
    }
    if (f_m.eof()) {
      is_material = false;
    }

    // maximum cell to read in is a tetra with 4 points
    std::vector<hrleVectorType<unsigned int, 4>> elements;
    elements.reserve(num_elems);

    auto materials = mesh.getMaterials();
    materials.clear();

    unsigned elems_fake;
    unsigned cell_type;
    unsigned cell_material;
    for (int i = 0; i < num_elems; i++) {
      f >> elems_fake;
      f_ct >> cell_type;
      if (is_material)
        f_m >> cell_material;
      else
        cell_material =
            1; // if there are no materials specified make all the same

      // check if the correct number of nodes for cell_type is given
      unsigned number_nodes = vtk_nodes_for_cell_type[cell_type];
      if (number_nodes == elems_fake || number_nodes == 0) {
        // check for different types to subdivide them into supported types
        switch (cell_type) {
        case 3: {
          hrleVectorType<unsigned, 2> elem;
          for (unsigned j = 0; j < number_nodes; ++j) {
            f >> elem[j];
          }
          mesh.getElements<2>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }
        case 5: // triangle for 2D
        {
          hrleVectorType<unsigned, 3> elem;
          for (unsigned j = 0; j < number_nodes; ++j) {
            f >> elem[j];
          }
          mesh.getElements<3>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }

        case 10: // tetra for 3D
        {
          hrleVectorType<unsigned, 4> elem;
          for (unsigned j = 0; j < number_nodes; ++j) {
            f >> elem[j];
          }
          mesh.getElements<4>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }

        case 9: // this is a quad, so just plit it into two triangles
        {
          hrleVectorType<unsigned, 3> elem;
          for (unsigned j = 0; j < 3; ++j) {
            f >> elem[j];
          }
          mesh.getElements<3>().push_back(
              elem); // push the first three nodes as a triangle
          materials.push_back(cell_material);

          f >> elem[1]; // replace middle element to create other triangle
          mesh.getElements<3>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }

        default:
          std::ostringstream oss;
          oss << "VTK Cell type " << cell_type
              << " is not supported. Cell ignored..." << std::endl;
          lsMessage::getInstance().addWarning(oss.str()).print();
        }
      } else {
        std::ostringstream oss;
        oss << "INVALID CELL TYPE! Expected number of nodes: " << number_nodes
            << ", Found number of nodes: " << elems_fake
            << "; Ignoring element...";
        lsMessage::getInstance().addError(oss.str());
        // ignore rest of lines
        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
    }

    f_ct.close();
    f_m.close();
    f.close();
  }
};

#endif // LS_VTK_READER_HPP
