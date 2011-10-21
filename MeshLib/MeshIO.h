
#pragma once
#include <iostream>
#include <fstream>
#include <string>

#include "Mesh.h"

namespace MeshLib
{

class MeshIOOGS
{
protected:
  static const std::string convertElementType2String(const ElementType::type t)
  {
    if (t == ElementType::LINE)			return "line";
    if (t == ElementType::QUAD)			return "quad";
    if (t == ElementType::HEXAHEDRON)	return "hex";
    if (t == ElementType::TRIANGLE)		return "tri";
    if (t == ElementType::TETRAHEDRON)	return "tet";
    if (t == ElementType::PRISM)		return "pris";
    if (t == ElementType::PYRAMID)		return "pyra";
    return "none";
  };

  static ElementType::type convertString2ElementType(const std::string &s)
  {
    if (s.compare("line") == 0) return ElementType::LINE;
    if (s.compare("quad") == 0) return ElementType::QUAD;
    if (s.compare("hex")  == 0) return ElementType::HEXAHEDRON;
    if (s.compare("tri")  == 0) return ElementType::TRIANGLE;
    if (s.compare("tet")  == 0) return ElementType::TETRAHEDRON;
    if (s.compare("pris") == 0) return ElementType::PRISM;
    if (s.compare("pyra") == 0) return ElementType::PYRAMID;
    return ElementType::INVALID;
  };

public:
  static void readMesh(std::string const& fileName, std::vector<IMesh*> &vec_mesh) {

    std::ifstream if_file (fileName.data(),std::ios::in);
    if (!if_file.is_open())
    {
      std::cout << "CFEMesh::FEMRead() - Could not open file...\n";
      return;
    }

    std::cout << "MSHRead:  ASCII file" << std::endl;
    std::string line_string ("");

    UnstructuredMesh *msh = NULL;

    while (!if_file.eof())
    {
      getline(if_file, line_string);

      if (line_string.find("#STOP") != std::string::npos)
        break;

      if (line_string.find("#FEM_MSH")!=std::string::npos) {
        msh = new UnstructuredMesh();
        vec_mesh.push_back(msh);
        continue;
      } 
      
      if (!msh) continue;

      if (line_string.find("$NODES") != std::string::npos)
      {
        double x, y, z;
        size_t no_nodes, idx;
        if_file >> no_nodes >> std::ws;
        std::string s;
        std::ios::pos_type position = if_file.tellg();
        for (size_t i = 0; i < no_nodes; i++)
        {
          if_file >> idx >> x >> y >> z;
          msh->setNode(idx, x, y, z);
          position = if_file.tellg();
          if_file >> s;
          if_file.seekg(position, std::ios::beg);
          if_file >> std::ws;
        }
      }
      else  if (line_string.find("$ELEMENTS") != std::string::npos)
      {
        size_t no_elements, idx, group_id, idummy;
        int grid_adaptation;
        if_file >> no_elements >> std::ws;
        for (size_t i = 0; i < no_elements; i++)
        {
          std::string buffer("");
          std::string name("");
          if_file >> idx >> group_id;
          if_file >> buffer;

          if (buffer.find("-1") != std::string::npos) {
            grid_adaptation = strtol(buffer.data(), NULL, 0);
            if_file >> name;
          } else {
            name = buffer;
          }

          ElementType::type ele_type = convertString2ElementType(name);
          IElement *newElem = ElemenetFactory::createNewElement(ele_type);
          newElem->setElementID(idx);
          newElem->setGroupID(group_id);
          for (size_t j=0; j < newElem->getNumberOfNodes(); j++) {
            if_file >> idummy;
            newElem->setNodeID(j, idummy);
          }
          msh->addElement(newElem);
        }
      }
    }

    if_file.close();
  };

  static void write(IMesh const * mesh, std::string const& fileName) 
  {
    std::cout << "Not implemented yet." << std::endl;
  };
};

class MeshIOLegacyVtk
{
public:
    static void WriteAsciiFile(std::string str_file_name, const IMesh &mesh) {
        std::fstream vtk_file (str_file_name.c_str(), std::ios::out);
        vtk_file.setf(std::ios::scientific,std::ios::floatfield);
        vtk_file.precision(12);
        if (!vtk_file.good())
        {
            std::cout << "Could not open file for writing: " << str_file_name << std::endl;
            return;
        }
        vtk_file.seekg(0L,std::ios::beg);
        WriteVTKHeader(vtk_file);
        WriteVTKPoints(vtk_file, mesh);
        WriteVTKCells(vtk_file, mesh);
        WriteVTKCellDataHeader(vtk_file, mesh);
        WriteVTKMaterialID(vtk_file, mesh);

        vtk_file.close();
    };

protected:
    static void WriteVTKHeader(std::fstream &vtk_file)
    {
        vtk_file << "# vtk DataFile Version 3.0" << std::endl;
        vtk_file << "Unstructured Grid from OpenGeoSys" << std::endl;
        vtk_file << "ASCII"  << std::endl;
        vtk_file << "DATASET UNSTRUCTURED_GRID"  << std::endl;
    }


    static void WriteVTKPoints(std::fstream &vtk_file, const IMesh &mesh)
    {
        const size_t n_nodes(mesh.getNumberOfNodes());
        vtk_file << "POINTS "<< n_nodes << " double" << std::endl;

        double pt[3]={};
        for(size_t i = 0; i < n_nodes ; i++)
        {
            mesh.getNodeCoordinates(i, pt);
            vtk_file << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
        }
    }

    static void WriteVTKCells(std::fstream &vtk_file, const IMesh &mesh)
    {
        size_t numCells = mesh.getNumberOfElements();

        // count overall length of element vector
        long numAllPoints =0;
        for(size_t i=0; i < numCells; i++)
        {
            IElement* ele = mesh.getElemenet(i);
            numAllPoints = numAllPoints + (ele->getNumberOfNodes()) + 1;
        }

        // write elements
        vtk_file << "CELLS " << numCells << " " << numAllPoints << std::endl;
        for(size_t i=0; i < numCells; i++)
        {
            IElement* ele = mesh.getElemenet(i);

            // Write number of points per cell
            switch(ele->getElementType())
            {
            case ElementType::LINE:
                vtk_file << "2"; break;
            case ElementType::QUAD:
                vtk_file << "4"; break;
            case ElementType::HEXAHEDRON:
                vtk_file << "8"; break;
            case ElementType::TRIANGLE:
                vtk_file << "3"; break;
            case ElementType::TETRAHEDRON:
                vtk_file << "4"; break;
            case ElementType::PRISM:
                vtk_file << "6"; break;
            case ElementType::PYRAMID:
                vtk_file << "5"; break;
            default:
                std::cerr << "COutput::WriteVTKElementData MshElemType not handled" << std::endl;
                break;
            }

            for(size_t j = 0; j < ele->getNumberOfNodes(); j++)
                vtk_file << " " << ele->getNodeID(j);

            vtk_file << std::endl;
        }
        vtk_file << std::endl;

        // write cell types
        vtk_file << "CELL_TYPES " << numCells << std::endl;

        for(size_t i=0; i < numCells; i++)
        {
            IElement* ele = mesh.getElemenet(i);

            // Write vtk cell type number (see vtkCellType.h)
            switch(ele->getElementType())
            {
            case ElementType::LINE:
                vtk_file << "3" << std::endl; break;
            case ElementType::QUAD:
                vtk_file << "9" << std::endl; break;
            case ElementType::HEXAHEDRON:
                vtk_file << "12" << std::endl; break;
            case ElementType::TRIANGLE:
                vtk_file << "5" << std::endl; break;
            case ElementType::TETRAHEDRON:
                vtk_file << "10" << std::endl; break;
            case ElementType::PRISM: // VTK_WEDGE
                vtk_file << "13" << std::endl; break;
            case ElementType::PYRAMID:
                vtk_file << "14" << std::endl; break;
            default:
                std::cerr << "COutput::WriteVTKElementData MshElemType not handled" << std::endl;
                break;
            }
        }
        vtk_file << std::endl;
    }

    static void WriteVTKCellDataHeader(std::fstream &vtk_file, const IMesh &mesh)
    {
        vtk_file << "CELL_DATA " << mesh.getNumberOfElements() << std::endl;
    }

    static void WriteVTKMaterialID(std::fstream &vtk_file, const IMesh &mesh)
    {
        vtk_file << "SCALARS " << "GroupID" << " int 1" << std::endl;
        vtk_file << "LOOKUP_TABLE default" << std::endl;
        for (size_t i = 0; i < mesh.getNumberOfElements(); i++) {
            IElement* e = mesh.getElemenet(i);
            vtk_file << e->getGroupID() << std::endl;
        }
    }
};

typedef std::pair<std::string, double*> NodalScalarValue;

class MeshIOLegacyVtk4Simulation : public MeshIOLegacyVtk 
{
public:
    static void WriteAsciiFile(std::string str_file_name, const IMesh &mesh, int time_step_number, double simulation_time, std::vector<NodalScalarValue> &nodalScalarValues) {
        std::fstream vtk_file (str_file_name.c_str(), std::ios::out);
        vtk_file.setf(std::ios::scientific,std::ios::floatfield);
        vtk_file.precision(12);
        if (!vtk_file.good())
        {
            std::cout << "Could not open file for writing: " << str_file_name << std::endl;
            return;
        }
        vtk_file.seekg(0L,std::ios::beg);
        WriteVTKHeader(vtk_file, time_step_number, simulation_time, nodalScalarValues);
        WriteVTKPoints(vtk_file, mesh);
        WriteVTKCells(vtk_file, mesh);
        WriteVTKDataArrays(vtk_file, mesh, nodalScalarValues);

        vtk_file.close();
    };

protected:
    static void WriteVTKHeader(std::fstream &vtk_file, int time_step_number, double simulation_time, std::vector<NodalScalarValue> &nodalScalarValues)
    {
        MeshIOLegacyVtk::WriteVTKHeader(vtk_file);
        // time information
        // see http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
        vtk_file << "FIELD TimesAndCycles 2" << std::endl;
        vtk_file << "TIME 1 1 double" << std::endl;
        vtk_file << simulation_time << std::endl;
        vtk_file << "CYLCE 1 1 long" << std::endl;
        vtk_file << time_step_number << std::endl;
    }

    static void WriteVTKDataArrays(std::fstream &vtk_file, const IMesh &mesh, std::vector<NodalScalarValue> &nodalScalarValues)
    {
        const size_t n_nodes = mesh.getNumberOfNodes();
        vtk_file << "POINT_DATA " << n_nodes << std::endl;
        for (size_t i=0; i<nodalScalarValues.size(); i++) {
            vtk_file << "SCALARS " << nodalScalarValues.at(i).first << " double 1" << std::endl;
            vtk_file << "LOOKUP_TABLE default" << std::endl;
            for (size_t j = 0; j<n_nodes; j++)
            {
                vtk_file << nodalScalarValues.at(i).second[j] << std::endl;
            }
        }
    }

};

} // end namespace

