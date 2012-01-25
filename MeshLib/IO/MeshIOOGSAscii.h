
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

    TemplateUnstructuredMesh *msh = NULL;

    while (!if_file.eof())
    {
      getline(if_file, line_string);

      if (line_string.find("#STOP") != std::string::npos)
        break;

      if (line_string.find("#FEM_MSH")!=std::string::npos) {
        msh = new TemplateUnstructuredMesh();
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
          newElem->setID(idx);
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
} // end namespace

