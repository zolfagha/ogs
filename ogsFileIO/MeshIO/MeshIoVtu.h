
#pragma once

#include <string>
#include <vector>

#include "CommonIO/VtuWriter.h"

class MeshIoVtu : public VtuWriter
{
public:
    explicit MeshIoVtu(bool binary_mode);
    virtual ~MeshIoVtu(void){}

public:
    bool WriteXMLUnstructuredGrid(const std::string &vtkfile,
                                  MeshLib::IMesh &msh);

protected:
    bool WriteElementGroupID(std::fstream &fin,
                                  bool output_data,
                                  MeshLib::IMesh& m_msh,
                                  long &offset);

};
