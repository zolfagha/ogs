
#pragma once

#include <string>
#include <vector>

#include "MeshIO/MeshIoVtu.h"


class OutputIoVtu : public MeshIoVtu
{
public:
    typedef std::pair<std::string, double*> NodalValue;
    typedef std::pair<std::string, double*> ElementalValue;

    explicit OutputIoVtu(bool binary_mode);
    virtual ~OutputIoVtu(void){}

public:
    bool Write(const std::string &vtkfile,
                                  MeshLib::IMesh &msh, std::vector<NodalValue> &nodal, std::vector<ElementalValue> &elemental);

protected:
    inline bool WriteNodalValue(std::fstream &fin,
                                bool output_data,
                                std::vector<NodalValue> &nodal,
                                MeshLib::IMesh &m_msh,
                                long &offset);
    inline bool WriteElementValue(std::fstream &fin,
                                  bool output_data,
                                  std::vector<ElementalValue> &elemental,
                                    MeshLib::IMesh &m_msh,
                                  long &offset);

};
