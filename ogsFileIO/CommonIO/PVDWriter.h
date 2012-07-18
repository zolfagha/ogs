
#pragma once

#include <string>
#include <fstream>

#include "PVDData.h"

class PVDWriter
{
public:
    PVDWriter(void){}
    virtual ~PVDWriter(void){}

    bool initialize(const std::string &file_base_namepcs_type_name);
    bool update(double time, const std::string &vtk_file);

private:
    bool writeHeader(std::fstream &fin);
    bool writeEnd(std::fstream &fin);
    bool writeDataset(std::fstream &fin, double timestep, const std::string &vtkfile);

private:
    PVDData _pvd_data;
    std::string _pvd_file_name;
};
