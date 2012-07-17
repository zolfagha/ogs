
#pragma once

#include <string>
#include <vector>


#include "PVDData.h"

class PVDWriter
{
public:
	PVDData pvd_data;
	std::string pvd_file_name;
	std::string pvd_vtk_file_name_base;
    std::string pvd_vtk_file_path_base;

public:
    PVDWriter(void){}
	virtual ~PVDWriter(void){}

protected:
	bool WriteHeaderOfPVD(std::fstream &fin);
	bool WriteEndOfPVD(std::fstream &fin);
	bool WriteDatasetOfPVD(std::fstream &fin, double timestep, const std::string &vtkfile);

public:
	bool InitializePVD(const std::string &file_base_namepcs_type_name);
	bool UpdatePVD(const std::string &pvdfile, const PVDData &pvd_data);

};
