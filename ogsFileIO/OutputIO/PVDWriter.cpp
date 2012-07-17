
#include "PVDWriter.h"

#include <fstream>

const std::string INDEX_STR = "  ";

bool PVDWriter::InitializePVD(const std::string &file_base_name)
{
	//PVD
	this->pvd_data.vec_dataset.clear();
	this->pvd_file_name = file_base_name;
	this->pvd_file_name += ".pvd";
//	//VTK
//	int ibs = (int)file_base_name.find_last_of("\\");
//	int is = (int)file_base_name.find_last_of("/");
//	if (ibs != (int)string::npos  || is != (int)string::npos)
//	{
//		int ibegin = ibs;
//		if (is > ibs)
//			ibegin = is;
//		ibegin += 1;
//		this->pvd_vtk_file_name_base = file_base_name.substr(ibegin);
//        this->pvd_vtk_file_path_base = file_base_name.substr(0, ibegin);
//	}
//	else
//    {
//		this->pvd_vtk_file_name_base = file_base_name;
//        this->pvd_vtk_file_path_base = "";
//    }
//	if (pcs_type_name.size() > 0)        // PCS
//		this->pvd_vtk_file_name_base += "_" + pcs_type_name;

	return true;
}

bool PVDWriter::WriteHeaderOfPVD(std::fstream &fin)
{
	fin << "<?xml version=\"1.0\"?>" << std::endl;
	fin <<
	"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"
	    << std::endl;
	fin << INDEX_STR << "<Collection>" << std::endl;
	return true;
}

bool PVDWriter::WriteEndOfPVD(std::fstream &fin)
{
	fin << INDEX_STR << "</Collection>" << std::endl;
	fin << "</VTKFile>" << std::endl;
	return true;
}

bool PVDWriter::WriteDatasetOfPVD(std::fstream &fin, double timestep, const std::string &vtkfile)
{
	fin.setf(std::ios::scientific,std::ios::floatfield);
	fin.precision(12);
	fin << INDEX_STR << INDEX_STR << "<DataSet timestep=\"" << timestep <<
	"\" group=\"\" part=\"0\" file=\"" << vtkfile << "\"/>" << std::endl;
	return true;
}

bool PVDWriter::UpdatePVD(const std::string &pvdfile, const PVDData &pvd_data)
{
	std::fstream fin(pvdfile.data(), std::ios::out);
	if (!fin.good())
		return false;

	PVDWriter::WriteHeaderOfPVD(fin);
	for (int i = 0; i < (int)pvd_data.vec_dataset.size(); i++)
		PVDWriter::WriteDatasetOfPVD(fin, pvd_data.vec_dataset[i].timestep, pvd_data.vec_dataset[i].vtk_file);
	PVDWriter::WriteEndOfPVD(fin);

	fin.close();

	return true;
}


