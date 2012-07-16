
#include "Ogs5FemIO.h"

#include "logog/include/logog.hpp"

#include "GeoIO/Ogs4GeoIO.h"
#include "MeshIO/MeshIoOgs5.h"


#include "BaseLib/CodingTools.h"

namespace ogs5
{

Ogs5FemData::~Ogs5FemData()
{
	BaseLib::releaseObjectsInStdVector(pcs_vector);
	BaseLib::releaseObjectsInStdVector(mfp_vector);
	BaseLib::releaseObjectsInStdVector(msp_vector);
	BaseLib::releaseObjectsInStdVector(mmp_vector);
	BaseLib::releaseObjectsInStdVector(cp_vector);
	BaseLib::releaseObjectsInStdVector(bc_vector);
	BaseLib::releaseObjectsInStdVector(st_vector);
	BaseLib::releaseObjectsInStdVector(ic_vector);
	BaseLib::releaseObjectsInStdVector(out_vector);
	BaseLib::releaseObjectsInStdVector(time_vector);
	BaseLib::releaseObjectsInStdVector(num_vector);
	//geo, msh objects are passed
}

bool Ogs5FemData::read(const std::string &proj_path)
{
	PCSRead(proj_path, pcs_vector);
	MFPRead(proj_path, mfp_vector);
	MSPRead(proj_path, msp_vector);
	MMPRead(proj_path, mmp_vector);
	CPRead(proj_path, cp_vector);
	BCRead(proj_path, bc_vector);
	STRead(proj_path, st_vector);
	ICRead(proj_path, ic_vector);
	OUTRead(proj_path, out_vector);
	TIMRead(proj_path, time_vector);
	NUMRead(proj_path, num_vector);

	// gli
	std::vector<std::string> geo_errors;
	this->geo_obj = new GeoLib::GEOObjects();
	Ogs4GeoIO::readGLIFileV4(proj_path+".gli", this->geo_obj, geo_unique_name, geo_errors);
	if (geo_errors.size()>0) {
		LOGOG_CERR << " Error when reading GLI" << std::endl;
		for (size_t i=0; i<geo_errors.size(); i++)
			LOGOG_CERR << " -> " << geo_errors[i] << std::endl;
		return false;
	}

	// mesh
	MeshIoOgs5::readMesh(proj_path+".msh", list_mesh);
	if (list_mesh.size()==0) {
		LOGOG_CERR << " Error: Cannot find mesh - " << proj_path+".msh" << std::endl;
		return false;
	}

	return true;


//	RCRead(proj_path);
//	KRRead(proj_path, geo_obj, unique_name);
//	KRWrite(proj_path);
//
//
//	PCTRead(proj_path);                   // PCH
//	FMRead(proj_path);                    // PCH
//	FCTRead(proj_path);                   //OK
//	CURRead(proj_path);                   //OK
}

}
