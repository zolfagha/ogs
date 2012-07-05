
#include "OGS5Data.h"

#include "ogs5/rf_pcs.h"
#include "ogs5/rf_mfp_new.h"
#include "ogs5/rfmat_cp.h"
#include "ogs5/rf_bc_new.h"
#include "ogs5/rf_st_new.h"
#include "ogs5/rf_ic_new.h"
#include "ogs5/rf_out_new.h"
#include "ogs5/rf_tim_new.h"
#include "ogs5/rf_msp_new.h"
#include "ogs5/rf_mmp_new.h"
#include "ogs5/rf_num_new.h"

#include "BaseLib/CodingTools.h"

namespace ogs6
{
OGS5Data::~OGS5Data()
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
}

void OGS5Data::read(const std::string &proj_path)
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

	//	std::vector<CFEMesh*> mesh_vec;
	//	FEMRead(proj_path, mesh_vec, &geo_obj, &unique_name);

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
