
#pragma once

#include <string>
#include <vector>

class CRFProcess;
class CFluidProperties;
class CompProperties;
class CBoundaryCondition;
class CSourceTerm;
class CInitialCondition;
class COutput;
class CTimeDiscretization;
class CSolidProperties;
class CMediumProperties;
class CNumerics;


namespace ogs6
{

class OGS5Data
{

public:
	OGS5Data() {};
	~OGS5Data();
	void read(const std::string &proj_path);

public:
	std::vector<CRFProcess*> pcs_vector;
	std::vector<CFluidProperties*> mfp_vector;
	std::vector<CompProperties*> cp_vector;
	std::vector<CBoundaryCondition*> bc_vector;
	std::vector<CSourceTerm*> st_vector;
	std::vector<CInitialCondition*> ic_vector;
	std::vector<COutput*> out_vector;
	std::vector<CTimeDiscretization*> time_vector;
	std::vector<CSolidProperties*> msp_vector;
	std::vector<CMediumProperties*> mmp_vector;
	std::vector<CNumerics*>num_vector;
};

}
