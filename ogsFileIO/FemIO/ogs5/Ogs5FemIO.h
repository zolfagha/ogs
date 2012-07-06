
#pragma once

#include <string>
#include <vector>

#include "rf_pcs.h"
#include "rf_mfp_new.h"
#include "rfmat_cp.h"
#include "rf_bc_new.h"
#include "rf_st_new.h"
#include "rf_ic_new.h"
#include "rf_out_new.h"
#include "rf_tim_new.h"
#include "rf_msp_new.h"
#include "rf_mmp_new.h"
#include "rf_num_new.h"

namespace ogs5
{

class Ogs5FemData
{

public:
	Ogs5FemData() {};
	~Ogs5FemData();
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
