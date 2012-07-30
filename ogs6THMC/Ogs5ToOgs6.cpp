
#include "Ogs5ToOgs6.h"

#include "logog/include/logog.hpp"

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"
#include "MaterialLib/Fluid.h"
#include "SolutionLib/FemProblem/FemDirichletBC.h"
#include "SolutionLib/FemProblem/FemNeumannBC.h"
#include "SolutionLib/FemProblem/AbstractFemIVBVProblem.h"
#include "GeoProcessBuilder.h"

using namespace ogs5;

namespace ogs6
{

namespace Ogs5ToOgs6
{

void convertTimeStepping(const std::vector<CTimeDiscretization*> &td_vector, std::vector<NumLib::ITimeStepFunction*> &tf_vector)
{
    bool done_shared = false;
    for (size_t i=0; i<td_vector.size(); i++)
    {
        const CTimeDiscretization* td = td_vector[i];
        if (td->time_independence || !done_shared) {
            NumLib::ITimeStepFunction* tf = new NumLib::TimeStepFunctionVector(td->time_start, td->time_end, td->time_step_vector);
            tf_vector.push_back(tf);
            if (!td->time_independence)
                done_shared = true;
        }
    }
}

void convertFluidProperty(const CFluidProperties &mfp, MaterialLib::Fluid &fluid)
{
    if (mfp.density_model==1) {
        fluid.density = new NumLib::TXFunctionConstant(mfp.rho_0);
    }

    if (mfp.viscosity_model==1) {
        fluid.dynamic_viscosity = new NumLib::TXFunctionConstant(mfp.my_0);
    }

}

void convertSolidProperty(const CSolidProperties &msp, MaterialLib::Solid &solid)
{
    if (msp.Density_mode==1) {
        solid.density = new NumLib::TXFunctionConstant((*msp.data_Density)(0));
    }

    solid.poisson_ratio = new NumLib::TXFunctionConstant(msp.PoissonRatio);

    if (msp.Youngs_mode==1) {
        solid.Youngs_modulus = new NumLib::TXFunctionConstant((*msp.data_Youngs)(0));
    }

}

void convertPorousMediumProperty(const CMediumProperties &mmp, MaterialLib::PorousMedia &pm)
{
    if (mmp.permeability_model==1) {
        pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(mmp.permeability_tensor[0]);
    }

    if (mmp.porosity_model==1) {
        pm.porosity = new NumLib::TXFunctionConstant(mmp.porosity);
    }

    if (mmp.storage_model==1) {
        pm.storage = new NumLib::TXFunctionConstant(mmp.storage_model_values[0]);
    }
}

std::string convertLinearSolverType(int ls_method)
{
    std::string ls_sol_type = "";
    switch (ls_method) {
    case 1:
        ls_sol_type = "GAUSS";
        break;
    case 2:
        ls_sol_type = "BICGSTAB";
        break;
    case 3:
        ls_sol_type = "BICG";
        break;
    case 5:
        ls_sol_type = "CG";
        break;
    case 805:
        ls_sol_type = "PARDISO";
        break;
    default:
        break;
    } 
    return ls_sol_type;
}

std::string convertLinearSolverPreconType(int ls_precon)
{
    std::string str = "";
    switch (ls_precon) {
    case 0:
        str = "NONE";
        break;
    case 1:
        str = "JACOBI";
        break;
    case 100:
        str = "ILU";
        break;
    default:
        break;
    } 
    return str;
}


void convert(const Ogs5FemData &ogs5fem, Ogs6FemData &ogs6fem, BaseLib::Options &option)
{

    // -------------------------------------------------------------------------
    // Materials
    // -------------------------------------------------------------------------
    // MFP
    const size_t n_fluid = ogs5fem.mfp_vector.size();
    for (size_t i=0; i<n_fluid; i++)
    {
        CFluidProperties* mfp = ogs5fem.mfp_vector[i];

        MaterialLib::Fluid* fluid = new MaterialLib::Fluid();
        ogs6fem.list_fluid.push_back(fluid);

        convertFluidProperty(*mfp, *fluid);
    }

    // MSP
    for (size_t i=0; i<ogs5fem.msp_vector.size(); i++)
    {
        CSolidProperties* msp = ogs5fem.msp_vector[i];

        MaterialLib::Solid* solid = new MaterialLib::Solid();
        ogs6fem.list_solid.push_back(solid);

        convertSolidProperty(*msp, *solid);
    }

    // MMP
    for (size_t i=0; i<ogs5fem.mmp_vector.size(); i++)
    {
        CMediumProperties* mmp = ogs5fem.mmp_vector[i];
        MaterialLib::PorousMedia* pm = new MaterialLib::PorousMedia();
        ogs6fem.list_pm.push_back(pm);
        convertPorousMediumProperty(*mmp, *pm);
    }

    // -------------------------------------------------------------------------
    // Geometry
    // -------------------------------------------------------------------------
    ogs6fem.geo_unique_name = ogs5fem.geo_unique_name;
    ogs6fem.geo = ogs5fem.geo_obj;

    // -------------------------------------------------------------------------
    // Mesh
    // -------------------------------------------------------------------------
    for (size_t i=0; i<ogs5fem.list_mesh.size(); i++) {
        ogs6fem.list_dis_sys.push_back(new DiscreteLib::DiscreteSystem(*ogs5fem.list_mesh[i]));
    }

    // -------------------------------------------------------------------------
    // Time group
    // -------------------------------------------------------------------------
    // TIM
    convertTimeStepping(ogs5fem.time_vector, ogs6fem.list_tim);

    // -------------------------------------------------------------------------
    // User-defined curves
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Individual process and IVBV
    // -------------------------------------------------------------------------
    // PCS
    BaseLib::Options* optPcsData = option.addSubGroup("ProcessData");
    for (size_t i=0; i<ogs5fem.pcs_vector.size(); i++)
    {
        CRFProcess* rfpcs = ogs5fem.pcs_vector[i];
        std::string pcs_name = FiniteElement::convertProcessTypeToString(rfpcs->getProcessType());
//        ProcessLib::Process* pcs6 = GeoProcessBuilder::getInstance()->create(pcs_name);
//        if (pcs6==0) {
//            LOGOG_CERR << " Error: Process not found - " << pcs_name << std::endl;
//            continue;
//        }
//        ogs6fem.list_pcs.insert(pcs_name, pcs6);

        BaseLib::Options* optPcs = optPcsData->addSubGroup(pcs_name);

        //Mesh
        optPcs->addOptionAsNum("MeshID", rfpcs->mesh_id);

        //Time
        optPcs->addOptionAsNum("TimeGroupID", rfpcs->timegroup_id);

        // IC
        BaseLib::Options* optIcList = optPcs->addSubGroup("ICList");
        for (size_t i=0; i<ogs5fem.ic_vector.size(); i++)
        {
            CInitialCondition* rfic = ogs5fem.ic_vector[i];
            std::string bc_pcs_name = FiniteElement::convertProcessTypeToString(rfic->getProcessType());
            if (bc_pcs_name.compare(pcs_name)==0) {
                BaseLib::Options* optIc = optIcList->addSubGroup("IC");
                optIc->addOption("Variable", rfic->primaryvariable_name);
                optIc->addOption("GeometryType", rfic->geo_type_name);
                optIc->addOption("GeometryName", rfic->geo_name);
                optIc->addOption("DistributionType", FiniteElement::convertDisTypeToString(rfic->getProcessDistributionType()));
                optIc->addOptionAsNum("DistributionValue", rfic->geo_node_value);
            }
        }

        // BC
        BaseLib::Options* optBcList = optPcs->addSubGroup("BCList");
        for (size_t i=0; i<ogs5fem.bc_vector.size(); i++)
        {
            CBoundaryCondition* rfbc = ogs5fem.bc_vector[i];
            std::string bc_pcs_name = FiniteElement::convertProcessTypeToString(rfbc->getProcessType());
            if (bc_pcs_name.compare(pcs_name)==0) {
                BaseLib::Options* optBc = optBcList->addSubGroup("BC");
                optBc->addOption("Variable", rfbc->primaryvariable_name);
                optBc->addOption("GeometryType", rfbc->geo_type_name);
                optBc->addOption("GeometryName", rfbc->geo_name);
                optBc->addOption("DistributionType", FiniteElement::convertDisTypeToString(rfbc->getProcessDistributionType()));
                optBc->addOptionAsNum("DistributionValue", rfbc->geo_node_value);
            }
        }

        //ST
        BaseLib::Options* optStList = optPcs->addSubGroup("STList");
        for (size_t i=0; i<ogs5fem.st_vector.size(); i++)
        {
            CSourceTerm* rfst = ogs5fem.st_vector[i];
            std::string st_pcs_name = FiniteElement::convertProcessTypeToString(rfst->getProcessType());
            if (st_pcs_name.compare(pcs_name)==0) {
                BaseLib::Options* optSt = optStList->addSubGroup("ST");
                optSt->addOption("Variable", rfst->primaryvariable_name);
                optSt->addOption("GeometryType", rfst->geo_type_name);
                optSt->addOption("GeometryName", rfst->geo_name);
                std::string ogs5distype = FiniteElement::convertDisTypeToString(rfst->getProcessDistributionType());
                std::string ogs6distype = ogs5distype;
                if (ogs5distype.find("_NEUMANN")!=std::string::npos) {
                    optSt->addOption("STType", "NEUMANN");
                    ogs6distype = BaseLib::replaceString("_NEUMANN", "", ogs6distype);
                } else {
                    optSt->addOption("STType", "SOURCESINK");
                }
                optSt->addOption("DistributionType", ogs6distype);
                optSt->addOptionAsNum("DistributionValue", rfst->geo_node_value);
            }
        }

        // NUM
        BaseLib::Options* optNum = optPcs->addSubGroup("Numerics");
        for (size_t i=0; i<ogs5fem.num_vector.size(); i++)
        {
            CNumerics* rfnum = ogs5fem.num_vector[i];
            if (rfnum->pcs_type_name.compare(pcs_name)==0) {
                // linear solver
                BaseLib::Options* optLS = optNum->addSubGroup("LinearSolver");
                optLS->addOption("solver_type", convertLinearSolverType(rfnum->ls_method));
                optLS->addOption("precon_type", convertLinearSolverPreconType(rfnum->ls_precond));
                //optLS->addOption("error_type", "NONE");
                optLS->addOptionAsNum("error_tolerance", rfnum->ls_error_tolerance);
                optLS->addOptionAsNum("max_iteration_step", rfnum->ls_max_iterations);

                // nonlinear solver
                BaseLib::Options* optNS = optNum->addSubGroup("NonlinearSolver");
                optNS->addOption("solver_type", rfnum->nls_method_name);
                optNS->addOptionAsNum("error_type", rfnum->nls_error_method);
                optNS->addOptionAsNum("error_tolerance", rfnum->nls_error_tolerance);
                optNS->addOptionAsNum("max_iteration_step", rfnum->nls_max_iterations);

                // other stuff
                optNum->addOptionAsNum("TimeTheta", rfnum->ls_theta);
                optNum->addOptionAsNum("GaussPoint", rfnum->ele_gauss_points);
                optNum->addOptionAsNum("FEM_FCT", rfnum->fct_method);

            }
        }

        // OUT
        BaseLib::Options* optOut = optPcs->addSubGroup("OutputList");
        for (size_t i=0; i<ogs5fem.out_vector.size(); i++)
        {
            COutput* rfout = ogs5fem.out_vector[i];
            std::string out_pcs_name = FiniteElement::convertProcessTypeToString(rfout->getProcessType());
            if (out_pcs_name.compare(pcs_name)==0) {
                BaseLib::Options* opt = optOut->addSubGroup("Output");
                opt->addOption("DataType", rfout->dat_type_name);
                opt->addOption("GeometryType", rfout->geo_type);
                opt->addOption("GeometryName", rfout->geo_name);
                opt->addOption("TimeType", rfout->tim_type_name);
                opt->addOptionAsNum("TimeSteps", rfout->nSteps);
                std::vector<std::string> values(rfout->_nod_value_vector);
                values.insert(values.end(), rfout->_ele_value_vector.begin(), rfout->_ele_value_vector.end());
                opt->addOptionAsArray("Variables", values);
                opt->addOptionAsArray("MMPValues", rfout->mmp_value_vector);
                opt->addOptionAsArray("MFPValues", rfout->mfp_value_vector);
            }
        }

    }


    // -------------------------------------------------------------------------
    // Coupling
    // -------------------------------------------------------------------------

}

};
} //end
