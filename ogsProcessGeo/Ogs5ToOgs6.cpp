
#include "Ogs5ToOgs6.h"

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
		pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(mmp.permeability);
	}

	if (mmp.porosity_model==1) {
		pm.porosity = new NumLib::TXFunctionConstant(mmp.porosity);
	}

	if (mmp.storage_model==1) {
		pm.storage = new NumLib::TXFunctionConstant(mmp.storage_model_values[0]);
	}
}

const GeoLib::GeoObject* searchGeoByName(const GeoLib::GEOObjects &geo_obj, const std::string &unique_geo_name, const std::string &geo_type_name, const std::string &geo_name)
{
	if (geo_type_name.find("POINT") != std::string::npos)
	{
		const GeoLib::PointVec* pnt_vec(geo_obj.getPointVecObj(unique_geo_name));
		if (pnt_vec)
		{
			const GeoLib::Point* pnt(pnt_vec->getElementByName(geo_name));
			if (pnt == NULL)
			{
				std::cerr << "ERROR in GeoIO::readGeoInfo: point name \""
						  << geo_name << "\" not found!" << std::endl;
				exit(1);
			}
			return pnt;
		}

		std::cerr << "Error in GeoIO::readGeoInfo: point vector not found!" <<std::endl;
		exit(1);
	}

	else if (geo_type_name.find("POLYLINE") != std::string::npos)
	{
		const GeoLib::PolylineVec* ply_vec = geo_obj.getPolylineVecObj(unique_geo_name);
		if (ply_vec)
		{
			const GeoLib::Polyline* ply(ply_vec->getElementByName(geo_name));
			if (ply == NULL)
			{
				std::cerr << "error in GeoIO::readGeoInfo: polyline name \""
						  << geo_name << "\" not found!" << std::endl;
				exit(1);
			}
			return ply;
		}

		std::cerr << "Error in GeoIO::readGeoInfo: polyline vector not found!" <<std::endl;
		exit(1);
	}

	else if (geo_type_name.find("SURFACE") != std::string::npos)
	{
	GeoLib::SurfaceVec const* sfc_vec (geo_obj.getSurfaceVecObj(unique_geo_name));
		if (sfc_vec)
		{
			const GeoLib::Surface* sfc(sfc_vec->getElementByName(geo_name));
			if (sfc == NULL)
			{
				std::cerr << "Error in GeoIO::readGeoInfo: surface name \""
						  << geo_name << "\" not found!" << std::endl;
				exit(1);
			}
			return sfc;
		}

		std::cerr << "Error in GeoIO::readGeoInfo: surface vector not found!" <<std::endl;
		exit(1);
	}

	else if (geo_type_name.find("VOLUME") != std::string::npos)
	{
		return NULL;
	}

	else if (geo_type_name.find("DOMAIN") != std::string::npos)
	{
		return NULL;
	}

	return NULL;
}

void convert(const Ogs5FemData &ogs5fem, const GeoLib::GEOObjects &geo, const std::string &unique_geo_name, const DiscreteLib::DiscreteSystem &dis)
{
	std::vector<ProcessLib::Process*> list_pcs;
	std::vector<SolutionLib::AbstractFemIVBVProblem*> list_problem;
	std::vector<MaterialLib::PorousMedia*> list_pm;
	std::vector<MaterialLib::Solid*> list_solid;
	std::vector<MaterialLib::Fluid*> list_fluid;
	std::vector<NumLib::ITimeStepFunction*> list_tim;

	// MFP
	const size_t n_fluid = ogs5fem.mfp_vector.size();
	for (size_t i=0; i<n_fluid; i++)
	{
		CFluidProperties* mfp = ogs5fem.mfp_vector[i];

		MaterialLib::Fluid* fluid = new MaterialLib::Fluid();
		list_fluid.push_back(fluid);

		convertFluidProperty(*mfp, *fluid);
	}

	// MMP, MSP
	const size_t n_mat = ogs5fem.mmp_vector.size();
	for (size_t i=0; i<n_mat; i++)
	{
		CSolidProperties* msp = ogs5fem.msp_vector[i];
		CMediumProperties* mmp = ogs5fem.mmp_vector[i];

		MaterialLib::Solid* solid = new MaterialLib::Solid();
		list_solid.push_back(solid);
		MaterialLib::PorousMedia* pm = new MaterialLib::PorousMedia();
		list_pm.push_back(pm);

		convertSolidProperty(*msp, *solid);
		convertPorousMediumProperty(*mmp, *pm);
	}


	// IVBV problem
    FemLib::LagrangianFeObjectContainer* _feObjects = new FemLib::LagrangianFeObjectContainer(*dis.getMesh());

    // PCS
	for (size_t i=0; i<ogs5fem.pcs_vector.size(); i++)
	{
		CRFProcess* rfpcs = ogs5fem.pcs_vector[i];
		std::string pcs_name = FiniteElement::convertProcessTypeToString(rfpcs->getProcessType());
		ProcessLib::Process* pcs6 = GeoProcessBuilder::getInstance()->create(pcs_name);
		list_pcs.push_back(pcs6);

		if (pcs_name.find("GROUNDWATER_FLOW")!=std::string::npos) {
			list_problem.push_back(0);
		}
	}

	SolutionLib::AbstractFemIVBVProblem* femProblem = 0;
	std::map<std::string, SolutionLib::FemVariable*> mapName2Var;

	// BC
	for (size_t i=0; i<ogs5fem.bc_vector.size(); i++)
	{
		CBoundaryCondition* rfbc = ogs5fem.bc_vector[i];
		std::string pcs_name = FiniteElement::convertProcessTypeToString(rfbc->getProcessType());
		SolutionLib::FemVariable* var = 0;
		if (mapName2Var.count(rfbc->primaryvariable_name)==0)
			var = femProblem->addVariable(rfbc->primaryvariable_name);
		else
			var = mapName2Var[rfbc->primaryvariable_name];

		const GeoLib::GeoObject* geo_obj = searchGeoByName(geo, unique_geo_name, rfbc->geo_type_name, rfbc->geo_name);
		if (rfbc->getProcessDistributionType()==ogs5::FiniteElement::DistributionType::CONSTANT)
		{
			var->addDirichletBC(new SolutionLib::FemDirichletBC(dis.getMesh(), geo_obj, new NumLib::TXFunctionConstant(.0)));
		}
	}

	//ST
	for (size_t i=0; i<ogs5fem.st_vector.size(); i++)
	{
		CSourceTerm* rfst = ogs5fem.st_vector[i];
		std::string pcs_name = FiniteElement::convertProcessTypeToString(rfst->getProcessType());
		SolutionLib::FemVariable* var = 0;
		if (mapName2Var.count(rfst->primaryvariable_name)==0)
			var = femProblem->addVariable(rfst->primaryvariable_name);
		else
			var = mapName2Var[rfst->primaryvariable_name];

		const GeoLib::GeoObject* geo_obj = searchGeoByName(geo, unique_geo_name, rfst->geo_type_name, rfst->geo_name);
		if (rfst->getProcessDistributionType()==ogs5::FiniteElement::DistributionType::CONSTANT)
		{
			var->addNeumannBC(new SolutionLib::FemNeumannBC(dis.getMesh(), _feObjects, geo_obj, new NumLib::TXFunctionConstant(.0)));
		}
	}

	// TIM
	convertTimeStepping(ogs5fem.time_vector, list_tim);

	// NUM

	// OUT


}

};
} //end
