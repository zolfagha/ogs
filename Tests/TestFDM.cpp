
#include <vector>

#include <gtest/gtest.h>

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "GeoLib/Core/Polyline.h"
#include "GeoLib/Shape/Line.h"
#include "GeoLib/Shape/Rectangle.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "MeshLib/Tools/MeshGenerator.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/Coupling/Algorithm/IConvergenceCheck.h"
#include "NumLib/Coupling/Algorithm/SerialStaggeredMethod.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientCoupling/AsyncPartitionedSystem.h"
#include "FemLib/Function/FemNorm.h"

#include "TestUtil.h"

#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Model/Head.h"
#include "Tests/Geo/Model/Velocity.h"
#include "Tests/Geo/Model/Concentration.h"
#include "Tests/ExactSolution/OgataBank.h"

#include "Tests/FDM/FdmIVBVProblem.h"
#include "Tests/FDM/SingleStepFDM.h"
#include "Tests/FDM/IStencil.h"
#include "Tests/FDM/FdmGw1D.h"
#include "Tests/FDM/FdmNorm.h"






using namespace MathLib;
using namespace GeoLib;
using namespace DiscreteLib;
using namespace NumLib;
using namespace FdmLib;
using namespace FemLib;

template <class T, class T_NORM>
class TemplateConvergenceCheck : public IConvergenceCheck
{
	T_NORM _norm;

public:
	TemplateConvergenceCheck(DiscreteLib::DiscreteSystem* dis) : _norm(dis) {};

	bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff)
	{
	    for (size_t i=0; i<vars_prev.size(); i++) {
	    	if (vars_prev.getName(i).compare("h")==0 || vars_prev.getName(i).compare("c")==0) {
		        const T* f_fem_prev = vars_prev.get<T>(i);
		        const T* f_fem_cur = vars_current.get<T>(i);
	    		//v_diff = f_fem_cur->norm_diff(*f_fem_prev);
	    		v_diff = _norm(*f_fem_prev, *f_fem_cur);
//	    	} else if (vars_prev.getName(i).compare("v")==0) {
//	    		const FEMIntegrationPointFunctionVector2d* f_fem_prev = vars_prev.get<FEMIntegrationPointFunctionVector2d>(i);
//	    		const FEMIntegrationPointFunctionVector2d* f_fem_cur = vars_current.get<FEMIntegrationPointFunctionVector2d>(i);
//	    		v_diff = f_fem_cur->norm_diff(*f_fem_prev);
	    	}
	        if (v_diff>eps) {
	            return false;
	        }
	    }
	    return true;
	}
};

typedef TemplateConvergenceCheck<FdmFunctionScalar, FdmLib::NormOfFdmNodalFunction<double> > FdmFunctionConvergenceCheck;
typedef TemplateConvergenceCheck<FemNodalFunctionScalar, FemLib::NormOfFemNodalFunction<double> > FemFunctionConvergenceCheck;


typedef FunctionHead<MathLib::CRSLisSolver> MyFunctionHead;
typedef Geo::FunctionHead<CRSLisSolver> MyFemFunctionHead;
typedef Geo::FunctionVelocity MyFemFunctionVelocity;
typedef Geo::FunctionConcentration<CRSLisSolver> MyFunctionConcentration;

GWFdmProblem* defineGWProblem4FDM(DiscreteSystem &dis, double h, GeoLib::Line &line, Geo::PorousMedia &pm)
{
    //equations
    GWFdmProblem::LinearAssemblerType* linear_assembler = new GWFdmProblem::LinearAssemblerType(h, pm);
    //IVBV problem
    GWFdmProblem* _problem = new GWFdmProblem(dis, linear_assembler);
    //BC
    size_t headId = _problem->createField();
    FdmFunctionScalar* _head = _problem->getField(headId);
    _problem->setIC(headId, *_head);
    TXFunctionConstant f1(.0);
    _problem->addDirichletBC(headId, *line.getPoint2(), false, f1);
    TXFunctionConstant f2(-1e-5);
    _problem->addNeumannBC(headId, *line.getPoint1(), false, f2);

    return _problem;
}



Geo::MassFemProblem* defineMassTransportProblem(DiscreteSystem &dis, GeoLib::Line &line, Geo::PorousMedia &pm, Geo::Compound &comp)
{
    LagrangianFeObjectContainer* _feObjects = new LagrangianFeObjectContainer(*dis.getMesh());
    //equations
    Geo::MassFemEquation::LinearAssemblerType* linear_assembler = new Geo::MassFemEquation::LinearAssemblerType(*_feObjects, pm, comp);
    Geo::MassFemEquation::ResidualAssemblerType* r_assembler = new Geo::MassFemEquation::ResidualAssemblerType(*_feObjects, pm, comp);
    Geo::MassFemEquation::JacobianAssemblerType* j_eqs = new Geo::MassFemEquation::JacobianAssemblerType(*_feObjects, pm, comp);
    Geo::MassFemEquation* eqs = new Geo::MassFemEquation(linear_assembler, r_assembler, j_eqs);
    //IVBV problem
    Geo::MassFemProblem* _problem = new Geo::MassFemProblem(&dis);
    _problem->setEquation(eqs);
    // var
    FemVariable* var = _problem->addVariable("concentration");
    //IC
    FemNodalFunctionScalar* c0 = new FemNodalFunctionScalar(dis, PolynomialOrder::Linear, 0);
    var->setIC(c0);
    //BC
    NumLib::TXFunctionConstant* f1 = new  NumLib::TXFunctionConstant(1.0);
    var->addDirichletBC(new FemDirichletBC(dis.getMesh(), line.getPoint1(), f1));

    return _problem;
}

TEST(Fdm, fdm1)
{
	try {
		const double len = 2.0;
		const size_t div = 20;
		const double h = len / div;
	    MeshLib::IMesh *msh = MeshLib::MeshGenerator::generateLineMesh(len, div, .0, .0, .0);
	    GeoLib::Line line(Point(0.0, .0, .0), Point(len, .0, .0));
	    Geo::PorousMedia pm;
	    pm.hydraulic_conductivity = new TXFunctionConstant(1.e-11);
	    DiscreteSystem dis(*msh);
	    GWFdmProblem* pGW = defineGWProblem4FDM(dis, h, line, pm);
        TimeStepFunctionConstant tim(.0, 10.0, 10.0);
        pGW->setTimeSteppingFunction(tim);
	    // options
	    BaseLib::Options options;
	    BaseLib::Options* op_lis = options.addSubGroup("Lis");
	    op_lis->addOption("solver_type", "BiCG");
	    op_lis->addOption("precon_type", "NONE");
	    op_lis->addOptionAsNum("error_tolerance", 1e-10);
	    op_lis->addOptionAsNum("max_iteration_step", 500);

		MyFunctionHead f_head;
		f_head.define(&dis, pGW, options);
        f_head.setOutputParameterName(0, "h");

	    TimeSteppingController timestepping;
	    timestepping.addTransientSystem(f_head);

	    //const double epsilon = 1.e-3;
	    timestepping.setBeginning(.0);
	    timestepping.solve(1.0);

	    const FdmFunctionScalar* r_f_head = (FdmFunctionScalar*) f_head.getOutput(0);
	    //const FEMIntegrationPointFunctionVector2d* r_f_v = apart1.getOutput<FEMIntegrationPointFunctionVector2d>(apart1.getOutputParameterID("v"));
	    const IDiscreteVector<double>* vec_h = r_f_head->getNodalValues();
	    //const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();

	    r_f_head->printout();
	    //r_f_v->printout();

	    std::vector<double> expected;
	    expected.resize(div+1);
        const double p_left = 2.e+6;
        const double p_right = .0;
        const double x_len = 2.0;
	    for (size_t i=0; i<expected.size(); i++) {
            double x = i*h;
            expected[i] = (p_right-p_left) / x_len * x + p_left;
	    }

	    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*vec_h)[0], vec_h->size());

	} catch (const char* e) {
		std::cout << "***Exception caught! " << e << std::endl;
	}

}


TEST(Fdm, fdm_fem1)
{
	try {
		const double len = 2.0;
		const size_t div = 20;
		const double h = len / div;
	    MeshLib::IMesh *msh = MeshLib::MeshGenerator::generateLineMesh(len, div, .0, .0, .0);
	    GeoLib::Line line(Point(0.0, .0, .0), Point(len, .0, .0));
	    Geo::PorousMedia pm;
	    pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(1.e-11);
	    pm.porosity = new NumLib::TXFunctionConstant(1.0);
	    Geo::Compound tracer;
	    tracer.molecular_diffusion = new NumLib::TXFunctionConstant(1.e-6);
	    DiscreteSystem dis(*msh);
	    GWFdmProblem* pGW = defineGWProblem4FDM(dis, h, line, pm);
	    Geo::MassFemProblem* pMass = defineMassTransportProblem(dis, line, pm, tracer);

        TimeStepFunctionConstant tim(.0, 1e+4, 1e+3);
        pGW->setTimeSteppingFunction(tim);
        pMass->setTimeSteppingFunction(tim);
	    // options
	    BaseLib::Options options;
	    BaseLib::Options* op_lis = options.addSubGroup("Lis");
	    op_lis->addOption("solver_type", "BiCG");
	    op_lis->addOption("precon_type", "NONE");
	    op_lis->addOptionAsNum("error_tolerance", 1e-10);
	    op_lis->addOptionAsNum("max_iteration_step", 500);
        BaseLib::Options optionsMT;
        op_lis = optionsMT.addSubGroup("Lis");
        op_lis->addOption("solver_type", "BiCG");
        op_lis->addOption("precon_type", "NONE");
        op_lis->addOptionAsNum("error_tolerance", 1e-10);
        op_lis->addOptionAsNum("max_iteration_step", 1000);

		MyFunctionHead f_head;
		f_head.define(&dis, pGW, options);
        f_head.setOutputParameterName(0, "h");
        FunctionFdmVelocity f_v;
		f_v.define(dis, pm);
        f_v.setInputParameterName(0, "h");
        f_v.setOutputParameterName(0, "v");
		MyFunctionConcentration f_c;
		f_c.define(&dis, pMass, optionsMT);
        f_c.setInputParameterName(0, "v");
        f_c.setOutputParameterName(0, "c");

        FdmFunctionConvergenceCheck checker(&dis);
	    SerialStaggeredMethod method(checker, 1e-5, 100);
	    AsyncPartitionedSystem apart1;
        apart1.setAlgorithm(method);
        apart1.resizeOutputParameter(3);
        apart1.setOutputParameterName(0, "h");
        apart1.setOutputParameterName(1, "v");
        apart1.setOutputParameterName(2, "c");
        apart1.addProblem(f_head);
        apart1.addProblem(f_v);
        apart1.addProblem(f_c);
        apart1.connectParameters();

	    TimeSteppingController timestepping;
	    timestepping.addTransientSystem(apart1);

	    //const double epsilon = 1.e-3;
	    timestepping.setBeginning(.0);
        timestepping.solve(tim.getEnd());

	    const FdmFunctionScalar* r_f_head = apart1.getOutput<FdmFunctionScalar>(apart1.getOutputParameterID("h"));
	    const FdmCellVectorFunction* r_f_v = apart1.getOutput<FdmCellVectorFunction>(apart1.getOutputParameterID("v"));
	    const IDiscreteVector<double>* vec_h = r_f_head->getNodalValues();
	    //const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();
        const FdmFunctionScalar* r_f_c = apart1.getOutput<FdmFunctionScalar>(apart1.getOutputParameterID("c"));
        //const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();
        const IDiscreteVector<double>* vec_c = r_f_c->getNodalValues();

	    r_f_head->printout();
	    r_f_v->printout();
        r_f_c->printout();

	    std::vector<double> expected;
	    expected.resize(div+1);
        const double p_left = 2.e+6;
        const double p_right = .0;
        const double x_len = 2.0;
	    for (size_t i=0; i<expected.size(); i++) {
            double x = i*h;
            expected[i] = (p_right-p_left) / x_len * x + p_left;
	    }

	    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*vec_h)[0], vec_h->size());

        std::vector<double> expectedC;
        expectedC.resize(div+1);

#ifdef OUTPUT_C
        std::cout << std::endl << "expected C=";
#endif
        for (size_t i=0; i<expectedC.size(); i++) {
            expectedC[i] = analyticalOgataBank(i*h, tim.getEnd(), 1e-5/1.0, 1e-6);
#ifdef OUTPUT_C
            std::cout << expectedC[i] << " ";
#endif
        }
#ifdef OUTPUT_C
        std::cout << std::endl;
#endif
        ASSERT_DOUBLE_ARRAY_EQ(&expectedC[0], &(*vec_c)[0], expectedC.size(), 5e-2);

	} catch (const char* e) {
		std::cout << "***Exception caught! " << e << std::endl;
	}

}

