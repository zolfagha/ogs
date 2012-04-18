
#include <cmath>
#include <gtest/gtest.h>

#include "Base/CodingTools.h"

#include "MathLib/Vector.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "MathLib/Coupling/Algorithm/TransientPartitionedAlgorithm.h"
#include "MeshLib/Tools/MeshGenerator.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/AsyncPartitionedSystem.h"

#include "SolutionLib/Nonlinear.h"

#include "Tests/Geo/Model/Head.h"
#include "Tests/Geo/Model/Velocity.h"
#include "Tests/Geo/Model/Concentration.h"

#include "Tests/ExactSolution/OgataBank.h"

#include "TestUtil.h"


class FemFunctionConvergenceCheck
{
public:
	bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff)
	{
	    for (size_t i=0; i<vars_prev.size(); i++) {
#if 1
	    	if (vars_prev.getName(i).compare("h")==0 || vars_prev.getName(i).compare("c")==0) {
		        const FemNodalFunctionScalar* f_fem_prev = vars_prev.get<FemNodalFunctionScalar>(i);
		        const FemNodalFunctionScalar* f_fem_cur = vars_current.get<FemNodalFunctionScalar>(i);
	    		v_diff = f_fem_cur->norm_diff(*f_fem_prev);
	    	} else if (vars_prev.getName(i).compare("v")==0) {
	    		const FEMIntegrationPointFunctionVector2d* f_fem_prev = vars_prev.get<FEMIntegrationPointFunctionVector2d>(i);
	    		const FEMIntegrationPointFunctionVector2d* f_fem_cur = vars_current.get<FEMIntegrationPointFunctionVector2d>(i);
	    		v_diff = f_fem_cur->norm_diff(*f_fem_prev);
	    	}
#endif
	        if (v_diff>eps) {
	            return false;
	        }
	    }
	    return true;
	}
};


Geo::GWFemProblem* defineGWProblem(DiscreteSystem &dis, Rectangle &_rec, Geo::PorousMedia &pm)
{
    LagrangianFeObjectContainer* _feObjects = new LagrangianFeObjectContainer(*dis.getMesh());
    //equations
    Geo::WeakFormGroundwaterFlow ele_eqs(*_feObjects, pm);
    //IVBV problem
    Geo::GWFemProblem* _problem = new Geo::GWFemProblem(dis, *dis.getMesh(), ele_eqs);
    //BC
    size_t headId = _problem->createField(PolynomialOrder::Linear);
    FemNodalFunctionScalar* _head = _problem->getField(headId);
    Polyline* poly_left = _rec.getLeft();
    Polyline* poly_right = _rec.getRight();
    _problem->setIC(headId, *_head);
    MathLib::SpatialFunctionConstant<double> f1(.0);
    _problem->addDirichletBC(headId, *poly_right, false, f1);
    MathLib::SpatialFunctionConstant<double> f2(-1e-5);
    _problem->addNeumannBC(headId, *poly_left, false, f2);

    return _problem;
}

Geo::MassFemProblem* defineMassTransportProblem(DiscreteSystem &dis, Rectangle &_rec, Geo::PorousMedia &pm, Geo::Compound &comp)
{
    LagrangianFeObjectContainer* _feObjects = new LagrangianFeObjectContainer(*dis.getMesh());
    //equations
    Geo::WeakFormMassTransport ele_eqs(*_feObjects, pm, comp) ;
    //IVBV problem
    Geo::MassFemProblem* _problem = new Geo::MassFemProblem(dis, *dis.getMesh(), ele_eqs);
    //BC
    size_t var_id = _problem->createField(PolynomialOrder::Linear);
    FemNodalFunctionScalar* _conc = _problem->getField(var_id);
    Polyline* poly_left = _rec.getLeft();
    //Polyline* poly_right = _rec.getRight();
    _problem->setIC(var_id, *_conc);
    MathLib::SpatialFunctionConstant<double> f1(1.0);
    _problem->addDirichletBC(var_id, *poly_left, false, f1);

    return _problem;
}


static void getGWExpectedHead(std::vector<double> &expected)
{
    expected.resize(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }
}

typedef Geo::FunctionHead<Linear,CRSLisSolver> MyFunctionHead;
typedef Geo::FunctionVelocity MyFunctionVelocity;
typedef Geo::FunctionConcentration<Linear,CRSLisSolver> MyFunctionConcentration;

#if 0
TEST(Solution, CouplingFem1)
{
	try {
	    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
	    Rectangle* _rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
	    Geo::PorousMedia pm;
	    pm.hydraulic_conductivity = new MathLib::SpatialFunctionConstant<double>(1.e-11);
	    pm.porosity = new MathLib::SpatialFunctionConstant<double>(0.2);
	    DiscreteSystem dis(*msh);
	    Geo::GWFemProblem* pGW = defineGWProblem(dis, *_rec, pm);
        TimeStepFunctionConstant tim(.0, 100.0, 10.0);
        pGW->setTimeSteppingFunction(tim);
	    // options
	    Base::Options options;
	    Base::Options* op_lis = options.addSubGroup("Lis");
	    op_lis->addOption("solver_type", "CG");
	    op_lis->addOption("precon_type", "NONE");
	    op_lis->addOptionAsNum("error_tolerance", 1e-10);
	    op_lis->addOptionAsNum("max_iteration_step", 500);

		MyFunctionHead f_head;
		f_head.define(dis, *pGW, options);
		MyFunctionVelocity f_vel;
		f_vel.define(dis, pm);


	    SerialStaggeredMethod<FemFunctionConvergenceCheck> method(1e-5, 100);
	    AsyncPartitionedSystem apart1(method);
	    apart1.addParameter("h", f_head, MyFunctionHead::Head);
	    apart1.addParameter("v", f_vel, MyFunctionVelocity::Velocity);
	    apart1.connectInput("h", f_vel, MyFunctionVelocity::Head);

	    TimeSteppingController timestepping;
	    timestepping.addTransientSystem(apart1);

	    const double epsilon = 1.e-3;
	    timestepping.setBeginning(.0);
	    timestepping.solve(1.0);

	    const FemNodalFunctionScalar* r_f_head = apart1.getOutput<FemNodalFunctionScalar>(apart1.getParameterID("h"));
	    const FEMIntegrationPointFunctionVector2d* r_f_v = apart1.getOutput<FEMIntegrationPointFunctionVector2d>(apart1.getParameterID("v"));
	    const DiscreteVector<double>* vec_h = r_f_head->getNodalValues();
	    const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();

	    r_f_head->printout();
	    r_f_v->printout();

	    std::vector<double> expected;
	    getGWExpectedHead(expected);

	    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*vec_h)[0], vec_h->size());

	} catch (const char* e) {
		std::cout << "***Exception caught! " << e << std::endl;
	}

}


TEST(Solution, CouplingFem2)
{
	try {
        //space
	    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 20, .0, .0, .0);
	    Rectangle* _rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        //time
        TimeStepFunctionConstant tim(.0, 1e+4, 1e+3);
        //TimeStepFunctionConstant tim(.0, 2e+5, .5e+5);
        //material
	    Geo::PorousMedia pm;
	    pm.hydraulic_conductivity = new MathLib::SpatialFunctionConstant<double>(1.e-11);
	    pm.porosity = new MathLib::SpatialFunctionConstant<double>(1.0);
	    Geo::Compound tracer;
	    tracer.molecular_diffusion = new MathLib::SpatialFunctionConstant<double>(1.e-6);
        //problems
	    DiscreteSystem dis(*msh);
	    Geo::GWFemProblem* pGW = defineGWProblem(dis, *_rec, pm);
	    Geo::MassFemProblem* pMass = defineMassTransportProblem(dis, *_rec, pm, tracer);
        pGW->setTimeSteppingFunction(tim);
        pMass->setTimeSteppingFunction(tim);
	    //options
	    Base::Options optionsGW;
	    Base::Options* op_lis = optionsGW.addSubGroup("Lis");
	    op_lis->addOption("solver_type", "CG");
	    op_lis->addOption("precon_type", "NONE");
	    op_lis->addOptionAsNum("error_tolerance", 1e-10);
	    op_lis->addOptionAsNum("max_iteration_step", 1000);
        Base::Options optionsMT;
        op_lis = optionsMT.addSubGroup("Lis");
        op_lis->addOption("solver_type", "BiCG");
        op_lis->addOption("precon_type", "NONE");
        op_lis->addOptionAsNum("error_tolerance", 1e-10);
        op_lis->addOptionAsNum("max_iteration_step", 1000);
        // unknowns
		MyFunctionHead f_head;
		f_head.define(dis, *pGW, optionsGW);
		MyFunctionVelocity f_vel;
		f_vel.define(dis, pm);
		MyFunctionConcentration f_c;
		f_c.define(dis, *pMass, optionsMT);


	    SerialStaggeredMethod<FemFunctionConvergenceCheck> method(1e-5, 100);
	    AsyncPartitionedSystem apart1(method);
	    apart1.addParameter("h", f_head, MyFunctionHead::Head);
	    apart1.addParameter("v", f_vel, MyFunctionVelocity::Velocity);
	    apart1.addParameter("c", f_c, MyFunctionConcentration::Concentration);
	    apart1.connectInput("h", f_vel, MyFunctionVelocity::Head);
	    apart1.connectInput("v", f_c, MyFunctionConcentration::Velocity);

	    TimeSteppingController timestepping;
	    timestepping.addTransientSystem(apart1);

	    const double epsilon = 1.e-3;
	    timestepping.setBeginning(.0);
	    timestepping.solve(tim.getEnd());

	    const FemNodalFunctionScalar* r_f_head = apart1.getOutput<FemNodalFunctionScalar>(apart1.getParameterID("h"));
	    const FEMIntegrationPointFunctionVector2d* r_f_v = apart1.getOutput<FEMIntegrationPointFunctionVector2d>(apart1.getParameterID("v"));
	    const FemNodalFunctionScalar* r_f_c = apart1.getOutput<FemNodalFunctionScalar>(apart1.getParameterID("c"));
	    const DiscreteVector<double>* vec_h = r_f_head->getNodalValues();
	    const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();
	    const DiscreteVector<double>* vec_c = r_f_c->getNodalValues();

	    //r_f_head->printout();
	    //r_f_v->printout();
#undef OUTPUT_C
//#define OUTPUT_C
#ifdef OUTPUT_C
        r_f_c->printout();
#endif
	    //std::vector<double> expectedHead;
	    //getGWExpectedHead(expectedHead);
        //ASSERT_DOUBLE_ARRAY_EQ(&expectedHead[0], &(*vec_h)[0], vec_h->size());


        std::vector<double> expectedC;
        expectedC.resize(21);

#ifdef OUTPUT_C
        std::cout << std::endl << "expected=";
#endif
        for (size_t i=0; i<expectedC.size(); i++) {
            expectedC[i] = analyticalOgataBank(i*0.1, tim.getEnd(), 1e-5/1.0, 1e-6);
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

#endif
