
#include <gtest/gtest.h>

#include "Base/CodingTools.h"

#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "MeshLib/Tools/MeshGenerator.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/Coupling/Algorithm/PartitionedAlgorithm.h"
#include "NumLib/TransientCoupling/TransientPartitionedMethod.h"
#include "NumLib/TransientCoupling/AsyncPartSolution.h"

#include "SolutionLib/Nonlinear.h"

#include "Head.h"
#include "Velocity.h"
#include "Concentration.h"

class MyConvergenceCheck
{
	typedef VariableContainer MyNamedVariableContainer;
public:
	bool isConverged(MyNamedVariableContainer& vars_prev, MyNamedVariableContainer& vars_current, double eps, double &v_diff)
	{
	    for (size_t i=0; i<vars_prev.size(); i++) {
	    	if (vars_prev.getName(i).compare("h")==0) {
		        FemNodalFunctionScalar* f_fem_prev = vars_prev.get<FemNodalFunctionScalar>(i);
		        FemNodalFunctionScalar* f_fem_cur = vars_current.get<FemNodalFunctionScalar>(i);
	  	        DiscreteLib::DiscreteVector<double>* vec_prev = f_fem_prev->getNodalValuesAsStdVec();
	  	        DiscreteLib::DiscreteVector<double>* vec_cur = f_fem_cur->getNodalValuesAsStdVec();
	  	        DiscreteLib::DiscreteVector<double> vec_diff(vec_prev->size());
	  	        vec_diff = *vec_cur;
	  	        vec_diff -= *vec_prev;
	  			v_diff = norm_max(vec_diff, vec_diff.size());

	    	} else {
	    		FEMIntegrationPointFunctionVector2d* f_fem_prev = vars_prev.get<FEMIntegrationPointFunctionVector2d>(i);
	    		FEMIntegrationPointFunctionVector2d* f_fem_cur = vars_current.get<FEMIntegrationPointFunctionVector2d>(i);

	    	}
	        if (v_diff>eps) {
	            return false;
	        }
	    }
	    return true;
	}
};


GWFemProblem* defineGWProblem(DiscreteSystem &dis, Rectangle &_rec, PorousMedia &pm)
{
    LagrangianFeObjectContainer* _feObjects = new LagrangianFeObjectContainer(*dis.getMesh());
    //equations
    GWAssembler ele_eqs(*_feObjects, *pm.K) ;
    //IVBV problem
    GWFemProblem* _problem = new GWFemProblem(dis, *dis.getMesh(), ele_eqs);
    //BC
    size_t headId = _problem->createField(PolynomialOrder::Linear);
    FemNodalFunctionScalar* _head = _problem->getField(headId);
    Polyline* poly_left = _rec.getLeft();
    Polyline* poly_right = _rec.getRight();
    _problem->setIC(headId, *_head);
    MathLib::FunctionConstant<GeoLib::Point, double> f1(.0);
    _problem->addDirichletBC(headId, *poly_right, false, f1);
    MathLib::FunctionConstant<GeoLib::Point, double> f2(-1e-5);
    _problem->addNeumannBC(headId, *poly_left, false, f2);
    //transient
    TimeStepFunctionConstant tim(.0, 100.0, 10.0);
    _problem->setTimeSteppingFunction(tim);

    return _problem;
}


typedef FunctionHead<Linear,CRSLisSolver> MyFunctionHead;
typedef FunctionVelocity MyFunctionVelocity;
typedef FunctionConcentration MyFunctionConcentration;


TEST(Solution, CouplingFem1)
{
	try {
	    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
	    Rectangle* _rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
	    PorousMedia pm;
	    pm.K = new MathLib::FunctionConstant<double*, double>(1.e-11);
	    pm.porosity = new MathLib::FunctionConstant<double*, double>(0.2);
	    DiscreteSystem dis(*msh);
	    GWFemProblem* pGW = defineGWProblem(dis, *_rec, pm);
	    // options
	    Base::Options options;
	    Base::Options* op_lis = options.addSubGroup("Lis");
	    op_lis->addOption("solver_type", "CG");
	    op_lis->addOption("precon_type", "NONE");
	    op_lis->addOptionAsNum("error_tolerance", 1e-10);
	    op_lis->addOptionAsNum("max_iteration_step", 500);

		MyFunctionHead f_head;
		f_head.define(dis, *pGW, options);
		f_head.setParameter(MyFunctionHead::Head, pGW->getIC(0));
		MyFunctionVelocity f_vel;
		f_vel.define(dis, pm);
		f_vel.setParameter(MyFunctionVelocity::Velocity, new FemLib::FEMIntegrationPointFunctionVector2d(dis, *msh));
		//MyFunctionConcentration f_c;


	    SerialStaggeredMethod<MyConvergenceCheck> method(1e-5, 100);
	    AsyncPartitionedSystem apart1(method);
	    apart1.addParameter("h", f_head, MyFunctionHead::Head);
	    apart1.addParameter("v", f_vel, MyFunctionVelocity::Velocity);
	    //apart1.addParameter("c", f_c, MyFunctionConcentration::Concentration);
	    apart1.connectInput("h", f_vel, MyFunctionVelocity::Head);
	    //apart1.connectInput("v", f_c, MyFunctionConcentration::Velocity);

	    TimeSteppingController timestepping;
	    timestepping.addTransientSystem(apart1);

	    const double epsilon = 1.e-3;
	    timestepping.setBeginning(.0);
	    timestepping.solve(1.0);
	//    double v1, v2, v3;
	//    apart1.getParameter(apart1.getParameterID("a"))->eval(0, v1);
	//    apart1.getParameter(apart1.getParameterID("b"))->eval(0, v2);
	//    apart1.getParameter(apart1.getParameterID("c"))->eval(0, v3);
	//    ASSERT_NEAR(1., v1, epsilon);
	//    ASSERT_NEAR(2., v2, epsilon);
	//    ASSERT_NEAR(3., v3, epsilon);

	} catch (const char* e) {
		std::cout << "***Exception caught! " << e << std::endl;
	}

}

