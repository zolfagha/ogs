
#include <vector>

#include <gtest/gtest.h>

#include "Base/CodingTools.h"
#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "MathLib/Coupling/Algorithm/IConvergenceCheck.h"
#include "MathLib/Coupling/Algorithm/SerialStaggeredMethod.h"
#include "GeoLib/Core/Polyline.h"
#include "GeoLib/Shape/Line.h"
#include "GeoLib/Shape/Rectangle.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "MeshLib/Core/IStencil.h"
#include "MeshLib/Tools/MeshGenerator.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientCoupling/AsyncPartitionedSystem.h"
#include "SolutionLib/Problem/FdmIVBVProblem.h"
#include "SolutionLib/Solution/SingleStepFDM.h"

#include "TestUtil.h"

#include "Tests/Geo/Material/PorousMedia.h"

class FdmGw1DLocalAssembler: public NumLib::IStencilWiseTransientLinearEQSLocalAssembler
{
private:
	Geo::PorousMedia* _pm;
	double _h;
	double _h2;
public:
	FdmGw1DLocalAssembler(double h, Geo::PorousMedia &pm)
	: _pm(&pm), _h(h), _h2(h*h)
	{
	};

    virtual void assembly(const NumLib::TimeStep &time,  MeshLib::IStencil &s, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalEquationType &eqs)
    {
    	const double dt = time.getTimeStepSize();
    	const size_t center_point_id = s.getCentralNodeID();
    	//double storage = .0;
    	//_pm->storage->eval(0, storage);
		double k = .0;
    	_pm->hydraulic_conductivity->eval(0, k);

    	if (center_point_id==0) {
            eqs.addA(0, 0, k/_h);
            const std::vector<size_t> &neighbor_points = s.getSurroundingNodes();
            for(size_t i=0; i<neighbor_points.size(); i++)
            {
                eqs.addA(0, i+1, -k/_h);
            }
    	} else if (center_point_id == 10) {
            eqs.addA(0, 0, 1.0);
    	} else {
        	eqs.addA(0, 0, 2.0*k/_h2);
    		const std::vector<size_t> &neighbor_points = s.getSurroundingNodes();

            for(size_t i=0; i<neighbor_points.size(); i++)
            {
    		  eqs.addA(0, i+1, -k/_h2);
            }
    	}

    }
};

class FdmGwResidualLocalAssembler : public NumLib::IStencilWiseTransientResidualLocalAssembler
{
public:
    virtual ~FdmGwResidualLocalAssembler() {};

    /// assemble a local residual for the given element
    /// @param time			time step
    /// @param e			element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param local_r		local residual
    virtual void assembly(const NumLib::TimeStep &time,  MeshLib::IStencil &e, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalVectorType &local_r)
    {

    }
};

class FdmGwJacobianLocalAssembler: public NumLib::IStencilWiseTransientJacobianLocalAssembler
{
private:
public:
	FdmGwJacobianLocalAssembler()
	{
	};

	void assembly(const NumLib::TimeStep &/*time*/, MeshLib::IStencil &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/,  LocalMatrixType &localJ)
	{
	}
};

typedef SolutionLib::FdmIVBVProblem
		<
		FdmGw1DLocalAssembler,
		FdmGwResidualLocalAssembler,
		FdmGwJacobianLocalAssembler
		> GWFdmProblem;



template <
	class T_LINEAR_SOLVER
	>
class FunctionHead : public NumLib::TemplateTransientMonolithicSystem
{
    enum Out { Head=0 };
public:
    typedef SolutionLib::SingleStepFDM
    		<
    			GWFdmProblem,
    			T_LINEAR_SOLVER
    		> SolutionForHead;

    FunctionHead()
    {
        TemplateTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteLib::DiscreteSystem* dis, GWFdmProblem* problem, Base::Options &option)
    {
        //solution algorithm
        _solHead = new SolutionForHead(dis, problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        typename SolutionForHead::LinearSolverType* linear_solver = _solHead->getLinearEquationSolver();
        linear_solver->setOption(option);
        _solHead->getNonlinearSolver()->setOption(option);
        this->setOutput(Head, problem->getIC(0));
    }

    int solveTimeStep(const NumLib::TimeStep &time)
    {
        _solHead->solveTimeStep(time);
        setOutput(Head, _solHead->getCurrentSolution(0));
        return 0;
    }

    double suggestNext(const NumLib::TimeStep &time_current) { return _solHead->suggestNext(time_current); }

    bool isAwake(const NumLib::TimeStep &time) { return _solHead->isAwake(time);  }

    void accept(const NumLib::TimeStep &time)
    {
        _solHead->accept(time);

        //std::cout << "Head=" << std::endl;
        //_solHead->getCurrentSolution(0)->printout();
    };

private:
    GWFdmProblem* _problem;
    SolutionForHead* _solHead;
    GeoLib::Rectangle *_rec;

    DISALLOW_COPY_AND_ASSIGN(FunctionHead);
};


using namespace MathLib;
using namespace GeoLib;
using namespace DiscreteLib;
using namespace NumLib;
using namespace FdmLib;

typedef FunctionHead<MathLib::CRSLisSolver> MyFunctionHead;

class FdmFunctionConvergenceCheck : public IConvergenceCheck
{
public:
	bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff)
	{
	    for (size_t i=0; i<vars_prev.size(); i++) {
	    	if (vars_prev.getName(i).compare("h")==0 || vars_prev.getName(i).compare("c")==0) {
		        const FdmFunctionScalar* f_fem_prev = vars_prev.get<FdmFunctionScalar>(i);
		        const FdmFunctionScalar* f_fem_cur = vars_current.get<FdmFunctionScalar>(i);
	    		v_diff = f_fem_cur->norm_diff(*f_fem_prev);
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


GWFdmProblem* defineGWProblem4FDM(DiscreteSystem &dis, double h, Line &line, Geo::PorousMedia &pm)
{
    //equations
    GWFdmProblem::LinearAssemblerType* linear_assembler = new GWFdmProblem::LinearAssemblerType(h, pm);
    GWFdmProblem::ResidualAssemblerType* r_assembler = new GWFdmProblem::ResidualAssemblerType();
    GWFdmProblem::JacobianAssemblerType* j_eqs = new GWFdmProblem::JacobianAssemblerType();
    //IVBV problem
    GWFdmProblem* _problem = new GWFdmProblem(dis, *dis.getMesh(), linear_assembler, r_assembler, j_eqs);
    //BC
    size_t headId = _problem->createField();
    FdmFunctionScalar* _head = _problem->getField(headId);
    _problem->setIC(headId, *_head);
    MathLib::SpatialFunctionConstant<double> f1(.0);
    _problem->addDirichletBC(headId, *line.getPoint2(), false, f1);
    MathLib::SpatialFunctionConstant<double> f2(-1e-5);
    _problem->addNeumannBC(headId, *line.getPoint1(), false, f2);

    return _problem;
}

TEST(Fdm, fdm1)
{
	try {
		const double len = 2.0;
		const double div = 10;
		const double h = len / div;
	    MeshLib::IMesh *msh = MeshLib::MeshGenerator::generateLineMesh(len, div, .0, .0, .0);
	    GeoLib::Line line(Point(0.0, .0, .0), Point(len, .0, .0));
	    Geo::PorousMedia pm;
	    pm.hydraulic_conductivity = new MathLib::SpatialFunctionConstant<double>(1.e-11);
	    DiscreteSystem dis(*msh);
	    GWFdmProblem* pGW = defineGWProblem4FDM(dis, h, line, pm);
        TimeStepFunctionConstant tim(.0, 10.0, 10.0);
        pGW->setTimeSteppingFunction(tim);
	    // options
	    Base::Options options;
	    Base::Options* op_lis = options.addSubGroup("Lis");
	    op_lis->addOption("solver_type", "BiCG");
	    op_lis->addOption("precon_type", "NONE");
	    op_lis->addOptionAsNum("error_tolerance", 1e-10);
	    op_lis->addOptionAsNum("max_iteration_step", 500);

		MyFunctionHead f_head;
		f_head.define(&dis, pGW, options);
        f_head.setOutputParameterName(0, "h");

        FdmFunctionConvergenceCheck checker;
	    SerialStaggeredMethod method(checker, 1e-5, 100);
	    AsyncPartitionedSystem apart1;
        apart1.setAlgorithm(method);
        apart1.resizeOutputParameter(1);
//        apart1.resizeOutputParameter(2);
        apart1.setOutputParameterName(0, "h");
//        apart1.setOutputParameterName(1, "v");
        apart1.addProblem(f_head);
        apart1.connectParameters();

	    TimeSteppingController timestepping;
	    timestepping.addTransientSystem(apart1);

	    //const double epsilon = 1.e-3;
	    timestepping.setBeginning(.0);
	    timestepping.solve(1.0);

	    const FdmFunctionScalar* r_f_head = apart1.getOutput<FdmFunctionScalar>(apart1.getOutputParameterID("h"));
	    //const FEMIntegrationPointFunctionVector2d* r_f_v = apart1.getOutput<FEMIntegrationPointFunctionVector2d>(apart1.getOutputParameterID("v"));
	    const DiscreteVector<double>* vec_h = r_f_head->getNodalValues();
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
