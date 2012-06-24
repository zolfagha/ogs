
#include <gtest/gtest.h>

#include <vector>

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"

#include "FemLib/Tools/FemElementObjectContainer.h"

#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"

#include "SolutionLib/FemProblem/FemDirichletBC.h"
#include "SolutionLib/FemProblem/FemNeumannBC.h"
#include "SolutionLib/FemProblem/FemIVBVProblem.h"
#include "SolutionLib/Solution/SingleStepFEM.h"

#include "TestUtil.h"

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace FemLib;
using namespace NumLib;
using namespace SolutionLib;
using namespace DiscreteLib;

template <class T>
class GWTimeODEAssembler: public T
{
public:
    typedef NumLib::LocalVector LocalVectorType;
    typedef NumLib::LocalMatrix LocalMatrixType;

    GWTimeODEAssembler(FemLib::LagrangianFeObjectContainer &feObjects, NumLib::ITXFunction &mat)
    : _matK(&mat), _feObjects(&feObjects)
    {
    };

    virtual ~GWTimeODEAssembler() {};

protected:
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/, LocalMatrixType &/*localM*/, LocalMatrixType &localK, LocalVectorType &/*localF*/)
    {
        IFiniteElement* fe = _feObjects->getFeObject(e);

        //localM = .0;
        //localF.resize(localF.size(), .0);
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);
        	LocalMatrixType k;
        	_matK->eval(real_x, k);
        	fe->integrateDWxDN(j, k, localK);
        }
    }
private:
    NumLib::ITXFunction* _matK;
    FemLib::LagrangianFeObjectContainer* _feObjects;
};

class GWJacobianAssembler : public NumLib::IElementWiseTransientJacobianLocalAssembler
{
private:
	NumLib::ITXFunction* _matK;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    typedef NumLib::LocalMatrix LocalMatrixType;
    typedef NumLib::LocalVector LocalVectorType;
public:
    GWJacobianAssembler(FemLib::LagrangianFeObjectContainer &feObjects, NumLib::ITXFunction &mat)
    : _matK(&mat), _feObjects(&feObjects)
    {
    };

    /// assemble a local Jacobian matrix for the given element
    /// @param time			time step
    /// @param e			element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param local_J		local Jacobian
    virtual void assembly(const TimeStep &/*time*/,  MeshLib::IElement &e, const LocalVectorType &local_u_n1, const LocalVectorType &/*local_u_n*/, LocalMatrixType &local_J)
    {
        IFiniteElement* fe = _feObjects->getFeObject(e);
        const size_t n = local_u_n1.size();

        //const double dt = time.getTimeStepSize();
        const double theta = 1.0;
        LocalVectorType localX0;
        LocalMatrixType localM(n,n), localK(n,n);
        LocalVectorType localF;
        localM *= .0;
        localK *= .0;
        //localF.resize(localF.size(), .0);
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);
        	NumLib::LocalMatrix k;
        	_matK->eval(real_x, k);
        	fe->integrateDWxDN(j, k, localK);
        }

        //localR = (1/dt * localM + theta * localK) * localX - (1/dt * localM + (1-theta) * localK) * localX0 - localF;
        LocalMatrixType tmpM(n, n);
        tmpM = localK;
        tmpM *= theta;
        //tmpM.axpy(1.0, &localX[0], 0.0, &localR[0]);

        local_J = localK;
    }
};

template <
	class T_LINEAR_SOLVER
	>
class GWFemTestSystem : public NumLib::ITransientSystem
{
	typedef TemplateFemEquation
            <
            	GWTimeODEAssembler<ElementWiseTimeEulerEQSLocalAssembler>,
            	GWTimeODEAssembler<ElementWiseTimeEulerResidualLocalAssembler>,
    			GWJacobianAssembler
    		> GWFemEquation;

    typedef FemIVBVProblem<GWFemEquation> GWFemProblem;

    typedef GWFemEquation::LinearEQSType MyLinearFunction;
    typedef GWFemEquation::ResidualEQSType MyResidualFunction;
    typedef GWFemEquation::DxEQSType MyDxFunction;

    typedef TemplateDiscreteNonlinearSolver
            <
                MyLinearFunction, 
                MyResidualFunction, 
                MyDxFunction
            > MyNonlinearFunction;

    typedef SingleStepFEM
    		<
    			GWFemProblem,
    			T_LINEAR_SOLVER
    		> SolutionForHead;

public:
    GWFemTestSystem()
    {
        BaseLib::zeroObject(_rec, _head, _problem);
    }

    virtual ~GWFemTestSystem()
    {
        BaseLib::releaseObject(_rec, _head, _problem);
    }

    double suggestNext(const TimeStep &time_current)
    {
        return _solHead->suggestNext(time_current);
    }
    bool isAwake(const TimeStep &time)
    {
        return _solHead->isAwake(time);
    }

    int solveTimeStep(const TimeStep &time)
    {
        _solHead->solveTimeStep(time);
        _head = _solHead->getCurrentSolution(0);
        return 0;
    }

    void accept(const TimeStep &time) 
    {
        _solHead->accept(time);
    };

    //#Define a problem
    void define(DiscreteSystem &dis, NumLib::ITXFunction &K, BaseLib::Options &option)
    {
        MeshLib::IMesh *msh = dis.getMesh();
        //size_t nnodes = msh->getNumberOfNodes();
        _feObjects = new LagrangianFeObjectContainer(*msh);
        //equations
        GWFemEquation::LinearAssemblerType* local_linear = new GWFemEquation::LinearAssemblerType(*_feObjects, K);
        GWFemEquation::ResidualAssemblerType* local_r = new GWFemEquation::ResidualAssemblerType(*_feObjects, K);
        GWFemEquation::JacobianAssemblerType* local_J = new GWFemEquation::JacobianAssemblerType(*_feObjects, K);
        GWFemEquation* eqs = new GWFemEquation(local_linear, local_r, local_J);
        //IVBV problem
        _problem = new GWFemProblem(&dis);
        _problem->setEquation(eqs);
        // var
        FemVariable* _head = _problem->addVariable("head");
        //IC
        FemNodalFunctionScalar* h0 = new FemNodalFunctionScalar(dis, PolynomialOrder::Linear, .0);
        _head->setIC(h0);
        //BC
        _rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        Polyline* poly_left = _rec->getLeft();
        Polyline* poly_right = _rec->getRight();
        //_head->setIC();
        _head->addDirichletBC(new SolutionLib::FemDirichletBC(msh, poly_right,  new NumLib::TXFunctionConstant(.0)));
        _head->addNeumannBC(new SolutionLib::FemNeumannBC(msh, _feObjects, poly_left, new NumLib::TXFunctionConstant(-1e-5)));
        //transient
        TimeStepFunctionConstant tim(.0, 100.0, 10.0);
        _problem->setTimeSteppingFunction(tim);
        //solution algorithm
        _solHead = new SolutionForHead(&dis, _problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        typename SolutionForHead::LinearSolverType* linear_solver = _solHead->getLinearEquationSolver();
        linear_solver->setOption(option);
        typename SolutionForHead::NonlinearSolverType* nonlinear_solver = _solHead->getNonlinearSolver();
        nonlinear_solver->setOption(option);


        //vel = new FEMIntegrationPointFunctionVector2d(msh);
    }

    FemNodalFunctionScalar* getCurrentHead()
    {
        return _head;
    }

private:
    GWFemProblem* _problem;
    SolutionForHead* _solHead; 
    Rectangle *_rec;
    FemNodalFunctionScalar *_head;
    LagrangianFeObjectContainer* _feObjects;

    DISALLOW_COPY_AND_ASSIGN(GWFemTestSystem);
};





TEST(Solution, Fem1_Linear)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    NumLib::TXFunctionConstant K(1.e-11);
    // options
    BaseLib::Options options;
    BaseLib::Options* op_lis = options.addSubGroup("Lis");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    // define problems
    GWFemTestSystem<MathLib::CRSLisSolver> gwProblem;
    gwProblem.define(dis, K, options);

    gwProblem.solveTimeStep(TimeStep(1.0, 1.0));

    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemNodalFunctionScalar* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h->getNodalValues())[0], h->getNumberOfNodes());

    BaseLib::releaseObject(msh);
}

TEST(Solution, Fem1_Picard)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    NumLib::TXFunctionConstant K(1.e-11);
    // options
    BaseLib::Options options;
    BaseLib::Options* op_lis = options.addSubGroup("Lis");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    BaseLib::Options* op_nl = options.addSubGroup("Nonlinear");
    op_nl->addOption("solver_type", "Picard");
    op_nl->addOptionAsNum("error_tolerance", 1e-6);
    op_nl->addOptionAsNum("max_iteration_step", 500);
    // define problems
    GWFemTestSystem<MathLib::CRSLisSolver> gwProblem;
    gwProblem.define(dis, K, options);

    gwProblem.solveTimeStep(TimeStep(1.0, 1.0));

    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemNodalFunctionScalar* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h->getNodalValues())[0], h->getNumberOfNodes());

    BaseLib::releaseObject(msh);
}

TEST(Solution, Fem1_Newton)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    NumLib::TXFunctionConstant K(1.e-11);
    // options
    BaseLib::Options options;
    BaseLib::Options* op_lis = options.addSubGroup("Lis");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    BaseLib::Options* op_nl = options.addSubGroup("Nonlinear");
    op_nl->addOption("solver_type", "Newton");
    op_nl->addOptionAsNum("error_tolerance", 1e-6);
    op_nl->addOptionAsNum("max_iteration_step", 10);
    // define problems
    GWFemTestSystem<MathLib::CRSLisSolver> gwProblem;
    gwProblem.define(dis, K, options);

    gwProblem.solveTimeStep(TimeStep(1.0, 1.0));

    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemNodalFunctionScalar* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h->getNodalValues())[0], h->getNumberOfNodes());

    BaseLib::releaseObject(msh);
}

TEST(Solution, Fem2)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    NumLib::TXFunctionConstant K(1.e-11);
    // options
    BaseLib::Options options;
    BaseLib::Options* op_lis = options.addSubGroup("Lis");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    // define problems
    GWFemTestSystem<MathLib::CRSLisSolver> gwProblem;
    gwProblem.define(dis, K, options);

    // start time stepping
    TimeSteppingController timeStepping;
    timeStepping.addTransientSystem(gwProblem);
    timeStepping.setBeginning(.0);
    timeStepping.solve(30.);


    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemNodalFunctionScalar* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h->getNodalValues())[0], h->getNumberOfNodes());


    BaseLib::releaseObject(msh);
}

