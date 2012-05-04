
#include <gtest/gtest.h>

#include <vector>

#include "Base/CodingTools.h"

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"

#include "FemLib/Tools/FemElementObjectContainer.h"

#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"

#include "SolutionLib/Problem/FemIVBVProblem.h"
#include "SolutionLib/Solution/SingleStepFEM.h"
#include "SolutionLib/Tools/TemplateTransientLinearFEMFunction.h"
#include "SolutionLib/Tools/TemplateTransientResidualFEMFunction.h"
#include "SolutionLib/Tools/TemplateTransientDxFEMFunction.h"

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
    typedef typename T::LocalVectorType LocalVectorType;
    typedef typename T::LocalMatrixType LocalMatrixType;

    GWTimeODEAssembler(FemLib::LagrangianFeObjectContainer &feObjects, MathLib::SpatialFunctionScalar &mat)
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
        	double k;
        	_matK->eval(real_x, k);
        	fe->integrateDWxDN(j, k, localK);
        }
    }
private:
    MathLib::SpatialFunctionScalar* _matK;
    FemLib::LagrangianFeObjectContainer* _feObjects;
};

class GWJacobianAssembler : public NumLib::IElementWiseTransientJacobianLocalAssembler
{
private:
    MathLib::SpatialFunctionScalar* _matK;
    FemLib::LagrangianFeObjectContainer* _feObjects;
public:
    GWJacobianAssembler(FemLib::LagrangianFeObjectContainer &feObjects, MathLib::SpatialFunctionScalar &mat)
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
        MathLib::DenseLinearEquations::VectorType localX0;
        MathLib::DenseLinearEquations::MatrixType localM(n,n), localK(n,n);
        MathLib::DenseLinearEquations::VectorType localF;
        localM = .0;
        localK = .0;
        //localF.resize(localF.size(), .0);
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);
        	double k;
        	_matK->eval(real_x, k);
        	fe->integrateDWxDN(j, k, localK);
        }

        //localR = (1/dt * localM + theta * localK) * localX - (1/dt * localM + (1-theta) * localK) * localX0 - localF;
        MathLib::DenseLinearEquations::MatrixType tmpM(n, n);
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
    typedef FemIVBVProblem
            <
            	GWTimeODEAssembler<ElementWiseTimeEulerEQSLocalAssembler>,
            	GWTimeODEAssembler<ElementWiseTimeEulerResidualLocalAssembler>,
    			GWJacobianAssembler
    		> GWFemProblem;

    typedef TemplateTransientLinearFEMFunction
            <
    			GWFemProblem,
    			typename GWFemProblem::LinearAssemblerType
			> MyLinearFunction;

    typedef TemplateTransientResidualFEMFunction
            <
    			GWFemProblem,
    			typename GWFemProblem::ResidualAssemblerType
			> MyResidualFunction;

    typedef TemplateTransientDxFEMFunction
            <
    			GWFemProblem,
    			typename GWFemProblem::JacobianAssemblerType
			> MyDxFunction;

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
        Base::zeroObject(_rec, _head, _problem);
    }

    virtual ~GWFemTestSystem()
    {
        Base::releaseObject(_rec, _head, _problem);
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
    void define(DiscreteSystem &dis, MathLib::SpatialFunctionScalar &K, Base::Options &option)
    {
        MeshLib::IMesh *msh = dis.getMesh();
        //size_t nnodes = msh->getNumberOfNodes();
        _feObjects = new LagrangianFeObjectContainer(*msh);
        //equations
        GWFemProblem::LinearAssemblerType* local_linear = new GWFemProblem::LinearAssemblerType(*_feObjects, K);
        GWFemProblem::ResidualAssemblerType* local_r = new GWFemProblem::ResidualAssemblerType(*_feObjects, K);
        GWFemProblem::JacobianAssemblerType* local_J = new GWFemProblem::JacobianAssemblerType(*_feObjects, K);
        //IVBV problem
        _problem = new GWFemProblem(dis, local_linear, local_r, local_J);
        //BC
        size_t headId = _problem->createField(PolynomialOrder::Linear);
        _head = _problem->getField(headId);
        _rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        Polyline* poly_left = _rec->getLeft();
        Polyline* poly_right = _rec->getRight();
        _problem->setIC(headId, *_head);
        MathLib::SpatialFunctionConstant<double> f1(.0);
        _problem->addDirichletBC(headId, *poly_right, false, f1);
        MathLib::SpatialFunctionConstant<double> f2(-1e-5);
        _problem->addNeumannBC(headId, *poly_left, false, f2);
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
    MathLib::SpatialFunctionConstant<double> K(1.e-11);
    // options
    Base::Options options;
    Base::Options* op_lis = options.addSubGroup("Lis");
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

    Base::releaseObject(msh);
}

TEST(Solution, Fem1_Picard)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    MathLib::SpatialFunctionConstant<double> K(1.e-11);
    // options
    Base::Options options;
    Base::Options* op_lis = options.addSubGroup("Lis");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    Base::Options* op_nl = options.addSubGroup("Nonlinear");
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

    Base::releaseObject(msh);
}

TEST(Solution, Fem1_Newton)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    MathLib::SpatialFunctionConstant<double> K(1.e-11);
    // options
    Base::Options options;
    Base::Options* op_lis = options.addSubGroup("Lis");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    Base::Options* op_nl = options.addSubGroup("Nonlinear");
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

    Base::releaseObject(msh);
}

TEST(Solution, Fem2)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    MathLib::SpatialFunctionConstant<double> K(1.e-11);
    // options
    Base::Options options;
    Base::Options* op_lis = options.addSubGroup("Lis");
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


    Base::releaseObject(msh);
}

