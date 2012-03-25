
#include <gtest/gtest.h>

#include <vector>

#include "Base/CodingTools.h"

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"

#include "FemLib/FemElementObjectContainer.h"

#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientAssembler/ElementLocalAssembler.h"

#include "SolutionLib/IProblem.h"
#include "SolutionLib/FemProblem.h"
#include "SolutionLib/SingleStepFEM.h"
#include "SolutionLib/Nonlinear.h"

#include "TestUtil.h"

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace FemLib;
using namespace NumLib;
using namespace SolutionLib;
using namespace DiscreteLib;


class GWAssembler: public NumLib::ITimeODEElementAssembler
{
private:
    MathLib::IFunction<double*, double>* _matK;
    FemLib::LagrangianFeObjectContainer* _feObjects;
public:
    GWAssembler(FemLib::LagrangianFeObjectContainer &feObjects, MathLib::IFunction<double*, double> &mat)
    : _matK(&mat), _feObjects(&feObjects)
    {
    };

    //protected:
    void assembly(const NumLib::TimeStep &time, MeshLib::IElement &e, MathLib::DenseLinearEquations::MatrixType &localM, MathLib::DenseLinearEquations::MatrixType &localK,  MathLib::DenseLinearEquations::VectorType &localF)
    {
        IFiniteElement* fe = _feObjects->getFeObject(e);

        //localM = .0;
        //localF.resize(localF.size(), .0);
        fe->integrateDWxDN(_matK, localK);
    }
};

class GWAssemblerNewton //: public NumLib::ITimeODEElementAssembler
{
private:
    MathLib::IFunction<double*, double>* _matK;
    FemLib::LagrangianFeObjectContainer* _feObjects;
public:
    GWAssemblerNewton(FemLib::LagrangianFeObjectContainer &feObjects, MathLib::IFunction<double*, double> &mat)
    : _matK(&mat), _feObjects(&feObjects)
    {
    };

    void assembly(const NumLib::TimeStep &time, MeshLib::IElement &e, MathLib::DenseLinearEquations::VectorType &localX, MathLib::DenseLinearEquations::MatrixType &localJ,  MathLib::DenseLinearEquations::VectorType &localR)
    {
        IFiniteElement* fe = _feObjects->getFeObject(e);
        const size_t n = localR.size();

        const double dt = time.getTimeStepSize();
        const double theta = 1.0;
        MathLib::DenseLinearEquations::VectorType localX0;
        MathLib::DenseLinearEquations::MatrixType localM, localK(n,n);
        MathLib::DenseLinearEquations::VectorType localF;
        //localM = .0;
        //localF.resize(localF.size(), .0);
        fe->integrateDWxDN(_matK, localK);
        
        //localR = (1/dt * localM + theta * localK) * localX - (1/dt * localM + (1-theta) * localK) * localX0 - localF;
        MathLib::DenseLinearEquations::MatrixType tmpM;
        tmpM = localK;
        tmpM *= theta;
        tmpM.axpy(1.0, &localX[0], 0.0, &localR[0]);
        
        localJ = localK;
    }
};

template <
	template <class> class T_NONLINEAR,
	class T_LINEAR_SOLVER
	>
class GWFemTestSystem : public NumLib::ITransientSystem
{
    typedef FemIVBVProblem<GWAssembler> GWFemProblem;
    typedef SingleStepFEM
    		<
    			TimeEulerElementAssembler,
    			T_LINEAR_SOLVER,
    			FemIVBVProblem<GWAssembler>,
    			T_NONLINEAR
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
    void define(DiscreteSystem &dis, MathLib::IFunction<double*, double> &K, Base::Options &option)
    {
        MeshLib::IMesh *msh = dis.getMesh();
        //size_t nnodes = msh->getNumberOfNodes();
        _feObjects = new LagrangianFeObjectContainer(*msh);
        //equations
        GWAssembler ele_eqs(*_feObjects, K) ;
        //IVBV problem
        _problem = new GWFemProblem(dis, *dis.getMesh(), ele_eqs);
        //BC
        size_t headId = _problem->createField(PolynomialOrder::Linear);
        _head = _problem->getField(headId);
        _rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        Polyline* poly_left = _rec->getLeft();
        Polyline* poly_right = _rec->getRight();
        _problem->setIC(headId, *_head);
        MathLib::FunctionConstant<GeoLib::Point, double> f1(.0);
        _problem->addDirichletBC(headId, *poly_right, false, f1);
        MathLib::FunctionConstant<GeoLib::Point, double> f2(-1e-5);
        _problem->addNeumannBC(headId, *poly_left, false, f2);
        //transient
        TimeStepFunctionConstant tim(.0, 100.0, 10.0);
        _problem->setTimeSteppingFunction(tim);
        //solution algorithm
        _solHead = new SolutionForHead(dis, *_problem);
        _solHead->getTimeODEAssembler()->setTheta(1.0);
        T_LINEAR_SOLVER* linear_solver = _solHead->getLinearEquationSolver();
        linear_solver->setOption(option);


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
    MathLib::FunctionConstant<double*, double> K(1.e-11);
    // options
    Base::Options options;
    Base::Options* op_lis = options.addSubGroup("Lis");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    // define problems
    GWFemTestSystem<SolutionLib::Linear, MathLib::CRSLisSolver> gwProblem;
    gwProblem.define(dis, K, options);

    gwProblem.solveTimeStep(TimeStep(1.0, 1.0));

    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemNodalFunctionScalar* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], h->getNodalValues(), h->getNumberOfNodes());

    Base::releaseObject(msh);
}

TEST(Solution, Fem1_Picard)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    MathLib::FunctionConstant<double*, double> K(1.e-11);
    // options
    Base::Options options;
    Base::Options* op_lis = options.addSubGroup("Lis");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    // define problems
    GWFemTestSystem<SolutionLib::Picard, MathLib::CRSLisSolver> gwProblem;
    gwProblem.define(dis, K, options);

    gwProblem.solveTimeStep(TimeStep(1.0, 1.0));

    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemNodalFunctionScalar* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], h->getNodalValues(), h->getNumberOfNodes());

    Base::releaseObject(msh);
}

TEST(Solution, Fem1_Newton)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    MathLib::FunctionConstant<double*, double> K(1.e-11);
    // options
    Base::Options options;
    Base::Options* op_lis = options.addSubGroup("Lis");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    // define problems
    GWFemTestSystem<SolutionLib::NewtonRaphson, MathLib::CRSLisSolver> gwProblem;
    gwProblem.define(dis, K, options);

    gwProblem.solveTimeStep(TimeStep(1.0, 1.0));

    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemNodalFunctionScalar* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], h->getNodalValues(), h->getNumberOfNodes());

    Base::releaseObject(msh);
}

TEST(Solution, Fem2)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    MathLib::FunctionConstant<double*, double> K(1.e-11);
    // options
    Base::Options options;
    Base::Options* op_lis = options.addSubGroup("Lis");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    // define problems
    GWFemTestSystem<SolutionLib::Linear, MathLib::CRSLisSolver> gwProblem;
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

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], h->getNodalValues(), h->getNumberOfNodes());


    Base::releaseObject(msh);
}

