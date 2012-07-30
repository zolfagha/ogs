
#pragma once

#include "MathLib/LinAlg/VectorNorms.h"
#include "GeoLib/Shape/Rectangle.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/Function/Function.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "SolutionLib/FemProblem/FemIVBVProblem.h"
#include "SolutionLib/Solution/SingleStepFEM.h"

#include "Tests/Geo/Equation/FemGroundwaterFlow.h"
#include "Tests/Geo/Material/PorousMedia.h"

using namespace GeoLib;
using namespace MathLib;
using namespace FemLib;
using namespace NumLib;
using namespace MeshLib;
using namespace SolutionLib;
using namespace DiscreteLib;

namespace Geo
{

typedef TemplateFemEquation<
        Geo::GroundwaterFlowTimeODELocalAssembler<ElementWiseTimeEulerEQSLocalAssembler>,
        Geo::GroundwaterFlowTimeODELocalAssembler<ElementWiseTimeEulerResidualLocalAssembler>,
        Geo::GroundwaterFlowJacobianLocalAssembler
        >
        GWFemEquation;

typedef FemIVBVProblem< GWFemEquation > GWFemProblem;



template <
    class T_LINEAR_SOLVER
    >
class FunctionHead : public AbstractTransientMonolithicSystem
{
    enum Out { Head=0 };
public:
    typedef SingleStepFEM
            <
                GWFemProblem,
                T_LINEAR_SOLVER
            > SolutionForHead;

    FunctionHead() 
    {
        AbstractTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteSystem* dis, GWFemProblem* problem, BaseLib::Options &option)
    {
        //MeshLib::IMesh *msh = dis.getMesh();
        //size_t nnodes = msh->getNumberOfNodes();
        //solution algorithm
        _solHead = new SolutionForHead(dis, problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        typename SolutionForHead::LinearSolverType* linear_solver = _solHead->getLinearEquationSolver();
        linear_solver->setOption(option);
        _solHead->getNonlinearSolver()->setOption(option);
        this->setOutput(Head, problem->getVariable(0)->getIC());
    }

    int solveTimeStep(const TimeStep &time)
    {
        _solHead->solveTimeStep(time);
        setOutput(Head, _solHead->getCurrentSolution(0));
        return 0;
    }

    double suggestNext(const TimeStep &time_current) { return _solHead->suggestNext(time_current); }

    bool isAwake(const TimeStep &time) { return _solHead->isAwake(time);  }

    void accept(const TimeStep &time)
    {
        _solHead->accept(time);

        //std::cout << "Head=" << std::endl;
        //_solHead->getCurrentSolution(0)->printout();
    };

private:
    GWFemProblem* _problem;
    SolutionForHead* _solHead;
    GeoLib::Rectangle *_rec;
    //FemNodalFunctionScalar *_head;
    LagrangianFeObjectContainer* _feObjects;

    DISALLOW_COPY_AND_ASSIGN(FunctionHead);
};


} //end


