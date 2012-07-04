
#pragma once

#if 1
#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/Function/Function.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "SolutionLib/FemProblem/FemIVBVProblem.h"
#include "SolutionLib/Solution/SingleStepFEM.h"
#include "MaterialLib/PorousMedia.h"

#include "FemGroundwaterFlow.h"

#include "ProcessLib/ProcessBuilder.h"

using namespace GeoLib;
using namespace MathLib;
using namespace FemLib;
using namespace NumLib;
using namespace MeshLib;
using namespace SolutionLib;
using namespace DiscreteLib;

namespace Geo
{

//--------------------------------------------------------------------------------------------------
// Equation definition
//--------------------------------------------------------------------------------------------------
typedef TemplateFemEquation<
		Geo::GroundwaterFlowTimeODELocalAssembler<ElementWiseTimeEulerEQSLocalAssembler>,
		Geo::GroundwaterFlowTimeODELocalAssembler<ElementWiseTimeEulerResidualLocalAssembler>,
		Geo::GroundwaterFlowJacobianLocalAssembler
		>
		GWFemEquation;


//--------------------------------------------------------------------------------------------------
// IVBV problem definition
//--------------------------------------------------------------------------------------------------
typedef FemIVBVProblem< GWFemEquation > GWFemProblem;



//--------------------------------------------------------------------------------------------------
// Function definition
//--------------------------------------------------------------------------------------------------

template <
	class T_PROBLEM,
	class T_LINEAR_SOLVER
	>
class AbstractTransientFemFunction : public TemplateTransientMonolithicSystem
{
public:
	typedef T_PROBLEM MyProblemType;
    typedef SingleStepFEM
    		<
    			MyProblemType,
    			T_LINEAR_SOLVER
    		> MySolutionType;

//    void define(DiscreteSystem* dis, MyProblemType* problem, BaseLib::Options &option)
//    {
//        _solution = new MySolutionType(dis, problem);
//        typename MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
//        linear_solver->setOption(option);
//        _solution->getNonlinearSolver()->setOption(option);
//        this->setOutput(Head, problem->getVariable(0)->getIC());
//    }
//
//    int solveTimeStep(const TimeStep &time)
//    {
//        _solution->solveTimeStep(time);
//        setOutput(Head, _solution->getCurrentSolution(0));
//        return 0;
//    }

    double suggestNext(const TimeStep &time_current) { return _solution->suggestNext(time_current); }

    bool isAwake(const TimeStep &time) { return _solution->isAwake(time);  }

    void accept(const TimeStep &time)
    {
        _solution->accept(time);
    };

private:
    MyProblemType* _problem;
    MySolutionType* _solution;
    LagrangianFeObjectContainer* _feObjects;
};


class FunctionHead : public TemplateTransientMonolithicSystem
{
    enum Out { Head=0 };
public:
    typedef SingleStepFEM
    		<
    			GWFemProblem,
    			MathLib::CRSLisSolver
    		> SolutionForHead;

    FunctionHead() 
    {
        TemplateTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteSystem* dis, GWFemProblem* problem, BaseLib::Options &option)
    {
        //solution algorithm
        _solHead = new SolutionForHead(dis, problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        SolutionForHead::LinearSolverType* linear_solver = _solHead->getLinearEquationSolver();
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
    LagrangianFeObjectContainer* _feObjects;

    DISALLOW_COPY_AND_ASSIGN(FunctionHead);

public:
    static ProcessLib::ProcessInfo* const _pcs_info;
};


} //end

OGS_PROCESS(GROUNDWATER_FLOW, Geo::FunctionHead);

#endif
