
#pragma once

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "GeoLib/Shape/Rectangle.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientAssembler/ElementLocalAssembler.h"
#include "NumLib/TransientAssembler/TimeEulerElementLocalAssembler.h"
#include "SolutionLib/FemProblem.h"
#include "SolutionLib/SingleStepFEM.h"

#include "Assembler.h"
#include "PorousMedia.h"

using namespace GeoLib;
using namespace MathLib;
using namespace FemLib;
using namespace NumLib;
using namespace MeshLib;
using namespace SolutionLib;
using namespace DiscreteLib;

typedef FemIVBVProblem<GWAssembler> GWFemProblem;

template <
	template <class> class T_NONLINEAR,
	class T_LINEAR_SOLVER
	>
class FunctionHead : public TemplateTransientMonolithicSystem<0,1>
{
public:
    typedef TemplateTransientLinearFEMFunction<
    			FemIVBVProblem,
    			TimeEulerElementLocalAssembler,
    			T_LINEAR_SOLVER,
    			GWAssembler
			> MyLinearFunction;

    typedef T_NONLINEAR<MyLinearFunction> MyNonlinearFunction;

    typedef SingleStepFEM
    		<
    			FemIVBVProblem<GWAssembler>,
    			//TimeEulerElementLocalAssembler,
    			MyLinearFunction,
    			MyNonlinearFunction,
    			T_LINEAR_SOLVER
    		> SolutionForHead;

    enum Parameters { Head=0 };

    FunctionHead() {};

    void define(DiscreteSystem &dis, GWFemProblem &problem, Base::Options &option)
    {
        //MeshLib::IMesh *msh = dis.getMesh();
        //size_t nnodes = msh->getNumberOfNodes();
        //solution algorithm
        _solHead = new SolutionForHead(dis, problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        T_LINEAR_SOLVER* linear_solver = _solHead->getLinearEquationSolver();
        linear_solver->setOption(option);
        this->setParameter(Head, problem.getIC(0));
    }

    int solveTimeStep(const TimeStep &time)
    {
        _solHead->solveTimeStep(time);
        _vec_parameters[Head] = _solHead->getCurrentSolution(0);
        return 0;
    }

    double suggestNext(const TimeStep &time_current) { return _solHead->suggestNext(time_current); }

    bool isAwake(const TimeStep &time) { return _solHead->isAwake(time);  }

    void accept(const TimeStep &time)
    {
        _solHead->accept(time);
    };

private:
    GWFemProblem* _problem;
    SolutionForHead* _solHead;
    Rectangle *_rec;
    FemNodalFunctionScalar *_head;
    LagrangianFeObjectContainer* _feObjects;

    DISALLOW_COPY_AND_ASSIGN(FunctionHead);
};





