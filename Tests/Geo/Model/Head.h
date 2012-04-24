
#pragma once

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "GeoLib/Shape/Rectangle.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "SolutionLib/Problem/FemIVBVProblem.h"
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

typedef FemIVBVProblem
		<
			ElementWiseTimeEulerEQSLocalAssembler<Geo::GroundwaterFlowTimeODELocalAssembler>,
			ElementWiseTimeEulerResidualLocalAssembler<Geo::GroundwaterFlowTimeODELocalAssembler>,
			Geo::GroundwaterFlowJacobianLocalAssembler
		> GWFemProblem;



template <
	class T_LINEAR_SOLVER
	>
class FunctionHead : public TemplateTransientMonolithicSystem
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
        TemplateTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteSystem &dis, GWFemProblem &problem, Base::Options &option)
    {
        //MeshLib::IMesh *msh = dis.getMesh();
        //size_t nnodes = msh->getNumberOfNodes();
        //solution algorithm
        _solHead = new SolutionForHead(dis, problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        T_LINEAR_SOLVER* linear_solver = _solHead->getLinearEquationSolver();
        linear_solver->setOption(option);
        this->setOutput(Head, problem.getIC(0));
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
    Rectangle *_rec;
    //FemNodalFunctionScalar *_head;
    LagrangianFeObjectContainer* _feObjects;

    DISALLOW_COPY_AND_ASSIGN(FunctionHead);
};


} //end


