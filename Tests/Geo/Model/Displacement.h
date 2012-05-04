
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

#include "Tests/Geo/Equation/FemLinearElasticity.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Material/Solid.h"

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
			Geo::FemLinearElasticLinearLocalAssembler,
			Geo::FemLinearElasticResidualLocalAssembler,
			Geo::FemLinearElasticJacobianLocalAssembler
		> FemLinearElasticProblem;



template <
	class T_LINEAR_SOLVER
	>
class FunctionDisplacement : public TemplateTransientMonolithicSystem
{
    enum Out { Displacement=0 };
public:
    typedef SingleStepFEM
    		<
    			FemLinearElasticProblem,
    			T_LINEAR_SOLVER
    		> MySolution;

    FunctionDisplacement()
    {
        TemplateTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteSystem* dis, FemLinearElasticProblem* problem, Base::Options &option)
    {
        //solution algorithm
        _sol = new MySolution(dis, problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        typename MySolution::LinearSolverType* linear_solver = _sol->getLinearEquationSolver();
        linear_solver->setOption(option);
        _sol->getNonlinearSolver()->setOption(option);
        this->setOutput(Displacement, problem->getIC(0));
    }

    int solveTimeStep(const TimeStep &time)
    {
        _sol->solveTimeStep(time);
        setOutput(Displacement, _sol->getCurrentSolution(0));
        return 0;
    }

    double suggestNext(const TimeStep &time_current) { return _sol->suggestNext(time_current); }

    bool isAwake(const TimeStep &time) { return _sol->isAwake(time);  }

    void accept(const TimeStep &time)
    {
        _sol->accept(time);

        //std::cout << "Head=" << std::endl;
        //_solHead->getCurrentSolution(0)->printout();
    };

private:
    FemLinearElasticProblem* _problem;
    MySolution* _sol;
    LagrangianFeObjectContainer* _feObjects;

    DISALLOW_COPY_AND_ASSIGN(FunctionDisplacement);
};


} //end


