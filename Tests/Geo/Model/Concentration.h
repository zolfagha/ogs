
#pragma once

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "GeoLib/Shape/Rectangle.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "SolutionLib/Problem/FemIVBVProblem.h"
#include "SolutionLib/Solution/SingleStepFEM.h"

#include "Tests/Geo/Equation/FemMassTransport.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Material/Compound.h"

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
			ElementWiseTimeEulerEQSLocalAssembler<Geo::MassTransportTimeODELocalAssembler>,
			ElementWiseTimeEulerResidualLocalAssembler<Geo::MassTransportTimeODELocalAssembler>,
			Geo::MassTransportJacobianLocalAssembler
		> MassFemProblem;

template <
	class T_LINEAR_SOLVER
	>
class FunctionConcentration : public TemplateTransientMonolithicSystem
{
    enum In { Velocity=0 };
    enum Out { Concentration = 0 };

public:
    typedef SingleStepFEM
    		<
    			MassFemProblem,
    			T_LINEAR_SOLVER
    		> SolutionForConc;

    FunctionConcentration() 
    {
        TemplateTransientMonolithicSystem::resizeInputParameter(1);
        TemplateTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteSystem* dis, MassFemProblem* problem, Base::Options &option)
    {
        //solution algorithm
        _solConc = new SolutionForConc(dis, problem);
        //_solConc->getTimeODEAssembler()->setTheta(1.0);
        SolutionForConc::LinearSolverType* linear_solver = _solConc->getLinearEquationSolver();
        linear_solver->setOption(option);
        this->setOutput(Concentration, problem->getIC(0));
    }

    int solveTimeStep(const TimeStep &time)
    {
        //input
        const MathLib::SpatialFunctionVector *vel = this->getInput<MathLib::SpatialFunctionVector>(Velocity);
        //TODO _solConc->getProblem()->getResidualAssembler().initialize(vel);

        // solve
        _solConc->solveTimeStep(time);

        // output
        setOutput(Concentration, _solConc->getCurrentSolution(0));

        return 0;
    }

    double suggestNext(const TimeStep &time_current) { return _solConc->suggestNext(time_current); }

    bool isAwake(const TimeStep &time) { return _solConc->isAwake(time);  }

    void accept(const TimeStep &time)
    {
    	_solConc->accept(time);

        //std::cout << "Concentration=" << std::endl;
        //_solConc->getCurrentSolution(0)->printout();
    };

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionConcentration);

    SolutionForConc* _solConc;
    //FemNodalFunctionScalar *_conc;
    LagrangianFeObjectContainer* _feObjects;
};

} //end
