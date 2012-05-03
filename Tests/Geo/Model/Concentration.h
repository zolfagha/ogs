
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

#include "Tests/Geo/Equation/FemMassTransport.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Material/Compound.h"

namespace Geo
{

typedef SolutionLib::FemIVBVProblem
		<
			Geo::MassTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler>,
			Geo::MassTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler>,
			Geo::MassTransportJacobianLocalAssembler
		> MassFemProblem;

template <
	class T_LINEAR_SOLVER
	>
class FunctionConcentration : public NumLib::TemplateTransientMonolithicSystem
{
    enum In { Velocity=0 };
    enum Out { Concentration = 0 };

public:
    typedef SolutionLib::SingleStepFEM
    		<
    			MassFemProblem,
    			T_LINEAR_SOLVER
    		> SolutionForConc;

    FunctionConcentration() 
    {
        TemplateTransientMonolithicSystem::resizeInputParameter(1);
        TemplateTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteLib::DiscreteSystem* dis, MassFemProblem* problem, Base::Options &option)
    {
        //solution algorithm
        _solConc = new SolutionForConc(dis, problem);
        //_solConc->getTimeODEAssembler()->setTheta(1.0);
        typename SolutionForConc::LinearSolverType* linear_solver = _solConc->getLinearEquationSolver();
        linear_solver->setOption(option);
        _solConc->getNonlinearSolver()->setOption(option);
        this->setOutput(Concentration, problem->getIC(0));
    }

    int solveTimeStep(const NumLib::TimeStep &time)
    {
        //input
        const MathLib::SpatialFunctionVector *vel = this->getInput<MathLib::SpatialFunctionVector>(Velocity);
        _solConc->getProblem()->getLinearAssembler()->initialize(vel);
        _solConc->getProblem()->getResidualAssembler()->initialize(vel);
        _solConc->getProblem()->getJacobianAssembler()->initialize(vel);

        // solve
        _solConc->solveTimeStep(time);

        // output
        setOutput(Concentration, _solConc->getCurrentSolution(0));

        return 0;
    }

    double suggestNext(const NumLib::TimeStep &time_current) { return _solConc->suggestNext(time_current); }

    bool isAwake(const NumLib::TimeStep &time) { return _solConc->isAwake(time);  }

    void accept(const NumLib::TimeStep &time)
    {
    	_solConc->accept(time);

        //std::cout << "Concentration=" << std::endl;
        //_solConc->getCurrentSolution(0)->printout();
    };

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionConcentration);

    SolutionForConc* _solConc;
    //FemNodalFunctionScalar *_conc;
    FemLib::LagrangianFeObjectContainer* _feObjects;
};

} //end
