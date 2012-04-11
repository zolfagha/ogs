
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


typedef FemIVBVProblem<Geo::WeakFormMassTransport> MassFemProblem;

template <
	template <class> class T_NONLINEAR,
	class T_LINEAR_SOLVER
	>
class FunctionConcentration : public TemplateTransientMonolithicSystem<2,1>
{
public:
    typedef TemplateTransientLinearFEMFunction<
    			FemIVBVProblem,
    			TimeEulerElementLocalAssembler,
    			T_LINEAR_SOLVER,
    			Geo::WeakFormMassTransport
			> MyLinearFunction;

    typedef T_NONLINEAR<MyLinearFunction> MyNonlinearFunction;

    typedef SingleStepFEM
    		<
    			FemIVBVProblem<Geo::WeakFormMassTransport>,
    			//TimeEulerElementLocalAssembler,
    			MyLinearFunction,
    			MyNonlinearFunction,
    			T_LINEAR_SOLVER
    		> SolutionForConc;

    enum Parameters { Velocity=0, Concentration = 1 };

    FunctionConcentration() {};

    void define(DiscreteSystem &dis, MassFemProblem &problem, Base::Options &option)
    {
        //MeshLib::IMesh *msh = dis.getMesh();
        //size_t nnodes = msh->getNumberOfNodes();
        //solution algorithm
        _solConc = new SolutionForConc(dis, problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        T_LINEAR_SOLVER* linear_solver = _solConc->getLinearEquationSolver();
        linear_solver->setOption(option);
        this->setOutput(Concentration, problem.getIC(0));
    }

    int solveTimeStep(const TimeStep &time)
    {
        const MathLib::SpatialFunctionVector *vel = this->getInput<MathLib::SpatialFunctionVector>(Velocity);
        _solConc->getProblem()->getElementAssemlby().initialize(vel);
    	_solConc->solveTimeStep(time);
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
