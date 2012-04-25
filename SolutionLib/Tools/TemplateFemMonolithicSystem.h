
#pragma once

#include "Base/CodingTools.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "SolutionLib/Problem/FemIVBVProblem.h"
#include "SolutionLib/Solution/SingleStepFEM.h"


namespace SolutionLib
{

template <
	class T_FEM_IVBV_PROBLEM,
	class T_LINEAR_SOLVER
	>
class TemplateFemMonolithicSystem : public NumLib::TemplateTransientMonolithicSystem
{
public:
    typedef T_FEM_IVBV_PROBLEM MyProblemType;
    typedef SingleStepFEM
    		<
    			MyProblemType,
    			T_LINEAR_SOLVER
    		> MySolutionType;

    TemplateFemMonolithicSystem(size_t n_in, size_t n_out) : _solution(0)
    {
        TemplateTransientMonolithicSystem::resizeInputParameter(n_in);
        TemplateTransientMonolithicSystem::resizeOutputParameter(n_out);
    };

    virtual ~TemplateFemMonolithicSystem()
    {
        Base::releaseObject(_solution);
    }

    void define(DiscreteLib::DiscreteSystem* dis, MyProblemType* problem, Base::Options &option)
    {
        _solution = new MySolutionType(dis, problem);
        T_LINEAR_SOLVER* linear_solver = _solution->getLinearEquationSolver();
        linear_solver->setOption(option);
        this->setOutput(0, problem->getIC(0));
    }

    int solveTimeStep(const NumLib::TimeStep &time)
    {
        _solution->solveTimeStep(time);
        setOutput(0, _solution->getCurrentSolution(0));
        return 0;
    }

    double suggestNext(const NumLib::TimeStep &time_current)
    {
    	return _solution->suggestNext(time_current);
    }

    bool isAwake(const NumLib::TimeStep &time) { return _solution->isAwake(time);  }

    void accept(const NumLib::TimeStep &time)
    {
        _solution->accept(time);

        //std::cout << "Head=" << std::endl;
        //_solHead->getCurrentSolution(0)->printout();
    };

private:
    DISALLOW_COPY_AND_ASSIGN(TemplateFemMonolithicSystem);

private:
    //MyProblemType* _problem;
    MySolutionType* _solution;
    //FemLib::LagrangianFeObjectContainer* _feObjects;

};
} //end
