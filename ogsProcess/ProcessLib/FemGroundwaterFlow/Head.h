
#pragma once

#include "BaseLib/CodingTools.h"
//#include "MathLib/LinAlg/VectorNorms.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
//#include "FemLib/Function/FemFunction.h"
//#include "NumLib/Function/Function.h"
//#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "SolutionLib/FemProblem/FemIVBVProblem.h"
#include "SolutionLib/FemProblem/AbstractTransientFemFunction.h"
#include "ProcessBuilder.h"


#include "GWTimeODELocalAssembler.h"
#include "GWJacobianLocalAssembler.h"


namespace Geo
{

//--------------------------------------------------------------------------------------------------
// Equation definition
//--------------------------------------------------------------------------------------------------
typedef SolutionLib::TemplateFemEquation <
		Geo::GroundwaterFlowTimeODELocalAssembler<
			NumLib::ElementWiseTimeEulerEQSLocalAssembler>,
		Geo::GroundwaterFlowTimeODELocalAssembler<
			NumLib::ElementWiseTimeEulerResidualLocalAssembler>,
		Geo::GroundwaterFlowJacobianLocalAssembler
		>
		GWFemEquation;


//--------------------------------------------------------------------------------------------------
// IVBV problem definition
//--------------------------------------------------------------------------------------------------
typedef SolutionLib::FemIVBVProblem< GWFemEquation > GWFemProblem;



//--------------------------------------------------------------------------------------------------
// Function definition
//--------------------------------------------------------------------------------------------------
class FunctionHead
: public SolutionLib::AbstractTransientFemFunction<0, 1>
{
    enum Out { Head=0 };
public:
    typedef SolutionLib::SingleStepFEM
    		<
    			GWFemProblem,
    			MathLib::CRSLisSolver
    		> MySolutionType;

    FunctionHead() 
    : _problem(0), _solution(0), _feObjects(0)
    {
    };

    void define(DiscreteLib::DiscreteSystem* dis, BaseLib::Options &option)
    {
    	MaterialLib::PorousMedia *pm = 0;

        _feObjects = new FemLib::LagrangianFeObjectContainer(*dis->getMesh());
        Geo::GWFemEquation::LinearAssemblerType* linear_assembler = new Geo::GWFemEquation::LinearAssemblerType(*_feObjects, *pm);
        Geo::GWFemEquation::ResidualAssemblerType* r_assembler = new Geo::GWFemEquation::ResidualAssemblerType(*_feObjects, *pm);
        Geo::GWFemEquation::JacobianAssemblerType* j_eqs = new Geo::GWFemEquation::JacobianAssemblerType(*_feObjects, *pm);
        Geo::GWFemEquation* eqs = new  Geo::GWFemEquation(linear_assembler, r_assembler, j_eqs);
        Geo::GWFemProblem* _problem = new Geo::GWFemProblem(dis);
        _problem->setEquation(eqs);
    	_solution = new MySolutionType(dis, _problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
    	MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
        linear_solver->setOption(option);
        _solution->getNonlinearSolver()->setOption(option);

        this->setOutput(Head, _problem->getVariable(0)->getIC());
    }


protected:
    virtual void updateOutput()
    {
        setOutput(Head, _solution->getCurrentSolution(0));
    }
    virtual MySolutionType* getSolution() {return _solution;};

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionHead);

private:
    GWFemProblem* _problem;
    MySolutionType* _solution;
    FemLib::LagrangianFeObjectContainer* _feObjects;


public:
    static ProcessLib::ProcessInfo* const _pcs_info;
};


} //end

OGS_PROCESS(GROUNDWATER_FLOW, Geo::FunctionHead);

