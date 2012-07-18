
#pragma once

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "SolutionLib/FemProblem/FemIVBVProblem.h"
#include "SolutionLib/FemProblem/AbstractTransientFemFunction.h"


#include "GWTimeODELocalAssembler.h"
#include "GWJacobianLocalAssembler.h"
#include "ProcessLib/TemplateTransientProcess.h"

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
: public ProcessLib::TemplateTransientProcess<0,1> // SolutionLib::AbstractTransientFemFunction<0, 1>
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

    virtual ~FunctionHead()
    {
        BaseLib::releaseObject(_problem, _solution, _feObjects);
    }

    virtual void initialize(const BaseLib::Options &option);

    virtual void finalize() {};


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


//public:
//    static ProcessLib::ProcessInfo* const _pcs_info;
};




