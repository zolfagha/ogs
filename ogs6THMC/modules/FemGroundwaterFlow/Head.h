
#pragma once

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "SolutionLib/FemProblem/FemIVBVProblem.h"
#include "SolutionLib/FemProblem/AbstractTransientFemFunction.h"
#include "ProcessLib/TemplateTransientProcess.h"

#include "GWTimeODELocalAssembler.h"
#include "GWJacobianLocalAssembler.h"

//--------------------------------------------------------------------------------------------------
// Equation definition
//--------------------------------------------------------------------------------------------------
typedef SolutionLib::TemplateFemEquation <
        GroundwaterFlowTimeODELocalAssembler<
            NumLib::ElementWiseTimeEulerEQSLocalAssembler>,
        GroundwaterFlowTimeODELocalAssembler<
            NumLib::ElementWiseTimeEulerResidualLocalAssembler>,
        GroundwaterFlowJacobianLocalAssembler
        >
        FemGWEquation;


//--------------------------------------------------------------------------------------------------
// IVBV problem definition
//--------------------------------------------------------------------------------------------------
typedef SolutionLib::FemIVBVProblem< FemGWEquation > FemGWProblem;



//--------------------------------------------------------------------------------------------------
// Function definition
//--------------------------------------------------------------------------------------------------
class FunctionHead
: public ProcessLib::TemplateTransientProcess<0,1>
{
    enum Out { Head=0 };
public:
    typedef SolutionLib::SingleStepFEM
            <
                FemGWProblem,
                MathLib::CRSLisSolver
            > MySolutionType;

    ///
    FunctionHead() 
    : TemplateTransientProcess<0,1>("GROUNDWATER_FLOW"), _problem(0), _solution(0), _feObjects(0)
    {
    };

    ///
    virtual ~FunctionHead()
    {
        BaseLib::releaseObject(_problem, _solution, _feObjects);
    }

    /// initialize this process
    virtual bool initialize(const BaseLib::Options &option);

    /// finalize but nothing to do here
    virtual void finalize() {};


protected:
    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual MySolutionType* getSolution() {return _solution;};

    virtual void output(const NumLib::TimeStep &time);


private:
    DISALLOW_COPY_AND_ASSIGN(FunctionHead);

private:
    FemGWProblem* _problem;
    MySolutionType* _solution;
    FemLib::LagrangianFeObjectContainer* _feObjects;
};




