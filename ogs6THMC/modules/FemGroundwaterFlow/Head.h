
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
#include "OutputIO/IOutput.h"

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
: public ProcessLib::TemplateTransientProcess<0,1>
{
    enum Out { Head=0 };
public:
    typedef SolutionLib::SingleStepFEM
            <
                GWFemProblem,
                MathLib::CRSLisSolver
            > MySolutionType;

    ///
    FunctionHead() 
    : _problem(0), _solution(0), _feObjects(0)
    {
    };

    ///
    virtual ~FunctionHead()
    {
        BaseLib::releaseObject(_problem, _solution, _feObjects);
        BaseLib::releaseObjectsInStdVector(_list_output);
    }

    /// initialize this process
    virtual void initialize(const BaseLib::Options &option);

    /// finalize but nothing to do here
    virtual void finalize() {};


protected:
    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual MySolutionType* getSolution() {return _solution;};

    virtual void output(const NumLib::TimeStep &time);


private:
    DISALLOW_COPY_AND_ASSIGN(FunctionHead);

private:
    GWFemProblem* _problem;
    MySolutionType* _solution;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    std::vector<IOutput*> _list_output;

//public:
//    static ProcessLib::ProcessInfo* const _pcs_info;
};




