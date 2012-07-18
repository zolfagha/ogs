
#pragma once

#include "FemLinearEQS.h"
#include "FemResidualEQS.h"
#include "FemDxEQS.h"

namespace SolutionLib
{

template
    <
    class T_LOCAL_ASSEMBLER_LINEAR,
    class T_LOCAL_ASSEMBLER_RESIDUAL,
    class T_LOCAL_ASSEMBLER_JACOBIAN
    >
class TemplateFemEquation
{
public:
    typedef T_LOCAL_ASSEMBLER_LINEAR LinearAssemblerType;
    typedef T_LOCAL_ASSEMBLER_RESIDUAL ResidualAssemblerType;
    typedef T_LOCAL_ASSEMBLER_JACOBIAN JacobianAssemblerType;
    typedef TemplateTransientLinearFEMFunction<LinearAssemblerType> LinearEQSType;
    typedef TemplateTransientResidualFEMFunction<ResidualAssemblerType> ResidualEQSType;
    typedef TemplateTransientDxFEMFunction<JacobianAssemblerType> DxEQSType;

    TemplateFemEquation(
            LinearAssemblerType *linear_assembly,
            ResidualAssemblerType *residual_assembly,
            JacobianAssemblerType *jacobian_assembly
    ) :
    _linear_assembler(linear_assembly),
    _residual_assembler(residual_assembly),
    _jacobian_assembler(jacobian_assembly)
    {};

    ///
    LinearAssemblerType* getLinearAssembler() const { return _linear_assembler; }

    ///
    ResidualAssemblerType* getResidualAssembler() const { return _residual_assembler; }

    ///
    JacobianAssemblerType* getJacobianAssembler() const { return _jacobian_assembler; }

private:
    LinearAssemblerType* _linear_assembler;
    ResidualAssemblerType* _residual_assembler;
    JacobianAssemblerType* _jacobian_assembler;
};




} //end
