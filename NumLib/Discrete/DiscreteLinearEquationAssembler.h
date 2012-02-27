
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/DiscreteLinearEquationAssembler.h"
#include "DiscreteLib/DoF.h"
#include "DiscreteLib/DiscreteVector.h"

#include "ElementLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"


namespace NumLib
{

/**
 * \brief Element-based discrete system assembler classes
 */
class ElementBasedTransientAssembler : public DiscreteLib::IDiscreteLinearEquationAssembler
{
public:
    ///
    ElementBasedTransientAssembler(const TimeStep &time, std::vector<DiscreteLib::DiscreteVector<double>*> &u0, ITransientElemenetLocalAssembler &a) 
        : _transient_e_assembler(&a), _timestep(&time), _u0(&u0)
    { };


    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(MeshLib::IMesh &msh, DiscreteLib::DofMapManager &dofManager, MathLib::ILinearEquations &eqs);

private:
    ITransientElemenetLocalAssembler* _transient_e_assembler;
    const TimeStep* _timestep;
    std::vector<DiscreteLib::DiscreteVector<double>*>* _u0;
};


}
