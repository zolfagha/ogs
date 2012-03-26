
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"


#include "DiscreteLib/EquationId/DofEquationIdTable.h"
#include "DiscreteLib/Core/DiscreteVector.h"
#include "DiscreteLib/Assembler/DiscreteLinearEquationAssembler.h"

#include "ElementLocalAssembler.h"

namespace MeshLib
{
class IMesh;
}

namespace NumLib
{

class TimeStep;


/**
 * \brief Element-based discrete system assembler classes
 */
class ElementBasedTransientAssembler : public DiscreteLib::IDiscreteLinearEquationAssembler
{
public:
    /// @param time
    /// @param u0
    /// @param u1
    /// @param a
    ElementBasedTransientAssembler(const TimeStep &time, const std::vector<DiscreteLib::DiscreteVector<double>*> &u0, const std::vector<DiscreteLib::DiscreteVector<double>*> &u1, ITransientElemenetLocalAssembler &a)
        : _transient_e_assembler(&a), _timestep(&time), _u0(&u0), _u1(&u1)
    { };


    /// Conduct the element by element assembly procedure
    ///
    /// @param msh 				Mesh
    /// @param dofManager 		Dof map manager
    /// @param list_dofId 		List of Dof IDs used in this problem
    /// @param eqs 				Linear equation solver
    void assembly(MeshLib::IMesh &msh, DiscreteLib::DofEquationIdTable &dofManager, MathLib::ILinearEquations &eqs);

private:
    ITransientElemenetLocalAssembler* _transient_e_assembler;
    const TimeStep* _timestep;
    const std::vector<DiscreteLib::DiscreteVector<double>*>* _u0;
    const std::vector<DiscreteLib::DiscreteVector<double>*>* _u1;
};


}
