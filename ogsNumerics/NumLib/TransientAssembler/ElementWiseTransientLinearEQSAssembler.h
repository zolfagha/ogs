
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"


#include "DiscreteLib/EquationId/DofEquationIdTable.h"
#include "DiscreteLib/Vector/DiscreteVector.h"
#include "DiscreteLib/Assembler/IDiscreteLinearEquationAssembler.h"

#include "IElementWiseTransientLinearEQSLocalAssembler.h"

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
class ElementWiseTransientLinearEQSAssembler : public DiscreteLib::IDiscreteLinearEquationAssembler
{
public:

    /// @param time
    /// @param u0
    /// @param u1
    /// @param a
    //ElementWiseTransientLinearEQSAssembler(const TimeStep* time, const std::vector<GlobalVector*>* u0, const std::vector<GlobalVector*>* u1, IElementWiseTransientLinearEQSLocalAssembler* a)
    ElementWiseTransientLinearEQSAssembler(const TimeStep* time, const GlobalVector* u0, const GlobalVector* u1, IElementWiseTransientLinearEQSLocalAssembler* a)
        : _transient_e_assembler(a), _timestep(time), _vec_u0(u0), _vec_u1(u1)
    { };


    /// Conduct the element by element assembly procedure
    ///
    /// @param msh                 Mesh
    /// @param dofManager         Dof map manager
    /// @param list_dofId         List of Dof IDs used in this problem
    /// @param eqs                 Linear equation solver
    void assembly(MeshLib::IMesh &msh, DiscreteLib::DofEquationIdTable &dofManager, MathLib::ILinearEquations &eqs);

private:
    IElementWiseTransientLinearEQSLocalAssembler* _transient_e_assembler;
    const TimeStep* _timestep;
    const GlobalVector* _vec_u0;
    const GlobalVector* _vec_u1;
    //const std::vector<GlobalVector*>* _vec_u0;
    //const std::vector<GlobalVector*>* _vec_u1;
};


}
