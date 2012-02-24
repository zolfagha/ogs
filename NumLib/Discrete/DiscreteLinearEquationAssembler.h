
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "NumLib/Discrete/ElementLocalAssembler.h"
#include "NumLib/Discrete/DoF.h"
#include "NumLib/TimeStepping/TimeStep.h"


namespace NumLib
{

/**
 * \brief Interface of discrete system assembler classes
 */
class IDiscreteLinearEquationAssembler
{
public:
    /// assembly
    virtual void assembly( MeshLib::IMesh &msh, DofMapManager &dofManager, MathLib::ILinearEquations &eqs) = 0;
};

/**
 * \brief Element-based discrete system assembler classes
 */
class ElementBasedAssembler : public IDiscreteLinearEquationAssembler
{
public:
    ///
    ElementBasedAssembler(IElemenetLocalAssembler &a) : _e_assembler(&a) {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(MeshLib::IMesh &msh, DofMapManager &dofManager, MathLib::ILinearEquations &eqs);

private:
    IElemenetLocalAssembler* _e_assembler;
};

/**
 * \brief Element-based discrete system assembler classes
 */
class ElementBasedTransientAssembler : public IDiscreteLinearEquationAssembler
{
public:
    ///
    ElementBasedTransientAssembler(const TimeStep &time, std::vector<std::vector<double>*> &u0, ITransientElemenetLocalAssembler &a) 
        : _transient_e_assembler(&a), _timestep(&time), _u0(&u0)
    { };


    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(MeshLib::IMesh &msh, DofMapManager &dofManager, MathLib::ILinearEquations &eqs);

private:
    ITransientElemenetLocalAssembler* _transient_e_assembler;
    const TimeStep* _timestep;
    std::vector<std::vector<double>*>* _u0;
};

class NodeBasedAssembler  : public IDiscreteLinearEquationAssembler
{
public:
    void assembly( MeshLib::IMesh& _msh, MathLib::ILinearEquations& _eqs) {};
};

class EdgeBasedAssembler  : public IDiscreteLinearEquationAssembler
{
public:
    void assembly( MeshLib::IMesh& _msh, MathLib::ILinearEquations& _eqs) {};
};

}
