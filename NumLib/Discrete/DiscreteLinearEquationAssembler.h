
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "NumLib/Discrete/ElementLocalAssembler.h"

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
    ElementBasedAssembler(const TimeStep &time, ITransientElemenetLocalAssembler &a) 
        : _element_assembler(&a), _timestep(&time)
    {

    }

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(MeshLib::IMesh &msh, DofMapManager &dofManager, MathLib::ILinearEquations &eqs)
    {
        const TimeStep &time = *_timestep;
        MathLib::DenseLinearEquations localEQS;
        std::vector<size_t> ele_node_ids, ele_node_size_order, local_dofmap;
        const size_t n_ele = msh.getNumberOfElements();

        for (size_t i=0; i<n_ele; i++) {
            MeshLib::IElement *e = msh.getElemenet(i);
            // get dof map
            e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);
            e->getListOfNumberOfNodesForAllOrders(ele_node_size_order);
            dofManager.getListOfEqsID(ele_node_ids, ele_node_size_order, local_dofmap);
            // local assembly
            localEQS.create(local_dofmap.size());
            _element_assembler->assembly(time, *e, localEQS);
            // update global
            eqs.addA(local_dofmap, *localEQS.getA());
            eqs.addRHS(local_dofmap, localEQS.getRHS());
        }

        //apply ST
    }

private:
    ITransientElemenetLocalAssembler* _element_assembler;
    const TimeStep* _timestep;
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
