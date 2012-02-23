
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
    ElementBasedAssembler(IElemenetLocalAssembler &a) : _e_assembler(&a) {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(MeshLib::IMesh &msh, DofMapManager &dofManager, MathLib::ILinearEquations &eqs)
    {
        MathLib::DenseLinearEquations localEQS;
        std::vector<size_t> ele_node_ids, ele_node_size_order, local_dofmap;
        const size_t n_ele = msh.getNumberOfElements();

        std::vector<double> local_u_n;
        for (size_t i=0; i<n_ele; i++) {
            MeshLib::IElement *e = msh.getElemenet(i);
            // get dof map
            e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);
            e->getListOfNumberOfNodesForAllOrders(ele_node_size_order);
            dofManager.getListOfEqsID(ele_node_ids, ele_node_size_order, local_dofmap);
            // local assembly
            localEQS.create(local_dofmap.size());
            _e_assembler->assembly(*e, localEQS);
            // update global
            eqs.addA(local_dofmap, *localEQS.getA());
            eqs.addRHS(local_dofmap, localEQS.getRHS());
        }

        //apply ST
    }

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
    void assembly(MeshLib::IMesh &msh, DofMapManager &dofManager, MathLib::ILinearEquations &eqs)
    {
        const TimeStep &time = *_timestep;
        MathLib::DenseLinearEquations localEQS;
        std::vector<size_t> ele_node_ids, ele_node_size_order, local_dofmap;
        const size_t n_ele = msh.getNumberOfElements();

        std::vector<double> local_u_n;
        for (size_t i=0; i<n_ele; i++) {
            MeshLib::IElement *e = msh.getElemenet(i);
            // get dof map
            e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);
            e->getListOfNumberOfNodesForAllOrders(ele_node_size_order);
            dofManager.getListOfEqsID(ele_node_ids, ele_node_size_order, local_dofmap);
            // get previous time step results
            dofManager.getLocalVector(ele_node_ids, ele_node_size_order, *_u0, local_u_n);
            // local assembly
            localEQS.create(local_dofmap.size());
            _transient_e_assembler->assembly(time, *e, local_u_n, localEQS);
            // update global
            eqs.addA(local_dofmap, *localEQS.getA());
            eqs.addRHS(local_dofmap, localEQS.getRHS());
        }

        //apply ST
    }

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
