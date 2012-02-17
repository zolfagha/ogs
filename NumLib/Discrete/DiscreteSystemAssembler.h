
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "NumLib/Discrete/ElementLocalAssembler.h"

namespace NumLib
{

class IDiscreteSystemAssembler
{
public:
    virtual void assembly( MeshLib::IMesh &msh, DofMapManager &dofManager, std::vector<size_t> &list_dofId, MathLib::ILinearEquations &eqs) = 0;
};


class ElementBasedAssembler : public IDiscreteSystemAssembler
{
private:
    IElemenetLocalAssembler* element_assembler;

public:
    ElementBasedAssembler(IElemenetLocalAssembler &a) : element_assembler(&a)
    {

    }

    void assembly( MeshLib::IMesh &msh, DofMapManager &dofManager, std::vector<size_t> &list_dofId, MathLib::ILinearEquations &eqs)
    {
        MathLib::DenseLinearEquations localEQS;
        std::vector<size_t> ele_node_ids, ele_node_size_order, local_dofmap;
        const size_t n_ele = msh.getNumberOfElements();

        for (size_t i=0; i<n_ele; i++) {
            MeshLib::IElement *e = msh.getElemenet(i);
            // get dof map
            e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);
            e->getListOfNumberOfNodesForAllOrders(ele_node_size_order);
            dofManager.getListOfEqsID(list_dofId, ele_node_ids, ele_node_size_order, local_dofmap);
            // local assembly
            localEQS.create(local_dofmap.size());
            element_assembler->assembly(*e, localEQS);
            // update global
            eqs.addA(local_dofmap, *localEQS.getA());
            eqs.addRHS(local_dofmap, localEQS.getRHS());
        }

        //apply ST
    }
};

class NodeBasedAssembler  : public IDiscreteSystemAssembler
{
public:
    void assembly( MeshLib::IMesh& _msh, MathLib::ILinearEquations& _eqs) {};
};

class EdgeBasedAssembler  : public IDiscreteSystemAssembler
{
public:
    void assembly( MeshLib::IMesh& _msh, MathLib::ILinearEquations& _eqs) {};
};

}
