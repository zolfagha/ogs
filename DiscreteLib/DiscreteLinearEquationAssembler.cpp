
#include "DiscreteLinearEquationAssembler.h"

#include <vector>

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

#include "MeshLib/Core/IMesh.h"



namespace DiscreteLib
{

void ElementBasedAssembler::assembly(MeshLib::IMesh &msh, DofMapManager &dofManager, MathLib::ILinearEquations &eqs)
{
    MathLib::DenseLinearEquations localEQS;
    std::vector<size_t> ele_node_ids, ele_node_size_order;
    std::vector<long> local_dofmap;
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
        eqs.addAsub(local_dofmap, *localEQS.getA());
        eqs.addRHSsub(local_dofmap, localEQS.getRHS());
    }

    //apply ST
};

} //end

