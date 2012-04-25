
#include "ElementWiseTransientResidualAssembler.h"

#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Utils/Tools.h"

namespace NumLib
{

void ElementWiseTransientResidualAssembler::assembly( const MeshLib::IMesh &msh, const DiscreteLib::DofEquationIdTable &dofManager, GlobalVectorType &globalVec)
{
    const TimeStep &time = *_timestep;
	LocalVectorType localVec;
    std::vector<size_t> ele_node_ids, ele_node_size_order;
    std::vector<size_t> local_dofmap;
    const size_t n_ele = msh.getNumberOfElements();

    LocalVectorType local_u_n1;
    LocalVectorType local_u_n;
    for (size_t i=0; i<n_ele; i++) {
        MeshLib::IElement *e = msh.getElemenet(i);
        // get dof map
        e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);
        e->getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        dofManager.mapEqsID(msh.getID(), ele_node_ids, local_dofmap);
        // previous and current results
        DiscreteLib::getLocalVector(dofManager, ele_node_ids, ele_node_size_order, *_u1, local_u_n1);
        DiscreteLib::getLocalVector(dofManager, ele_node_ids, ele_node_size_order, *_u0, local_u_n);
        // local assembly
        localVec.resize(local_dofmap.size(), .0);
        _transient_e_assembler->assembly(time, *e, local_u_n1, local_u_n, localVec);
        // update global
        globalVec.addSubvector(local_dofmap, &localVec[0]);
    }

}

} //end
