
#include "ElementWiseTransientLinearEQSAssembler.h"

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Utils/Tools.h"
#include "NumLib/TimeStepping/TimeStep.h"

namespace NumLib
{

void ElementWiseTransientLinearEQSAssembler::assembly(MeshLib::IMesh &msh, DiscreteLib::DofEquationIdTable &dofManager, MathLib::ILinearEquations &eqs)
{
    const TimeStep &time = *_timestep;
    MathLib::DenseLinearEquations localEQS;
    std::vector<size_t> ele_node_ids, ele_node_size_order;
    std::vector<long> local_dofmap;
    const size_t n_ele = msh.getNumberOfElements();

    LocalVectorType local_u_n1;
    LocalVectorType local_u_n;
    for (size_t i=0; i<n_ele; i++) {
        MeshLib::IElement *e = msh.getElemenet(i);
        // get dof map
        e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);
        e->getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        dofManager.mapEqsID(msh.getID(), ele_node_ids, local_dofmap);
        //dofManager.mapEqsID(ele_node_ids, ele_node_size_order, local_dofmap);
        // previous and current results
        DiscreteLib::getLocalVector(dofManager, ele_node_ids, ele_node_size_order, *_u1, local_u_n1);
        DiscreteLib::getLocalVector(dofManager, ele_node_ids, ele_node_size_order, *_u0, local_u_n);
        // local assembly
        localEQS.create(local_dofmap.size());
        _transient_e_assembler->assembly(time, *e, local_u_n1, local_u_n, localEQS);
        // update global
        eqs.addAsub(local_dofmap, *localEQS.getA());
        eqs.addRHSsub(local_dofmap, localEQS.getRHS());
    }

    //apply ST
}


} //end

