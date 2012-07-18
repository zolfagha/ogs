
#pragma once

#include <vector>
#include <valarray>

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/EquationId/DofEquationIdTable.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "DiscreteLib/Assembler/IDiscreteVectorAssembler.h"
#include "DiscreteLib/Utils/Tools.h"
#include "NumLib/TimeStepping/TimeStep.h"

#include "IElementWiseTransientJacobianLocalAssembler.h"

namespace NumLib
{

//class TimeStep;


/**
 * \brief Element-based discrete system assembler classes
 */
class ElementWiseTransientDxEQSAssembler : public DiscreteLib::IDiscreteLinearEquationAssembler
{
public:
    /// @param u0
    /// @param u1
    /// @param a
    ElementWiseTransientDxEQSAssembler(const TimeStep* time, const std::vector<GlobalVector*>* u0, const std::vector<GlobalVector*>* u1, IElementWiseTransientJacobianLocalAssembler* a)
        : _transient_e_assembler(a), _timestep(time), _u0(u0), _u1(u1)
    { };


    /// Conduct the element by element assembly procedure
    ///
    /// @param msh                 Mesh
    /// @param dofManager         Dof map manager
    /// @param list_dofId         List of Dof IDs used in this problem
    /// @param J                 Jacobian matrix
    virtual void assembly( MeshLib::IMesh &msh, DiscreteLib::DofEquationIdTable &dofManager, MathLib::ILinearEquations &eqs)
    {
        const TimeStep &time = *_timestep;
        LocalEquation localEQS;
        std::vector<size_t> ele_node_ids, ele_node_size_order;
        std::vector<size_t> local_dofmap;
        const size_t n_ele = msh.getNumberOfElements();

        LocalVector local_u_n1;
        LocalVector local_u_n;
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
            localEQS.create(local_dofmap.size());
            _transient_e_assembler->assembly(time, *e, local_u_n1, local_u_n, *localEQS.getA());
            // update global
            eqs.addAsub(local_dofmap, *localEQS.getA());
        }
    }

private:
    IElementWiseTransientJacobianLocalAssembler* _transient_e_assembler;
    const TimeStep* _timestep;
    const std::vector<GlobalVector*>* _u0;
    const std::vector<GlobalVector*>* _u1;
};


}
