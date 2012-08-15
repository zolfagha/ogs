/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTransientDxEQSAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <valarray>

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "DiscreteLib/Core/IDiscreteVectorAssembler.h"
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
    ElementWiseTransientDxEQSAssembler(const TimeStep* time, const GlobalVector* u0, const GlobalVector* u1, IElementWiseTransientJacobianLocalAssembler* a)
        : _transient_e_assembler(a), _timestep(time), _vec_u0(u0), _vec_u1(u1)
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
            dofManager.mapEqsID(msh.getID(), ele_node_ids, local_dofmap, DiscreteLib::DofNumberingType::BY_POINT);
            // previous and current results
            DiscreteLib::getLocalVector(local_dofmap, *_vec_u1, local_u_n1);
            DiscreteLib::getLocalVector(local_dofmap, *_vec_u0, local_u_n);
            // local assembly
            localEQS.create(local_dofmap.size());
            _transient_e_assembler->assembly(time, *e, local_u_n1, local_u_n, *localEQS.getA());
            // update global
            eqs.addAsub(local_dofmap, *localEQS.getA());

//            if (i<2)
//                std::cout << "local A = \n" << *localEQS.getA() << std::endl;
        }
    }

private:
    IElementWiseTransientJacobianLocalAssembler* _transient_e_assembler;
    const TimeStep* _timestep;
    const GlobalVector* _vec_u0;
    const GlobalVector* _vec_u1;
};


}
