/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTransientLinearEQSAssembler.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "ElementWiseTransientLinearEQSAssembler.h"

#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Utils/Tools.h"
#include "NumLib/TimeStepping/TimeStep.h"

namespace NumLib
{

void ElementWiseTransientLinearEQSAssembler::assembly(MeshLib::IMesh &msh, DiscreteLib::DofEquationIdTable &dofManager, MathLib::ILinearEquation &eqs)
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
        _transient_e_assembler->assembly(time, *e, local_u_n1, local_u_n, localEQS);

//        if (i<3) {
//            std::cout << "local A = \n" << *localEQS.getA();
//            std::cout << "local b = \n" << *localEQS.getRHS();
//        }

        // update global
        eqs.addAsub(local_dofmap, *localEQS.getA());
        eqs.addRHSsub(local_dofmap, localEQS.getRHS());
    }

    //apply ST
}


} //end

