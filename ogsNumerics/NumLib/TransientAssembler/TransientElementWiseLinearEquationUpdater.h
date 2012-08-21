/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TransientElementWiseLinearEquationUpdater.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "MeshLib/Core/IElement.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Utils/Tools.h"

namespace NumLib
{

/**
 * \brief Element-wise updater for transient case
 *
 * \tparam T_LOCAL  Local assembler
 */
template <class T_LOCAL>
class TransientElementWiseLinearEquationUpdater
{
public:
    typedef T_LOCAL LocalAssemblerType;
    typedef DiscreteLib::IDiscreteVector<double> GlobalVector;

    /**
     *
     * @param time          Time step
     * @param msh           Mesh
     * @param dofManager    Dof manager
     * @param u0            Previous time step data
     * @param u1            Current time step data (guess)
     * @param a             Local assembler
     */
    TransientElementWiseLinearEquationUpdater   (
            const TimeStep* time, MeshLib::IMesh* msh,
            DiscreteLib::DofEquationIdTable* dofManager,
            const GlobalVector* u0, const GlobalVector* u1,
            LocalAssemblerType* a   )
    : _msh(msh), _dofManager(dofManager),
      _transient_e_assembler(new LocalAssemblerType(*a)), _timestep(time),
      _vec_u0(u0), _vec_u1(u1)
    {
    }

    /**
     * Copy constructor
     *
     * @param obj
     */
    TransientElementWiseLinearEquationUpdater(const TransientElementWiseLinearEquationUpdater<LocalAssemblerType> &obj)
    {
        _msh = obj._msh;
        _dofManager = obj._dofManager;
        _transient_e_assembler = new LocalAssemblerType(*obj._transient_e_assembler);
        _timestep = obj._timestep;
        _vec_u0 = obj._vec_u0;
        _vec_u1 = obj._vec_u1;
    }

    /**
     *
     * @param obj
     * @return
     */
    TransientElementWiseLinearEquationUpdater<LocalAssemblerType> &operator=(const TransientElementWiseLinearEquationUpdater<LocalAssemblerType> &obj)
    {
        return * new TransientElementWiseLinearEquationUpdater<LocalAssemblerType>(obj);
    }

    /**
     *
     */
    virtual ~TransientElementWiseLinearEquationUpdater()
    {
        BaseLib::releaseObject(_transient_e_assembler);
    }

    /**
     * update a global equation for the given element
     *
     * @param e     mesh element
     * @param eqs   global equation
     */
    template <class T_LINEAR_EQS>
    void update(const MeshLib::IElement &e, T_LINEAR_EQS &eqs)
    {
        std::vector<size_t> ele_node_ids;
        std::vector<size_t> local_dofmap_row, local_dofmap_column;
        LocalEquation localEQS;
        LocalVector local_u_n1;
        LocalVector local_u_n;

        // get dof map
        e.getNodeIDList(e.getMaximumOrder(), ele_node_ids);
        _dofManager->mapEqsID(_msh->getID(), ele_node_ids, local_dofmap_row, local_dofmap_column);
        // previous and current results
        DiscreteLib::getLocalVector(local_dofmap_row, *_vec_u1, local_u_n1);
        DiscreteLib::getLocalVector(local_dofmap_row, *_vec_u0, local_u_n);

        // local assembly
        localEQS.create(local_dofmap_row.size());
        _transient_e_assembler->assembly(*_timestep, e, local_u_n1, local_u_n, localEQS);

//        if (i<3) {
//            std::cout << "local A = \n" << *localEQS.getA();
//            std::cout << "local b = \n" << *localEQS.getRHS();
//        }

        // update global
        eqs.addAsub(local_dofmap_row, local_dofmap_column, *localEQS.getA());
        eqs.addRHSsub(local_dofmap_row, localEQS.getRHS());
    }

#if 0 // in progress for coupled problem. could be removed later
    template <class T_LINEAR_EQS>
    void update2(const MeshLib::IElement &e, T_LINEAR_EQS &eqs)
    {
        LocalEquation localEQS;
        size_t n_var;

        std::vector<std::vector<size_t> > local_dofmap_row(n_var), local_dofmap_column(n_var);
        std::vector<size_t> vec_var_order(n_var);
        std::vector<size_t> vec_var_id(n_var);

        std::vector<LocalVector> vec_local_u0(n_var), vec_local_u1(n_var);
        for (size_t i=0; i<n_var; i++) {
            std::vector<size_t> ele_node_ids;
            e.getNodeIDList(vec_var_order[i], ele_node_ids);
            _dofManager->mapEqsID(vec_var_id[i], _msh->getID(), ele_node_ids, local_dofmap_row[i]);
            DiscreteLib::getLocalVector(local_dofmap_row[i], *_vec_u1, vec_local_u1[i]);
            DiscreteLib::getLocalVector(local_dofmap_row[i], *_vec_u0, vec_local_u0[i]);
        }

        // local assembly
        localEQS.create(local_dofmap_row.size());
        _transient_e_assembler->assembly(*_timestep, e, vec_local_u1, vec_local_u0, localEQS);

        // update global
        eqs.addAsub(local_dofmap_row, local_dofmap_column, *localEQS.getA());
        eqs.addRHSsub(local_dofmap_row, localEQS.getRHS());
    }
#endif


private:
    MeshLib::IMesh* _msh;
    DiscreteLib::DofEquationIdTable* _dofManager;
    LocalAssemblerType* _transient_e_assembler;
    const TimeStep* _timestep;
    const GlobalVector* _vec_u0;
    const GlobalVector* _vec_u1;
};

} //end
