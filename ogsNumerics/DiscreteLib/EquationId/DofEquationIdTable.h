/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DofEquationIdTable.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <map>
#include <cassert>
#include <algorithm>

#include "BaseLib/CodingTools.h"

#include "EquationIdStorage.h"
#include "SequentialEquationIdStorage.h"
#include "RandomEquationIdStorage.h"


namespace DiscreteLib
{

/**
 * \brief Type of numbering DoF's equation index
 */
struct DofNumberingType
{
    enum type
    {
        BY_DOF,
        BY_POINT
    };
};


/**
 * \brief DoF and equation index table
 *
 * This class contains a table which relates DoFs and index in a equation.
 * A DoF is implicitly defined as a unique combination of the followings 
 * - variable type
 * - mesh id
 * - discrete point id in the mesh (e.g. node id, element id)
 * 
 */
class DofEquationIdTable
{
public:

    ///
    DofEquationIdTable() : _total_dofs(0), _total_dofs_without_ghosts(0), _is_seq_address(false) {};

    ///
    virtual ~DofEquationIdTable()
    {
        for (size_t i=0; i<_map_var2dof.size(); i++) {
            BaseLib::releaseObjectsInStdMap(_map_var2dof[i]);
        }
        _map_var2dof.clear();
    }

    //---------------------------------------------------------------------------------------------
    // setup DoFs
    //---------------------------------------------------------------------------------------------
    /// add DoFs
    /// @param var_id       variable type
    /// @param mesh_id      mesh id
    /// @param list_dof_pt_id   a list of discrete point id defined as DoFs
    size_t addVariableDoFs(size_t mesh_id, const std::vector<size_t> &list_dof_pt_id)
    {
        RandomEquationIdStorage* address = new RandomEquationIdStorage(list_dof_pt_id);
        size_t var_id = _map_var2dof.size();
        _map_var2dof.resize(_map_var2dof.size()+1);
        _map_var2dof[var_id][mesh_id] = address;
        return var_id;
    }

    /// add DoFs
    /// @param var_id           variable type
    /// @param mesh_id          mesh id
    /// @param dof_pt_id_begin  the beginning id of discrete points
    /// @param dof_pt_count     the number of discrete points
    size_t addVariableDoFs(size_t mesh_id, size_t dof_pt_id_begin, size_t dof_pt_count)
    {
        RandomEquationIdStorage* address = new RandomEquationIdStorage(dof_pt_id_begin, dof_pt_count);
        size_t var_id = _map_var2dof.size();
        _map_var2dof.resize(_map_var2dof.size()+1);
        _map_var2dof[var_id][mesh_id] = address;
        _map_msh2var[mesh_id].push_back(var_id);
        return var_id;
    }

    /// add DoFs
    void addVariableDoFs(size_t n_var, size_t mesh_id, size_t dof_pt_id_begin, size_t dof_pt_count)
    {
        _is_seq_address = true;
        for (size_t i=0; i<n_var; i++) {
            SequentialEquationIdStorage* address = new SequentialEquationIdStorage(dof_pt_id_begin, dof_pt_count);
            size_t var_id = _map_var2dof.size();
            _map_var2dof.resize(_map_var2dof.size()+1);
            _map_var2dof[var_id][mesh_id] = address;
            _map_msh2var[mesh_id].push_back(var_id);
        }
    }

    /// deactivate DoFs
    void deactivateDoFs(size_t var_id, size_t mesh_id, size_t pt_id)
    {
        getPointEquationIdTable(var_id, mesh_id)->activate(pt_id, false);
    }

    /// set ghost points for overlapped regions appeared with domain-decomposition methods
    void setGhostPoints(size_t mesh_id, std::vector<size_t> &list_pt_id)
    {
        _ghost_pt[mesh_id].assign(list_pt_id.begin(), list_pt_id.end());
    }

    //---------------------------------------------------------------------------------------------
    // construct the table
    //---------------------------------------------------------------------------------------------
    /// numbering DoFs
    virtual void construct(DofNumberingType::type num=DofNumberingType::BY_DOF, long offset=0);

    //---------------------------------------------------------------------------------------------
    // get summary of the table
    //---------------------------------------------------------------------------------------------
    //# Mesh
    /// get the number of meshes
    size_t getNumberOfMeshes() const {return _map_msh2var.size();};

    //# Variables
    /// get the number of registered variables
    size_t getNumberOfVariables() const { return _map_var2dof.size(); }
    /// get the number of variables associated with the given mesh
    size_t getNumberOfVariables(size_t mesh_id) {return _map_msh2var[mesh_id].size();};

    //# Ghost 
    /// is this point a ghost?
    bool isGhostPoint(size_t mesh_id, size_t pt_id) const
    {
        if (_ghost_pt.size()==0) return false;
        assert(_ghost_pt.count(mesh_id)>0);
        const std::vector<size_t>& ghosts = _ghost_pt.find(mesh_id)->second;
        return std::find(ghosts.begin(), ghosts.end(), pt_id)!=ghosts.end();
    }
    /// get the number of ghost points
    size_t getNumberOfGhostPoints(size_t mesh_id) {return _ghost_pt[mesh_id].size();}
    /// get the total number of ghost points
    size_t getTotalNumberOfGhostPoints() const
    {
        size_t count = 0;
        for (std::map<size_t, std::vector<size_t> >::const_iterator itr=_ghost_pt.begin(); itr!=_ghost_pt.end(); ++itr) {
            count += itr->second.size();
        }
        return count;
    }

    //# DoFs
    /// is this DoF active?
    bool isActiveDoF(size_t var_id, size_t mesh_id, size_t pt_id) const
    {
        return getPointEquationIdTable(var_id, mesh_id)->isActive(pt_id);
    }
    /// get the total number of active DoFs
    size_t getTotalNumberOfActiveDoFs() const { return _total_dofs; }
    /// get the total number of active DoFs excluding ghost points
    size_t getTotalNumberOfActiveDoFsWithoutGhost() const {return _total_dofs_without_ghosts;};


    //---------------------------------------------------------------------------------------------
    // equation address mapping
    //---------------------------------------------------------------------------------------------
    /// map dof to equation id
    size_t mapEqsID(size_t var_id, size_t mesh_id, size_t pt_id) const
    {
        const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
        return add->address(pt_id);
    }
    void mapEqsID(size_t var_id, size_t mesh_id, const std::vector<size_t> &pt_id, std::vector<size_t> &eqs_id) const
    {
        eqs_id.resize(pt_id.size());
        const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
        for (size_t i=0; i<pt_id.size(); i++)
            eqs_id[i] = add->address(pt_id[i]);
    }
    void mapEqsID(size_t mesh_id, const std::vector<size_t> &pt_id, std::vector<size_t> &eqs_id, DofNumberingType::type numbering = DofNumberingType::BY_DOF) const
    {
        if (_map_msh2var.count(mesh_id)==0) return;

        const std::vector<size_t> &list_var = _map_msh2var.find(mesh_id)->second;
        const size_t n_pt = pt_id.size();
        const size_t n_dof_per_pt = list_var.size();
        const size_t n_total_dof = n_pt * n_dof_per_pt;
        eqs_id.resize(n_total_dof);

        if (numbering==DofNumberingType::BY_DOF) {
            for (size_t i=0; i<list_var.size(); i++) {
                size_t var_id = list_var[i];
                const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
                for (size_t j=0; j<pt_id.size(); j++) {
                    eqs_id[j+i*n_pt] = add->address(pt_id[j]);
                }
            }
        } else {
            for (size_t i=0; i<list_var.size(); i++) {
                size_t var_id = list_var[i];
                const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
                for (size_t j=0; j<pt_id.size(); j++) {
                    eqs_id[i + j*n_dof_per_pt] = add->address(pt_id[j]);
                }
            }
        }
    }
    void mapEqsID(size_t mesh_id, const std::vector<size_t> &pt_id, std::vector<size_t> &eqs_id, std::vector<size_t> &eqs_id_without_ghost) const
    {
        if (_map_msh2var.count(mesh_id)==0) return;

        const std::vector<size_t> &list_var = _map_msh2var.find(mesh_id)->second;
        eqs_id.resize(list_var.size()*pt_id.size());
        eqs_id_without_ghost.resize(list_var.size()*pt_id.size());
        for (size_t i=0; i<list_var.size(); i++) {
            size_t var_id = list_var[i];
            const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
            for (size_t i=0; i<pt_id.size(); i++) {
                eqs_id[i] = add->address(pt_id[i]);
                eqs_id_without_ghost[i] = isGhostPoint(mesh_id, pt_id[i]) ? BaseLib::index_npos : eqs_id[i];
            }
        }
    }
    void mapDoF(size_t eqs_id, size_t &var_id, size_t &mesh_id, size_t &pt_id) const
    {
        var_id = BaseLib::index_npos;
        mesh_id = BaseLib::index_npos;
        pt_id = BaseLib::index_npos;
        for (size_t i=0; i<_map_var2dof.size(); i++) {
            const std::map<size_t, IEquationIdStorage*> &obj = _map_var2dof[i];
            for (std::map<size_t, IEquationIdStorage*>::const_iterator itr=obj.begin(); itr!=obj.end(); ++itr) {
                const IEquationIdStorage* pt2dof = itr->second;
                if (!pt2dof->hasValue(eqs_id)) continue;
                pt_id = pt2dof->key(eqs_id);
                var_id = i;
                mesh_id = itr->first;
                return;
            }
        }
    }

    const IEquationIdStorage* getPointEquationIdTable(size_t var_id, size_t mesh_id) const 
    {
        std::map<size_t, IEquationIdStorage*>::const_iterator itr = _map_var2dof[var_id].find(mesh_id);
        if (itr!=_map_var2dof[var_id].end()) {
            return itr->second;
        } else {
            return 0;
        }
    }
    IEquationIdStorage* getPointEquationIdTable(size_t var_id, size_t mesh_id) { return _map_var2dof[var_id][mesh_id]; }

protected:
    void setTotalDoFs(size_t n) {_total_dofs = n;};
    //std::map<size_t, std::vector<DofMap*> >& getMapMsh2Dof() {return _map_msh2dof;}
    /// get the DoF mapping for the given variable

private:
    DISALLOW_COPY_AND_ASSIGN(DofEquationIdTable);


private:
    std::vector<std::map<size_t, IEquationIdStorage*> > _map_var2dof; // var -> (mesh, address)
    //std::vector<Base::IMappedAddress*> _map_var2dof;
    std::map<size_t, std::vector<size_t> > _map_msh2var; // mesh -> (var)
    std::map<size_t, std::vector<size_t> > _ghost_pt; // mesh -> ghost points

    size_t _total_dofs;
    size_t _total_dofs_without_ghosts;
    bool _is_seq_address;
};

} //end
