
#pragma once

#include "Base/BidirectionalMap.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteSystem.h"
#include "DomainDecomposition.h"
#include "MeshBasedDiscreteLinearEquation.h"
#include "DDCDiscreteVector.h"

namespace DiscreteLib
{

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
class NodeDDCSerialSharedLinearEquation : public IDiscreteLinearEquation
{
public:
    NodeDDCSerialSharedLinearEquation(DDCGlobal &ddc_global, T_LINEAR_SOLVER &shared_linear_solver, DofEquationIdTable &global_dofManager)
    {
        _ddc_global = &ddc_global;
        for (size_t i=0; i<ddc_global.getNumberOfSubDomains(); i++) {
            DDCSubDomain* dom = ddc_global.getSubDomain(i);
            DofEquationIdTable *local_dofManager = new DofEquationIdTable();
            _list_local_eq.push_back(new TemplateMeshBasedDiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*dom->getLoalMesh(), shared_linear_solver, *local_dofManager));
        }
        _do_create_eqs = true;
        _sheared_eqs = &shared_linear_solver;
        _global_dofManager = &global_dofManager;
    };

    virtual ~NodeDDCSerialSharedLinearEquation() 
    {
        Base::releaseObjectsInStdVector(_list_local_eq);
    };

    /// initialize EQS
    void initialize()
    {
        DofEquationIdTable* global_dofManager = getDofMapManger();
        if (_do_create_eqs) {
            _do_create_eqs = false;
            //create global linear equation
            MathLib::RowMajorSparsity global_sparse;
            SparsityBuilderFromDDC<T_SPARSITY_BUILDER> sp_builder(*_ddc_global, *global_dofManager, global_sparse);
            _sheared_eqs->create(global_dofManager->getTotalNumberOfActiveDoFs(), &global_sparse);
            //local dof
            size_t mesh_id = _ddc_global->getID();
            for (size_t i=0; i<_list_local_eq.size(); i++) {
                DDCSubDomain* dom = _ddc_global->getSubDomain(i);
                MeshLib::IMesh* local_msh = dom->getLoalMesh();
                DofEquationIdTable *local_dofManager = _list_local_eq[i]->getDofMapManger();
                // for each var
                for (size_t j=0; j<global_dofManager->getNumberOfVariables(); j++) {
                    local_dofManager->addVariableDoFs(mesh_id, 0, local_msh->getNumberOfNodes());
                    IMappedAddress* local_pt2dof = local_dofManager->getPointEquationIdTable(j, mesh_id);
                    for (size_t k=0; k<local_msh->getNumberOfNodes(); k++) {
                        size_t dof_id = global_dofManager->mapEqsID(mesh_id, j, dom->getGlobalLocalIdMap()->local2global(k));
                        local_pt2dof->set(k, dof_id);
                    }
                }
                // ghost
                if (dom->getGhostList()!=0) {
                    std::set<size_t>* dom_ghosts = dom->getGhostList();
                    std::vector<size_t> vec_ghosts(dom_ghosts->begin(), dom_ghosts->end());
                    local_dofManager->setGhostPoints(mesh_id, vec_ghosts);
                }
            }
        } else {
            _sheared_eqs->reset();
        }
    }

    /// set prescribed dof
    void setPrescribedDoF(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_prescribed_values)
    {
        //assert(_list_prescribed_dof_id.size()==0);
        _list_prescribed_dof_id.clear();
        _list_prescribed_values.clear();

        const size_t n = list_discrete_pt_id.size();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (_global_dofManager->isActiveDoF(dofId, 0, pt_id)) {
                _list_prescribed_dof_id.push_back(_global_dofManager->mapEqsID(dofId, 0, pt_id));
                _list_prescribed_values.push_back(list_prescribed_values[i]);
            }
        }
    }

    /// construct 
    void construct(IDiscreteLinearEquationAssembler& assemler) 
    {
        DofEquationIdTable* dofManager = getDofMapManger();

        // local eqs
        for (size_t i=0; i<_list_local_eq.size(); i++) {
            _list_local_eq[i]->construct(assemler);
            //_global_eqs->printout();
        }

        //apply 1st bc
        _sheared_eqs->setKnownX(_list_prescribed_dof_id, _list_prescribed_values);
        //_global_eqs->printout();
    }

    /// solve
    void solve()
    {
        _sheared_eqs->solve();
    }

    /// get solution
    double* getLocalX()
    {
        return _sheared_eqs->getX();
    }

    void getGlobalX(std::vector<double> &x)
    {
        x.resize(_sheared_eqs->getDimension());
        double *tmp_x = _sheared_eqs->getX();
        for (size_t i=0; i<x.size(); i++)
            x[i] = tmp_x[i];
    }
    DiscreteVector<double>* getX() 
    {
        return 0;
    };
    /// get a Dof map manager
    DofEquationIdTable* getDofMapManger() const
    {
        return _global_dofManager;
    }
    /// set additional RHS values
    void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_rhs_values, double fkt)
    {
        const size_t n = list_discrete_pt_id.size();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (_global_dofManager->isActiveDoF(dofId, 0, pt_id)) {
                _sheared_eqs->addRHS(_global_dofManager->mapEqsID(dofId, 0, pt_id), list_rhs_values[i]*fkt);
            }
        }
    }

private:
    DDCGlobal* _ddc_global;
    std::vector<AbstractMeshBasedDiscreteLinearEquation*> _list_local_eq;
    T_LINEAR_SOLVER* _sheared_eqs;
    bool _do_create_eqs;
    DofEquationIdTable* _global_dofManager;
    std::vector<size_t> _list_prescribed_dof_id;
    std::vector<double> _list_prescribed_values;

    DISALLOW_COPY_AND_ASSIGN(NodeDDCSerialSharedLinearEquation);
};


/**
 * 
 */
class NodeDDCSerialSharedDiscreteSystem : public IDiscreteSystem
{
public:
    NodeDDCSerialSharedDiscreteSystem(DDCGlobal &ddc_global) : _ddc_global(&ddc_global)
    {

    }
    virtual ~NodeDDCSerialSharedDiscreteSystem()
    {

    }

    /// create a new linear equation
    template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    IDiscreteLinearEquation* createLinearEquation(T_LINEAR_SOLVER &linear_solver, DofEquationIdTable &dofManager)
    {
        _vec_linear_sys.push_back(new NodeDDCSerialSharedLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_ddc_global, linear_solver, dofManager));
        //return _vec_linear_sys.size()-1;
        return _vec_linear_sys.back();
    }

    //IDiscreteLinearEquation* getLinearEquation(size_t i)
    //{
    //    return _vec_linear_sys[i];
    //}

    template<typename T>
    IDiscreteVector<T>* createVector(const size_t &n) 
    {
        std::vector<size_t> list_range_begin;
        for (size_t i=0; i<_ddc_global->getNumberOfSubDomains(); i++) {
            size_t cnt = (size_t)((double)n / (double)_ddc_global->getNumberOfSubDomains() * (double)i);
            list_range_begin.push_back(cnt);
        }

        DDCDiscreteVector<T>* v = new DDCDiscreteVector<T>(n, list_range_begin);
        _vec_vectors.push_back(v);
        return v;
    };
private:
    DISALLOW_COPY_AND_ASSIGN(NodeDDCSerialSharedDiscreteSystem);

private:
    DDCGlobal* _ddc_global;

    // linear equations
    std::vector<IDiscreteLinearEquation*> _vec_linear_sys;
    // vector
    std::vector<IDiscreteVectorBase*> _vec_vectors;};

} //end
