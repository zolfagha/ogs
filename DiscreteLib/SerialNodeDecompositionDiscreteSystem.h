
#pragma once

#include "Base/BidirectionalMap.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteSystem.h"
#include "DomainDecomposition.h"

namespace DiscreteLib
{




template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
class SerialDDCLinearEquation : public IDiscreteLinearEquation
{
public:
    SerialDDCLinearEquation(DDCGlobal &ddc_global, T_LINEAR_SOLVER &linear_solver, DofMapManager &dofManager)
    {
        _ddc_global = &ddc_global;
        for (size_t i=0; i<ddc_global.getNumberOfSubDomains(); i++) {
            DDCSubDomain* dom = ddc_global.getSubDomain(i);
            _list_local_eq.push_back(new TemplateMeshBasedDiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*dom->getLoalMesh(), linear_solver, dofManager));
        }
        _do_create_eqs = true;
        _global_eqs = &linear_solver;
        _global_dofManager = &dofManager;
    };

    virtual ~SerialDDCLinearEquation() 
    {
        Base::releaseObjectsInStdVector(_list_local_eq);
    };

    void initialize()
    {
        DofMapManager* dofManager = getDofMapManger();
        if (_do_create_eqs) {
            _do_create_eqs = false;
            //create global linear equation
            MathLib::RowMajorSparsity global_sparse;
            std::vector<MathLib::RowMajorSparsity*> list_local_sparse(_list_local_eq.size());
            for (size_t i=0; i<_list_local_eq.size(); i++) {
                list_local_sparse[i] = _list_local_eq[i]->getSparsity();
            }
            SparsityBuilderFromLocalSparsity sp_builder(list_local_sparse, *dofManager, global_sparse);
            _global_eqs->create(dofManager->getTotalNumberOfActiveDoFs(), &global_sparse);
        } else {
            _global_eqs->reset();
        }
    }

    /// construct 
    void construct(IDiscreteLinearEquationAssembler& assemler) 
    {
        DofMapManager* dofManager = getDofMapManger();

        // local eqs
        for (size_t i=0; i<_list_local_eq.size(); i++)
            _list_local_eq[i]->construct(assemler);

        // global eqs
        for (size_t i=0; i<_list_local_eq.size(); i++) {
            
        }
    }

    /// solve
    void solve()
    {
        _global_eqs->solve();
    }

    /// get solution
    double* getLocalX()
    {
        return _global_eqs->getX();
    }

    void getGlobalX(std::vector<double> &x)
    {
        x.resize(_global_eqs->getDimension());
        double *tmp_x = _global_eqs->getX();
        for (size_t i=0; i<x.size(); i++)
            x[i] = tmp_x[i];
    }
    DiscreteVector<double>* getX() 
    {
        return 0;
    };
    /// get a Dof map manager
    DofMapManager* getDofMapManger() const
    {
        return _global_dofManager;
    }
    /// set prescribed dof
    void setPrescribedDoF(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> list_prescribed_values)
    {
        //assert(_list_prescribed_dof_id.size()==0);
        _list_prescribed_dof_id.clear();
        _list_prescribed_values.clear();

        const DofMap* dofMap = getDofMapManger()->getDofMap(dofId);
        const size_t n = list_discrete_pt_id.size();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (dofMap->isActiveDoF(pt_id)) {
                _list_prescribed_dof_id.push_back(dofMap->getEqsID(pt_id));
                _list_prescribed_values.push_back(list_prescribed_values[i]);
            }
        }
    }
    /// set additional RHS values
    void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> list_rhs_values, double fkt)
    {
        const DofMap* dofMap = getDofMapManger()->getDofMap(dofId);
        const size_t n = list_discrete_pt_id.size();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (dofMap->isActiveDoF(pt_id)) {
                _global_eqs->addRHS(dofMap->getEqsID(pt_id), list_rhs_values[i]*fkt);
            }
        }
    }

private:
    DDCGlobal* _ddc_global;
    std::vector<AbstractMeshBasedDiscreteLinearEquation*> _list_local_eq;
    T_LINEAR_SOLVER* _global_eqs;
    bool _do_create_eqs;
    DofMapManager* _global_dofManager;
    std::vector<size_t> _list_prescribed_dof_id;
    std::vector<double> _list_prescribed_values;

    DISALLOW_COPY_AND_ASSIGN(SerialDDCLinearEquation);
};


/**
 * 
 */
class SerialNodeDecompositionDiscreteSystem : public IDiscreteSystem
{
public:
    SerialNodeDecompositionDiscreteSystem(DDCGlobal &ddc_global) : _ddc_global(&ddc_global)
    {

    }
    virtual ~SerialNodeDecompositionDiscreteSystem()
    {

    }

    /// create a new linear equation
    template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    IDiscreteLinearEquation* createLinearEquation(T_LINEAR_SOLVER &linear_solver, DofMapManager &dofManager)
    {
        _vec_linear_sys.push_back(new SerialDDCLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_ddc_global, linear_solver, dofManager));
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
        DDCDiscreteVector<T>* v = new DDCDiscreteVector<T>(n, _ddc_global->getNumberOfSubDomains());
        _vec_vectors.push_back(v);
        return v;
    };
private:
    DISALLOW_COPY_AND_ASSIGN(SerialNodeDecompositionDiscreteSystem);

private:
    DDCGlobal* _ddc_global;

    // linear equations
    std::vector<IDiscreteLinearEquation*> _vec_linear_sys;
    // vector
    std::vector<IDiscreteVectorBase*> _vec_vectors;};

} //end
