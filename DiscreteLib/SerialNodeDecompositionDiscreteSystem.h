
#pragma once

#include "Base/BidirectionalMap.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteSystem.h"
#include "DomainDecomposition.h"
#include "MeshBasedDiscreteLinearEquation.h"
#include "DDCDiscreteVector.h"

namespace DiscreteLib
{

void createLocalDofManager(const DofEquationIdTable &global, DDCGlobal &ddc_global, DDCSubDomain &dom, DofEquationIdTable &local)
{
    MeshLib::IMesh* local_msh = dom.getLoalMesh();

    for (size_t i=0; i<global.getNumberOfVariables(); i++) {
        //const DofMap* dofmap = global.getVariableDoF(i);
        //local.addVariableDoF(local_msh->getNumberOfNodes(), dom.getGhostList(), 0, dofmap->getOrder(), 0);
    }
    size_t offset = dom.getGlobalLocalIdMap()->local2global(0)*global.getNumberOfVariables();
    local.construct(DofNumberingType::BY_POINT, offset);
}


template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
class DDCSerialSharedLinearEquation : public IDiscreteLinearEquation
{
public:
    DDCSerialSharedLinearEquation(DDCGlobal &ddc_global, T_LINEAR_SOLVER &global_linear_solver, DofEquationIdTable &global_dofManager)
    {
        _ddc_global = &ddc_global;
        for (size_t i=0; i<ddc_global.getNumberOfSubDomains(); i++) {
            DDCSubDomain* dom = ddc_global.getSubDomain(i);
            DofEquationIdTable *local_dofManager = new DofEquationIdTable();
            createLocalDofManager(global_dofManager,  ddc_global, *dom, *local_dofManager);
            _list_local_eq.push_back(new TemplateMeshBasedDiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*dom->getLoalMesh(), global_linear_solver, *local_dofManager));
        }
        _do_create_eqs = true;
        _global_eqs = &global_linear_solver;
        _global_dofManager = &global_dofManager;
    };

    virtual ~DDCSerialSharedLinearEquation() 
    {
        Base::releaseObjectsInStdVector(_list_local_eq);
    };

    /// initialize EQS
    void initialize()
    {
        DofEquationIdTable* dofManager = getDofMapManger();
        if (_do_create_eqs) {
            _do_create_eqs = false;
            //create global linear equation
            MathLib::RowMajorSparsity global_sparse;
            SparsityBuilderFromDDC<T_SPARSITY_BUILDER> sp_builder(*_ddc_global, *dofManager, global_sparse);
            _global_eqs->create(dofManager->getTotalNumberOfActiveDoFs(), &global_sparse);
        } else {
            _global_eqs->reset();
        }
    }

    /// set prescribed dof
    void setPrescribedDoF(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> list_prescribed_values)
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
    DofEquationIdTable* getDofMapManger() const
    {
        return _global_dofManager;
    }
    /// set additional RHS values
    void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> list_rhs_values, double fkt)
    {
        const size_t n = list_discrete_pt_id.size();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (_global_dofManager->isActiveDoF(dofId, 0, pt_id)) {
                _global_eqs->addRHS(_global_dofManager->mapEqsID(dofId, 0, pt_id), list_rhs_values[i]*fkt);
            }
        }
    }

private:
    DDCGlobal* _ddc_global;
    std::vector<AbstractMeshBasedDiscreteLinearEquation*> _list_local_eq;
    T_LINEAR_SOLVER* _global_eqs;
    bool _do_create_eqs;
    DofEquationIdTable* _global_dofManager;
    std::vector<size_t> _list_prescribed_dof_id;
    std::vector<double> _list_prescribed_values;

    DISALLOW_COPY_AND_ASSIGN(DDCSerialSharedLinearEquation);
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
        _vec_linear_sys.push_back(new DDCSerialSharedLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_ddc_global, linear_solver, dofManager));
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
