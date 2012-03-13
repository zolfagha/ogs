
#pragma once

#include <vector>

#include "Base/CodingTools.h"

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"

#include "MeshLib/Core/IMesh.h"

#include "DoF.h"
#include "DiscreteLinearEquationAssembler.h"

namespace DiscreteLib
{

/** 
 * \brief Interface of discrete linear equations
 */
class IDiscreteLinearEquation
{
public:
    /// 
    virtual void initialize() = 0;
    /// construct 
    virtual void construct(IDiscreteLinearEquationAssembler& assemler) = 0;
    /// solve
    virtual void solve() = 0;
    /// get solution
    virtual void getGlobalX(std::vector<double> &x) = 0;
    virtual double* getLocalX() = 0;
    DiscreteVector<double>* getX() {return 0;};
//    virtual DiscreteVector<double>* getX() = 0;
    ///// get RHS 
    //virtual double* getRHS() = 0;
    ///// get a linear equation object
    //virtual MathLib::ILinearEquations* getLinearEquation() const = 0;
    /// get a Dof map manager
    virtual DofMapManager* getDofMapManger() const = 0;
    /// set prescribed dof
    virtual void setPrescribedDoF(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> list_prescribed_values) = 0;
    /// set additional RHS values
    virtual void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> list_rhs_values, double fkt) = 0;
};

/**
 * \brief Abstract class for mesh-based discrete linear equations
 */
class AbstractMeshBasedDiscreteLinearEquation : public IDiscreteLinearEquation
{
public:
    ///
    AbstractMeshBasedDiscreteLinearEquation(MeshLib::IMesh &msh, DofMapManager &dofManager) : _msh(&msh), _dofManager(&dofManager), _sparsity(0)
    {
    }

    virtual ~AbstractMeshBasedDiscreteLinearEquation()
    {
        //Base::releaseObject(_dofManager);
        Base::releaseObject(_sparsity);
    }

    MeshLib::IMesh* getMesh() const
    {
        return _msh;
    }

    DofMapManager* getDofMapManger() const 
    {
        return _dofManager;
    }

    MathLib::RowMajorSparsity* getSparsity() const
    {
        return _sparsity;
    }

protected:
    void setSparsity(MathLib::RowMajorSparsity* sp)
    {
        _sparsity = sp;;
    }

private:
    MeshLib::IMesh* _msh;
    DofMapManager* _dofManager;
    MathLib::RowMajorSparsity* _sparsity;

    DISALLOW_COPY_AND_ASSIGN(AbstractMeshBasedDiscreteLinearEquation);
};


/**
 * \brief Implementation of mesh based discrete linear equation classes combined with the specific linear solver class
 *
 * @tparam T_LINEAR_SOLVER Linear solver class
 */
template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
class TemplateMeshBasedDiscreteLinearEquation : public AbstractMeshBasedDiscreteLinearEquation
{
public:
    ///
    TemplateMeshBasedDiscreteLinearEquation(MeshLib::IMesh &msh, T_LINEAR_SOLVER &linear_solver, DofMapManager &dofManager) 
        : AbstractMeshBasedDiscreteLinearEquation(msh, dofManager), _eqs(&linear_solver), _do_create_eqs(true)
    {
    };

    ///// get the linear solver 
    //MathLib::ILinearEquations* getLinearEquation() const
    //{
    //    return _eqs;
    //}

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

    void initialize()
    {
        DofMapManager* dofManager = getDofMapManger();
        assert(dofManager->getNumberOfDof()>0);

        if (_do_create_eqs && !_eqs->isCreated()) {
            _do_create_eqs = false;
            MathLib::RowMajorSparsity* sparse = new MathLib::RowMajorSparsity();
            T_SPARSITY_BUILDER sp_builder(*getMesh(), *dofManager, *sparse);
            AbstractMeshBasedDiscreteLinearEquation::setSparsity(sparse);
            _eqs->create(dofManager->getTotalNumberOfActiveDoFs(), sparse);
        } else {
            _eqs->reset();
        }
    }

    /// construct the linear equation
    void construct(IDiscreteLinearEquationAssembler& assemler)
    {
        DofMapManager* dofManager = getDofMapManger();
        assert(dofManager->getNumberOfDof()>0);

        assemler.assembly(*getMesh(), *dofManager, *_eqs);

        //apply 1st bc
        _eqs->setKnownX(_list_prescribed_dof_id, _list_prescribed_values);
    }

    /// set additional RHS values
    void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> list_rhs_values, double fkt=1.0)
    {
        const DofMap* dofMap = getDofMapManger()->getDofMap(dofId);
        const size_t n = list_discrete_pt_id.size();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (dofMap->isActiveDoF(pt_id)) {
                _eqs->addRHS(dofMap->getEqsID(pt_id), list_rhs_values[i]*fkt);
            }
        }
    }


    /// solve
    void solve()
    {
        _eqs->solve();
    }


    /// get the solution vector
    double* getLocalX()
    {
        return _eqs->getX();
    }

    void getGlobalX(std::vector<double> &x)
    {
        x.resize(_eqs->getDimension());
        double *tmp_x = _eqs->getX();
        for (size_t i=0; i<x.size(); i++)
            x[i] = tmp_x[i];
    }

    ///// get the RHS vector
    //double* getRHS()
    //{
    //    return _eqs->getRHS();
    //}
private:
    T_LINEAR_SOLVER* _eqs;
    bool _do_create_eqs;
    std::vector<size_t> _list_prescribed_dof_id;
    std::vector<double> _list_prescribed_values;

    DISALLOW_COPY_AND_ASSIGN(TemplateMeshBasedDiscreteLinearEquation);
};

}
