
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"

#include "MeshLib/Core/IMesh.h"

#include "NumLib/Discrete/DoF.h"
#include "NumLib/Discrete/DiscreteLinearEquationAssembler.h"

namespace NumLib
{

/** 
 * \brief Interface of discrete linear equations
 */
class IDiscreteLinearEquation
{
public:
    /// construct 
    virtual void construct(IDiscreteLinearEquationAssembler& assemler) = 0;
    /// solve
    virtual void solve() = 0;
    /// get solution
    virtual double* getX() = 0;
    /// get RHS 
    virtual double* getRHS() = 0;
    /// get a linear equation object
    virtual MathLib::ILinearEquations* getLinearEquation() const = 0;
    /// get a Dof map manager
    virtual DofMapManager* getDofMapManger() const = 0;

};

/**
 * \brief Abstract class for mesh-based discrete linear equations
 */
class AbstractMeshBasedDiscreteLinearEquation : public IDiscreteLinearEquation
{
public:
    ///
    AbstractMeshBasedDiscreteLinearEquation(MeshLib::IMesh &msh) : _msh(&msh)
    {
        _dofManager = new DofMapManager();
    }

    virtual ~AbstractMeshBasedDiscreteLinearEquation()
    {
        Base::releaseObject(_dofManager);
    }

    MeshLib::IMesh* getMesh() const
    {
        return _msh;
    }

    DofMapManager* getDofMapManger() const 
    {
        return _dofManager;
    }

private:
    MeshLib::IMesh* _msh;
    DofMapManager* _dofManager;

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
    TemplateMeshBasedDiscreteLinearEquation(MeshLib::IMesh &msh, T_LINEAR_SOLVER &linear_solver) 
        : AbstractMeshBasedDiscreteLinearEquation(msh), _eqs(&linear_solver), _do_create_eqs(true)
    {
    };

    /// get the linear solver 
    MathLib::ILinearEquations* getLinearEquation() const
    {
        return _eqs;
    }

    /// construct the linear equation
    void construct(IDiscreteLinearEquationAssembler& assemler)
    {
        DofMapManager* dofManager = getDofMapManger();
        assert(dofManager->getNumberOfDof()>0);

        if (_do_create_eqs) {
            _do_create_eqs = false;
            MathLib::RowMajorSparsity sparse;
            T_SPARSITY_BUILDER sp_builder(*getMesh(), *dofManager, sparse);
            getLinearEquation()->create(dofManager->getTotalNumberOfDiscretePoints(), &sparse);
        } else {
            getLinearEquation()->reset();
        }
        assemler.assembly(*getMesh(), *dofManager, *getLinearEquation());
    }

    /// solve
    void solve()
    {
        _eqs->solve();
    }


    /// get the solution vector
    double* getX()
    {
        return _eqs->getX();
    }

    /// get the RHS vector
    double* getRHS()
    {
        return _eqs->getRHS();
    }
private:
    T_LINEAR_SOLVER* _eqs;
    bool _do_create_eqs;

    DISALLOW_COPY_AND_ASSIGN(TemplateMeshBasedDiscreteLinearEquation);
};


}
