
#pragma once

#include "Base/CodingTools.h"

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"

#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/Core/IDiscreteLinearEquation.h"
#include "DiscreteLib/Assembler/IDiscreteLinearEquationAssembler.h"
#include "DiscreteLib/Utils/SparsityBuilder.h"

#include "AbstractMeshBasedDiscreteLinearEquation.h"

namespace DiscreteLib
{

/**
 * \brief Implementation of mesh based discrete linear equation classes combined with the specific linear solver class
 *
 * @tparam T_LINEAR_SOLVER Linear solver class
 */
template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER=SparsityBuilderDummy>
class TemplateMeshBasedDiscreteLinearEquation : public AbstractMeshBasedDiscreteLinearEquation
{
public:
    ///
    TemplateMeshBasedDiscreteLinearEquation(MeshLib::IMesh &msh, T_LINEAR_SOLVER &linear_solver, DofEquationIdTable &dofManager) 
        : AbstractMeshBasedDiscreteLinearEquation(msh, dofManager), _eqs(&linear_solver), _do_create_eqs(true)
    {
    };

    ///// get the linear solver 
    //MathLib::ILinearEquations* getLinearEquation() const
    //{
    //    return _eqs;
    //}

    /// set prescribed dof
    void setPrescribedDoF(size_t varId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_prescribed_values)
    {
        //assert(_list_prescribed_dof_id.size()==0);
        _list_prescribed_dof_id.clear();
        _list_prescribed_values.clear();

        const size_t n = list_discrete_pt_id.size();
        const DofEquationIdTable* dofmap = getDofMapManger(); 
        const size_t msh_id = getMesh()->getID();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (dofmap->isActiveDoF(varId, msh_id, pt_id)) {
                _list_prescribed_dof_id.push_back(dofmap->mapEqsID(varId, msh_id, pt_id));
                _list_prescribed_values.push_back(list_prescribed_values[i]);
            }
        }
    }

    /// initialize
    void initialize()
    {
        DofEquationIdTable* dofManager = getDofMapManger();
        assert(dofManager->getNumberOfVariables()>0);

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
        DofEquationIdTable* dofManager = getDofMapManger();
        assert(dofManager->getNumberOfVariables()>0);

        assemler.assembly(*getMesh(), *dofManager, *_eqs);

        //apply 1st bc
        _eqs->setKnownX(_list_prescribed_dof_id, _list_prescribed_values);
    }

    /// set additional RHS values
    void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_rhs_values, double fkt=1.0)
    {
        const size_t n = list_discrete_pt_id.size();
        const DofEquationIdTable* dofmap = getDofMapManger(); 
        const size_t msh_id = getMesh()->getID();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (dofmap->isActiveDoF(dofId, msh_id, pt_id)) {
                _eqs->addRHS(dofmap->mapEqsID(dofId, msh_id, pt_id), list_rhs_values[i]*fkt);
            }
        }
    }

    virtual void addRHS(const GlobalVectorType &v, double fkt=1.0)
    {
        const size_t n = v.size();
        for (size_t i=0; i<n; i++) {
			_eqs->addRHS(i, v[i]*fkt);
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

    void getX(GlobalVectorType &x)
    {
        double *tmp_x = _eqs->getX();
        for (size_t i=x.getRangeBegin(); i<x.getRangeEnd(); i++)
            x[i] = tmp_x[i];
    };

    ///// get the RHS vector
    //double* getRHS()
    //{
    //    return _eqs->getRHS();
    //}

    T_LINEAR_SOLVER* getLinearSolver() {return _eqs;};
private:
    T_LINEAR_SOLVER* _eqs;
    bool _do_create_eqs;
    std::vector<size_t> _list_prescribed_dof_id;
    std::vector<double> _list_prescribed_values;

    DISALLOW_COPY_AND_ASSIGN(TemplateMeshBasedDiscreteLinearEquation);
};

} //end
