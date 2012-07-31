
#pragma once

#include <vector>

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "NumLib/Function/IFunction.h"
#include "NumLib/TransientAssembler/ElementWiseTransientResidualAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "FemLib/Function/FemFunction.h"

#include "FemVariable.h"
#include "SolutionLib/DataType.h"

namespace SolutionLib
{

/**
 * \brief Template class for FEM residual functions
 *
 * \tparam T_LOCAL_RESIDUAL_ASSEMBLER
 */
template <
    class T_LOCAL_RESIDUAL_ASSEMBLER
    >
class TemplateTransientResidualFEMFunction
    : public NumLib::TemplateFunction<SolutionVector, SolutionVector>
{
public:
    typedef T_LOCAL_RESIDUAL_ASSEMBLER UserLocalResidualAssembler;

    /// constructor
    /// @param problem        Fem problem
    /// @param linear_eqs    Discrete linear equation
    TemplateTransientResidualFEMFunction(MeshLib::IMesh* msh, const std::vector<FemVariable*> &list_var, DiscreteLib::DofEquationIdTable* dofManager, UserLocalResidualAssembler* asssembler)
        : _local_assembler(asssembler), _dofManager(dofManager),
          _t_n1(0), _u_n0(0), _st(0), _list_var(list_var), _msh(msh)
    {
    };

    ///
    virtual ~TemplateTransientResidualFEMFunction() {};

    ///
    NumLib::TemplateFunction<SolutionVector,SolutionVector>* clone() const
    {
        return new TemplateTransientResidualFEMFunction
                    <
                        UserLocalResidualAssembler
                    >(_msh, _list_var, _dofManager, _local_assembler);
    }

    /// reset property
    void reset(const NumLib::TimeStep* t, SolutionVector* u_n0, SolutionVector* st)
    {
        this->_t_n1 = const_cast<NumLib::TimeStep*>(t);
        this->_u_n0 = u_n0;
        this->_st = st;
    };

    /// evaluate residual
    /// @param u_n1    current results
    /// @param r residual
    void eval(const SolutionVector &u_n1, SolutionVector &r)
    {
        // input, output
        const NumLib::TimeStep &t_n1 = *this->_t_n1;
        SolutionVector *u_n = this->_u_n0;

#if 0
        // setup BC1 for solution increment
        for (size_t i=0; i<_list_var.size(); i++) {
            FemVariable* var = _list_var[i];
            std::vector<size_t> var_bc_id;
            std::vector<double> var_bc_val;
            for (size_t j=0; j<var->getNumberOfDirichletBC(); j++) {
                FemDirichletBC* bc1 = var->getDirichletBC(j);
                bc1->setup();
                std::vector<double> bc_value_for_dx(bc1->getListOfBCNodes().size(), .0);
                for (size_t k=0; k<bc_value_for_dx.size(); k++) {
                    size_t id = bc1->getListOfBCNodes()[k];
                    double v =  bc1->getListOfBCValues()[k];
                    (*u_n)[id] = v; 
                }
            }
        }
#endif

        //TODO temporally
        std::vector<SolutionVector*> vec_un;
        vec_un.push_back(const_cast<SolutionVector*>(u_n));
        std::vector<SolutionVector*> vec_un1;
        vec_un1.push_back(const_cast<SolutionVector*>(&u_n1));

        // assembly
        NumLib::ElementWiseTransientResidualAssembler assembler(&t_n1, &vec_un, &vec_un1, _local_assembler);
        r = .0;
        assembler.assembly(*_msh, *_dofManager, r);

        // add source/sink term (i.e. RHS in linear equation)
        if (_st!=0)
            r += *_st;

    }

private:
    UserLocalResidualAssembler *_local_assembler;
    DiscreteLib::DofEquationIdTable* _dofManager;
    NumLib::TimeStep* _t_n1;
    SolutionVector* _u_n0;
    SolutionVector* _st;
    std::vector<FemVariable*> _list_var;
    MeshLib::IMesh* _msh;
};


} //end
