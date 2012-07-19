
#pragma once

#include <string>

#include "FemLib/Function/FemNorm.h"
#include "NumLib/Coupling/Algorithm/IConvergenceCheck.h"
#include "NumLib/Function/FunctionConstant.h"

class MyConvergenceCheck : public NumLib::IConvergenceCheck
{
public:
    typedef NumLib::FunctionConstant<double,double> MyFunction;

    virtual ~MyConvergenceCheck() {};

    virtual bool isConverged(NumLib::UnnamedParameterSet& vars_prev, NumLib::UnnamedParameterSet& vars_current, double eps, double &v_diff)
    {
        for (size_t i=0; i<vars_prev.size(); i++) {
            double v_prev = .0;
            vars_prev.get<MyFunction>(i)->eval(v_prev);
            double v_cur = .0;
            vars_current.get<MyFunction>(i)->eval(v_cur);
            v_diff = std::abs(v_cur - v_prev);
            if (v_diff>eps) {
                return false;
            }
        }
        return true;
    }
};

class FemFunctionConvergenceCheck : public NumLib::IConvergenceCheck
{
    //FemLib::NormOfFemNodalFunction<double> _norm;
public:
    //explicit FemFunctionConvergenceCheck(DiscreteLib::DiscreteSystem *dis) : _norm(dis)
    //{

    //}
    FemFunctionConvergenceCheck()
    {

    }

    bool isConverged(NumLib::UnnamedParameterSet& vars_prev, NumLib::UnnamedParameterSet& vars_current, double eps, double &v_diff)
    {

        for (size_t i=0; i<vars_prev.size(); i++) {
#if 1
            if (vars_prev.getName(i).compare("h")==0 || vars_prev.getName(i).compare("c")==0) {
                const FemLib::FemNodalFunctionScalar* f_fem_prev = vars_prev.get<FemLib::FemNodalFunctionScalar>(i);
                const FemLib::FemNodalFunctionScalar* f_fem_cur = vars_current.get<FemLib::FemNodalFunctionScalar>(i);
                //v_diff = f_fem_cur->norm_diff(*f_fem_prev);
                FemLib::NormOfFemNodalFunction<double> _norm(f_fem_cur->getDiscreteSystem());
                v_diff = _norm(*f_fem_prev, *f_fem_cur);
            } else if (vars_prev.getName(i).compare("v")==0) {
                const FemLib::FEMIntegrationPointFunctionVector* f_fem_prev = vars_prev.get<FemLib::FEMIntegrationPointFunctionVector>(i);
                const FemLib::FEMIntegrationPointFunctionVector* f_fem_cur = vars_current.get<FemLib::FEMIntegrationPointFunctionVector>(i);
                //v_diff = f_fem_cur->norm_diff(*f_fem_prev);
                FemLib::NormOfFemNodalFunction<double> _norm(f_fem_cur->getDiscreteSystem());
                v_diff = _norm(*f_fem_prev, *f_fem_cur);
            }
#endif
            if (v_diff>eps) {
                return false;
            }
        }
        return true;
    }
};

class MyConvergenceCheckerFactory
{
public:
    NumLib::IConvergenceCheck* create(const std::string &name)
    {
        if (name.compare("MyConvergenceCheck")==0) {
            return new MyConvergenceCheck();
        } else if (name.compare("FemFunctionConvergenceCheck")==0) {
            return new FemFunctionConvergenceCheck();
        }
        return NULL;
    };
};
