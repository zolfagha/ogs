/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MyConvergenceCheckerFactory.h
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

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

template <class T_DIS_SYS>
class FemFunctionConvergenceCheck : public NumLib::IConvergenceCheck
{
    //FemLib::NormOfFemNodalFunction<double> _norm;
public:
    typedef typename FemLib::FemNodalFunctionScalar<T_DIS_SYS>::type MyNodalFunctionScalar;
    typedef typename FemLib::FEMIntegrationPointFunctionVector<T_DIS_SYS>::type MyIntegrationPointFunctionVector;
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
                const MyNodalFunctionScalar* f_fem_prev = vars_prev.get<MyNodalFunctionScalar>(i);
                const MyNodalFunctionScalar* f_fem_cur = vars_current.get<MyNodalFunctionScalar>(i);
                //v_diff = f_fem_cur->norm_diff(*f_fem_prev);
                FemLib::NormOfFemNodalFunction<T_DIS_SYS, double> _norm(f_fem_cur->getDiscreteSystem());
                v_diff = _norm(*f_fem_prev, *f_fem_cur);
            } else if (vars_prev.getName(i).compare("v")==0) {
                const MyIntegrationPointFunctionVector* f_fem_prev = vars_prev.get<MyIntegrationPointFunctionVector>(i);
                const MyIntegrationPointFunctionVector* f_fem_cur = vars_current.get<MyIntegrationPointFunctionVector>(i);
                //v_diff = f_fem_cur->norm_diff(*f_fem_prev);
                FemLib::NormOfFemNodalFunction<T_DIS_SYS, double> _norm(f_fem_cur->getDiscreteSystem());
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

template <class T_DIS_SYS>
class MyConvergenceCheckerFactory
{
public:
    NumLib::IConvergenceCheck* create(const std::string &name)
    {
        if (name.compare("MyConvergenceCheck")==0) {
            return new MyConvergenceCheck();
        } else if (name.compare("FemFunctionConvergenceCheck")==0) {
            return new FemFunctionConvergenceCheck<T_DIS_SYS>();
        }
        return NULL;
    };
};
