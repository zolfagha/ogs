
#pragma once

#include <string>

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

class MyConvergenceCheckerFactory
{
public:
    NumLib::IConvergenceCheck* create(const std::string &)
    {
        return new MyConvergenceCheck();
    };
};
