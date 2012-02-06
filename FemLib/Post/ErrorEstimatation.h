
#pragma once

namespace FemLib
{
/**
 * \brief Error estimator
 */    
class IErrorEstimator {};

/**
 * \brief Zienkiewicz-Zhu error estimator 
 */
class FemErrorZienkiewiczZhu : IErrorEstimator
{
public:
    void estimate(double* duh, double* dwh)
    {

    }
};

}
