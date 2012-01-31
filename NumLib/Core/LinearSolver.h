
#pragma once

#include "LinearEquation.h"

namespace NumLib
{
template<typename Tmat, typename Tvec>
class LinearSolver
{
public:

    void initialize() 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    void solve( LinearEquation<Tmat, Tvec> * _linearEQS ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }



};
}
