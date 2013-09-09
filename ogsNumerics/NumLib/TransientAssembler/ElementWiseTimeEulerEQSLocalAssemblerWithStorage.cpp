/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTimeEulerEQSLocalAssemblerWithStorage.cpp
 *
 * Created on 2013-09-09 by Haibing Shao
 */

#include "ElementWiseTimeEulerEQSLocalAssemblerWithStorage.h"

namespace NumLib
{

/// assemble a local linear equation for the given element
/// @param time            time step
/// @param e            element
/// @param local_u_n1    guess of current time step value
/// @param local_u_n    previous time step value
/// @param eqs            local algebraic equation
void ElementWiseTimeEulerEQSLocalAssemblerWithStorage::assembly(const TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n, MathLib::LocalEquation &eqs)
{
    const double delta_t = time.getTimeStepSize();
    const size_t n_dof   = eqs.getDimension();
    const double       t = time.getTime(); 

    // if current time is not the same as _stored_t, 
    // then we need to calculate mass and conductance matrix and force vector
    // otherwise we can use the same values as stored
    if ( t != _stored_t )
    {
        _M = MathLib::LocalMatrix::Zero(n_dof, n_dof);
        _K = MathLib::LocalMatrix::Zero(n_dof, n_dof);
        _F = MathLib::LocalVector::Zero(n_dof);
        
        assembleODE(time, e, local_u_n1, local_u_n, _M, _K, _F);
        
        //std::cout << "M="; M.write(std::cout); std::cout << std::endl;
        //std::cout << "K="; K.write(std::cout); std::cout << std::endl;
        _stored_t = t; // TODO: this will not be for each element. 
    }

    MathLib::LocalMatrix *localA = eqs.getA();
    MathLib::LocalVector *localRHS = eqs.getRHSAsVec();

    // A = 1/dt M + theta K
    (*localA) = 1.0/delta_t * _M + _theta * _K;
    // RHS = (1/dt M - (1-theta) K) u0 + F
    (*localRHS) = (1.0/delta_t * _M - (1.-_theta) * _K) * local_u_n;
    (*localRHS) += _F;
}


} //end
