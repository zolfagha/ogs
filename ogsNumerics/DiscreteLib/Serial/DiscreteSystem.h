/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"

#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteLinearEquation.h"

#include "DiscreteVector.h"
#include "DiscreteLinearEquation.h"

namespace MeshLib
{
class IMesh;
}

namespace DiscreteLib
{

class DofEquationIdTable;

/**
 * \brief Discrete system based on a single mesh
 */
class DiscreteSystem : public IDiscreteSystem
{
public:
   
    /// constructor
    /// \param msh  a mesh to represent discrete systems by nodes or elements or edges
    explicit DiscreteSystem(MeshLib::IMesh& msh) : _msh(&msh) {};

    ///
    virtual ~DiscreteSystem() {};

    /// get this mesh
    MeshLib::IMesh* getMesh() const { return _msh; };

    /// create a new linear equation
    /// @tparam T_LINEAR_EQUATION
    /// @tparam T_LINEAR_SOLVER
    /// @tparam T_SPARSITY_BUILDER
    /// @param linear_solver         Linear solver
    /// @param dofManager            Equation index table
    template <class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    DiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>* createLinearEquation(T_LINEAR_SOLVER* linear_solver, DofEquationIdTable* dofManager)
    {
        DiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>* eq;
        eq = DiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::createInstance(*this, _msh, linear_solver, dofManager);
        return eq;
    }

    /// create a new vector
    /// @param n    vector length
    /// @return vector object
    template<typename T>
    DiscreteVector<T>* createVector(const size_t n)
    {
        DiscreteVector<T>* v = DiscreteVector<T>::createInstance(*this, n);
        return v;
    };

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteSystem);

protected:
    MeshLib::IMesh* _msh;
};

} //end
