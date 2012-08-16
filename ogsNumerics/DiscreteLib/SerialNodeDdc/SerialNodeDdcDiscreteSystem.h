/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NodeDDCSerialSharedDiscreteSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/DDC/DecomposedDomain.h"
#include "DiscreteLib/DDC/DecomposedMesh.h"
#include "DecomposedVector.h"
#include "SerialNodeDdcSharedLinearEquation.h"
#include "SerialNodeDdcDistributedLinearEquation.h"

namespace DiscreteLib
{

/**
 * \brief template class for serial node-based decomposition system
 * 
 * \tparam T_EQS Linear equation class
 */
template <template <typename, typename> class T_EQS>
class TemplateSerialNodeDdcDiscreteSystem : public IDiscreteSystem
{
public:
    template<typename T_LINEAR_SOLVER, typename T_SPARSITY_BUILDER>
    struct MyLinearEquation
    {
        typedef T_EQS<T_LINEAR_SOLVER, T_SPARSITY_BUILDER> type;
    };
    template<typename T>
    struct MyVector
    {
        typedef DecomposedVector<T> type;
    };

    /// \param ddc_global   DDC object
    explicit TemplateSerialNodeDdcDiscreteSystem(DecomposedDomain* ddc_global)
    : _ddc_global(ddc_global), _msh(NULL)
    { };
    
    /// 
    virtual ~TemplateSerialNodeDdcDiscreteSystem() {};

    ///
    virtual MeshLib::IMesh* getMesh() const { return (MeshLib::IMesh*)_msh; };

    /// create a new linear equation
    /// \tparam T_LINEAR_SOLVER     Linear solver class
    /// \tparam T_SPARSITY_BUILDER  Sparse pattern builder class
    /// \param linear_solver
    /// \param dofManager
    template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    typename MyLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::type* createLinearEquation(T_LINEAR_SOLVER* linear_solver, DofEquationIdTable* dofManager)
    {
        return MyLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::type::createInstance(*this, _ddc_global, linear_solver, dofManager);
    }

    //// create a vector
    /// \tparam T   data type
    /// \param n    vector size
    template<typename T>
    typename MyVector<T>::type* createVector(const size_t &n)
    {
        std::vector<size_t> list_range_begin;
        for (size_t i=0; i<_ddc_global->getNumberOfSubDomains(); i++) {
            size_t cnt = (size_t)((double)n / (double)_ddc_global->getNumberOfSubDomains() * (double)i);
            list_range_begin.push_back(cnt);
        }

        return MyVector<T>::type::createInstance(*this, n, list_range_begin);
    };

private:
    DISALLOW_COPY_AND_ASSIGN(TemplateSerialNodeDdcDiscreteSystem);

private:
    DecomposedDomain* _ddc_global;
    DecomposedMesh* _msh;
};

/// \brief Discrete system for node-based decomposition using shared memory
typedef TemplateSerialNodeDdcDiscreteSystem<SerialNodeDdcSharedLinearEquation> SerialNodeDdcSharedDiscreteSystem;

/// \brief Discrete system for node-based decomposition using distributed memory
typedef TemplateSerialNodeDdcDiscreteSystem<SerialNodeDdcDistributedLinearEquation> SerialNodeDdcDistributedDiscreteSystem;

} //end
