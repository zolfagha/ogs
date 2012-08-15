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
#include "DiscreteLib/DDC/DDCGlobal.h"
#include "DDCDiscreteVector.h"
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
    /// \param ddc_global   DDC object
    explicit TemplateSerialNodeDdcDiscreteSystem(DDCGlobal &ddc_global)
    : _ddc_global(&ddc_global)
    { };
    
    /// 
    virtual ~TemplateSerialNodeDdcDiscreteSystem() {};

    /// create a new linear equation
    /// \tparam T_LINEAR_SOLVER     Linear solver class
    /// \tparam T_SPARSITY_BUILDER  Sparse pattern builder class
    /// \param linear_solver
    /// \param dofManager
    template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    IDiscreteLinearEquation* createLinearEquation(T_LINEAR_SOLVER &linear_solver, DofEquationIdTable &dofManager)
    {
        T_EQS<T_LINEAR_SOLVER, T_SPARSITY_BUILDER> *eqs = new T_EQS<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_ddc_global, linear_solver, dofManager);
        _data_container.addLinearEquation(eqs);
        return eqs;
    }

    //// create a vector
    /// \tparam T   data type
    /// \param n    vector size
    template<typename T>
    IDiscreteVector<T>* createVector(const size_t &n) 
    {
        std::vector<size_t> list_range_begin;
        for (size_t i=0; i<_ddc_global->getNumberOfSubDomains(); i++) {
            size_t cnt = (size_t)((double)n / (double)_ddc_global->getNumberOfSubDomains() * (double)i);
            list_range_begin.push_back(cnt);
        }

        DDCDiscreteVector<T>* v = new DDCDiscreteVector<T>(n, list_range_begin);
        _data_container.addVector(v);
        return v;
    };

private:
    DISALLOW_COPY_AND_ASSIGN(TemplateSerialNodeDdcDiscreteSystem);

private:
    DDCGlobal* _ddc_global;
    DiscreteDataContainer _data_container;
};

/// \brief Discrete system for node-based decomposition using shared memory
typedef TemplateSerialNodeDdcDiscreteSystem<SerialNodeDdcSharedLinearEquation> SerialNodeDdcSharedDiscreteSystem;

/// \brief Discrete system for node-based decomposition using distributed memory
typedef TemplateSerialNodeDdcDiscreteSystem<SerialNodeDdcDistributedLinearEquation> SerialNodeDdcDistributedDiscreteSystem;

} //end
