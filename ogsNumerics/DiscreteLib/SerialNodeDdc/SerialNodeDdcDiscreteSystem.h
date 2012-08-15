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

#include "BaseLib/BidirectionalMap.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/LinearEquation/MeshBasedDiscreteLinearEquation.h"
#include "DiscreteLib/DDC/DDCGlobal.h"
#include "DDCDiscreteVector.h"
#include "SerialNodeDdcSharedLinearEquation.h"
#include "SerialNodeDdcDistributedLinearEquation.h"

namespace DiscreteLib
{

/**
 * 
 */
template <template <typename, typename> class T_EQS>
class TemplateSerialNodeDdcDiscreteSystem : public IDiscreteSystem
{
public:

    TemplateSerialNodeDdcDiscreteSystem(DDCGlobal &ddc_global)
    : _ddc_global(&ddc_global)
    {

    }

    virtual ~TemplateSerialNodeDdcDiscreteSystem()
    {

    }

    /// create a new linear equation
    template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    IDiscreteLinearEquation* createLinearEquation(T_LINEAR_SOLVER &linear_solver, DofEquationIdTable &dofManager)
    {
        _vec_linear_sys.push_back(new T_EQS<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_ddc_global, linear_solver, dofManager));
        //return _vec_linear_sys.size()-1;
        return _vec_linear_sys.back();
    }

    //IDiscreteLinearEquation* getLinearEquation(size_t i)
    //{
    //    return _vec_linear_sys[i];
    //}

    template<typename T>
    IDiscreteVector<T>* createVector(const size_t &n) 
    {
        std::vector<size_t> list_range_begin;
        for (size_t i=0; i<_ddc_global->getNumberOfSubDomains(); i++) {
            size_t cnt = (size_t)((double)n / (double)_ddc_global->getNumberOfSubDomains() * (double)i);
            list_range_begin.push_back(cnt);
        }

        DDCDiscreteVector<T>* v = new DDCDiscreteVector<T>(n, list_range_begin);
        _vec_vectors.push_back(v);
        return v;
    };
private:
    DISALLOW_COPY_AND_ASSIGN(TemplateSerialNodeDdcDiscreteSystem);

private:
    DDCGlobal* _ddc_global;

    // linear equations
    std::vector<IDiscreteLinearEquation*> _vec_linear_sys;
    // vector
    std::vector<IDiscreteResource*> _vec_vectors;
};

typedef TemplateSerialNodeDdcDiscreteSystem<SerialNodeDdcSharedLinearEquation> SerialNodeDdcSharedDiscreteSystem;
typedef TemplateSerialNodeDdcDiscreteSystem<SerialNodeDdcDistributedLinearEquation> SerialNodeDdcDistributedDiscreteSystem;

} //end
