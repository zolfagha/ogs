/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OMPGlobalDiscreteVector.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"

#include "DiscreteLib/Core/IDiscreteResource.h"
#include "OMPLocalDiscreteVector.h"

namespace DiscreteLib
{

/**
 * \brief Discrete vector container for OpenMP
 *
 * This vector container utilizes shared memory systems and makes local vectors keep actual memory
 */
template<typename T>
class OMPGlobalDiscreteVector : public IDiscreteResource
{
public:
    /// @param n_global the size of vector
    /// @param n_divide the number of local vectors
    OMPGlobalDiscreteVector(size_t n_global, size_t n_divide)
    {
        _global_n = n_global;
        _local_v.resize(n_divide);
    };

    /// destructor
    virtual ~OMPGlobalDiscreteVector() 
    {
        BaseLib::releaseObjectsInStdVector(_local_v);
    };

    /// get the size of this vector
    size_t size() const { return _global_n; }

    /// access data
    T& operator[] (size_t idx) 
    {
        return local(idx).global(idx);
    }

    /// access data
    const T& operator[] (size_t idx) const
    {
        return local(idx).global(idx);
    }

    /// create a local vector
    /// @param local_id     local vector id
    /// @param i_begin
    /// @param i_end
    /// @param ghost_id
    /// @return Local vector
    OMPLocalDiscreteVector<T>* createLocal(size_t local_id, size_t i_begin, size_t i_end, std::vector<size_t> &ghost_id)
    {
        OMPLocalDiscreteVector<T> *v = new OMPLocalDiscreteVector<T>(local_id, _global_n, i_begin, i_end, ghost_id);
        _local_v[local_id] = v;
        return v;
    }

    /// create a local vector
    /// @param local_id     local vector id
    /// @param i_begin
    /// @param i_end
    /// @return Local vector
    OMPLocalDiscreteVector<T>* createLocal(size_t local_id, size_t i_begin, size_t i_end)
    {
        OMPLocalDiscreteVector<T> *v = new OMPLocalDiscreteVector<T>(local_id, _global_n, i_begin, i_end);
        _local_v[local_id] = v;
        return v;
    }

    /// This function should be called whenever local vectors are updated. It update ghost element from other local vectors.
    void finishUpdate()
    {
        for (size_t i=0; i<_local_v.size(); i++) {
            OMPLocalDiscreteVector<T> *v = _local_v[i];
            for (size_t j=0; j<v->getNumberOfGhostElements(); i++) {
                size_t id = v->getGhostElementId(j);
                v->global(id) = (*this)[id];
            }
        }
    }

private:
    /// return a local vector which contain the given element
    OMPLocalDiscreteVector<T>& local(size_t idx)
    {
        return *_local_v[localId(idx)];
    };

    /// return a local vector which contain the given element
    const OMPLocalDiscreteVector<T>& local(size_t idx) const
    {
        return *_local_v[localId(idx)];
    };

    /// return an index of a local vector which contain the given element
    size_t localId(size_t idx) const
    {
        size_t pos = 0;
        for (size_t i=0; i<_local_v.size(); ++i) {
            pos = i;
            if (idx < _local_v[i]->getRangeBegin()) {
                pos = i - 1;
                break;
            }
        }
        return pos;
    }

private:
    /// vector length
    size_t _global_n;
    /// a list of local vectors
    std::vector<OMPLocalDiscreteVector<T>*> _local_v;
};



} //end
