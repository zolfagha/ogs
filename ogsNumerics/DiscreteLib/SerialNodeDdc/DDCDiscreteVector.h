/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DDCDiscreteVector.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "DiscreteLib/Core/IDiscreteVector.h"
#include "DdcLocalDiscreteVector.h"

namespace DiscreteLib
{

/**
 * \brief Discrete vector container for OpenMP
 *
 * This vector container utilizes shared memory systems and makes local vectors keep actual memory
 */
template<typename T>
class DDCDiscreteVector : public IDiscreteVector<T>
{
private:
    /// @param n_global the size of vector
    /// @param n_divide the number of local vectors
    DDCDiscreteVector(size_t n_global, std::vector<size_t> &list_range_start)
    {
        _global_n = n_global;
        _local_v.resize(list_range_start.size());
        for (size_t i=0; i<list_range_start.size(); i++) {
            size_t i_begin = list_range_start[i];
            size_t i_end = i+1<list_range_start.size() ? list_range_start[i+1] : n_global;
            _local_v[i] = new DDCLocalDiscreteVector<T>(i, _global_n, i_begin, i_end);
        }
    };

public:
    static DDCDiscreteVector<T>* createInstance(IDiscreteSystem &sys, size_t n_global, std::vector<size_t> &list_range_start)
    {
        DDCDiscreteVector<T>* v = new DDCDiscreteVector<T>(n_global, list_range_start);
        sys.addVector(v);
        return v;
    }

    /// destructor
    virtual ~DDCDiscreteVector() 
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

    /// get the range begin
    size_t getRangeBegin() const {return 0;};
    /// get the range end
    size_t getRangeEnd() const {return _global_n;};


    ///// create a local vector
    ///// @param local_id     local vector id
    ///// @param i_begin
    ///// @param i_end
    ///// @param ghost_id
    ///// @return Local vector
    //DDCLocalDiscreteVector<T>* createLocal(size_t local_id, size_t i_begin, size_t i_end, std::vector<size_t> &ghost_id)
    //{
    //    DDCLocalDiscreteVector<T> *v = new DDCLocalDiscreteVector<T>(local_id, _global_n, i_begin, i_end, ghost_id);
    //    _local_v[local_id] = v;
    //    return v;
    //}

    ///// create a local vector
    ///// @param local_id     local vector id
    ///// @param i_begin
    ///// @param i_end
    ///// @return Local vector
    //DDCLocalDiscreteVector<T>* createLocal(size_t local_id, size_t i_begin, size_t i_end)
    //{
    //    DDCLocalDiscreteVector<T> *v = new DDCLocalDiscreteVector<T>(local_id, _global_n, i_begin, i_end);
    //    _local_v[local_id] = v;
    //    return v;
    //}

    /// This function should be called whenever local vectors are updated. It update ghost element from other local vectors.
    void finishUpdate()
    {
        for (size_t i=0; i<_local_v.size(); i++) {
            DDCLocalDiscreteVector<T> *v = _local_v[i];
            for (size_t j=0; j<v->getNumberOfGhostElements(); i++) {
                size_t id = v->getGhostElementId(j);
                v->global(id) = (*this)[id];
            }
        }
    }

private:
    /// vector length
    size_t _global_n;
    /// a list of local vectors
    std::vector<DDCLocalDiscreteVector<T>*> _local_v;

    /// return a local vector which contain the given element
    DDCLocalDiscreteVector<T>& local(size_t idx)
    {
        return *_local_v[localId(idx)];
    };

    /// return a local vector which contain the given element
    const DDCLocalDiscreteVector<T>& local(size_t idx) const
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

};


} //end

