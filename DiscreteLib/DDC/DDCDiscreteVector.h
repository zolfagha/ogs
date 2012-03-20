
#pragma once

#include "DiscreteLib/Core/DiscreteVector.h"

namespace DiscreteLib
{
/**
 * \brief Local discrete vector for OpenMP
 */
template<typename T>
class DDCLocalDiscreteVector : public DiscreteVector<T>
{     
 public:
     /// Constructor without ghost elements
     /// @param local_id    id of this local vector
     /// @param n_global    size of the global vector 
     /// @param i_begin     global index starting this local vector
     /// @param i_end       global index ending this local vector + 1 
    DDCLocalDiscreteVector(size_t local_id, size_t n_global, size_t i_begin, size_t i_end) : DiscreteVector<T>(i_end - i_begin)
    {
        _local_id = local_id;
        _global_n = n_global;
        _i_start = i_begin;
        _i_end = i_end;
    };

    /// Constructor with ghost elements
    /// @param local_id    id of this local vector
    /// @param n_global    size of the global vector 
    /// @param i_begin     global index starting this local vector
    /// @param i_end       global index ending this local vector + 1 
    /// @param ghost_id    a list of ghost id
    DDCLocalDiscreteVector(size_t local_id, size_t n_global, size_t i_begin, size_t i_end, std::vector<size_t> &ghost_id) : DiscreteVector<T>(i_end - i_begin + ghost_id.size()), _ghost_id(ghost_id)
    {
        _local_id = local_id;
        _global_n = n_global;
        _i_start = i_begin;
        _i_end = i_end;
    };

    /// destructor
    virtual ~DDCLocalDiscreteVector() {};

    /// get the global vector size
    size_t getGlobalSize() const {return _global_n;};

    /// get the range begin
    size_t getRangeBegin() const {return _i_start;};
    /// get the range end
    size_t getRangeEnd() const {return _i_end;};

    /// get the number of ghost elements in this vector
    size_t getNumberOfGhostElements() const {return _ghost_id.size();};
    /// get the ghost element id (global) in this vector
    size_t getGhostElementId(size_t i) const {return _ghost_id[i];};
    /// check if the give id can be a ghost element in this vector. This function doesn't check if actually this vector contain the element.
    bool isGhost(size_t global_idx) const {return !(_i_start<=global_idx && global_idx<_i_end);};

    /// check if this vector contain an element with the given global id
    bool has(size_t global_idx) const
    {
        if (!isGhost(global_idx)) {
            return (_i_start<=global_idx && global_idx<_i_end);
        } else {
            return (_ghost_id.end()!=std::find(_ghost_id.begin(), _ghost_id.end(), global_idx));
        }
    }

    /// access with local index
    T& operator[] (size_t i_local) { return DiscreteVector<T>::_data[i_local]; }
    /// access with local index
    const T& operator[] (size_t i_local) const { return DiscreteVector<T>::_data[i_local]; }

    /// access with global index
    T& global (size_t i) { return DiscreteVector<T>::_data[access(i)]; }
    /// access with global index
   const T& global (size_t i) const { return DiscreteVector<T>::_data[access(i)]; }


private:
    size_t _global_n;
    size_t _i_start;
    size_t _i_end;
    size_t _local_id;
    /// a list of global id of ghost elements
    std::vector<size_t> _ghost_id;

    /// get local id
    inline size_t access(size_t global_idx) const 
    {
        assert (has(global_idx));
        if (!isGhost(global_idx)) {
            return global_idx - _i_start;
        } else {
            return _i_end + (std::find(_ghost_id.begin(), _ghost_id.end(), global_idx) - _ghost_id.begin());
        }
    };

};

/**
 * \brief Discrete vector container for OpenMP
 *
 * This vector container utilizes shared memory systems and makes local vectors keep actual memory
 */
template<typename T>
class DDCDiscreteVector : public IDiscreteVector<T>
{
public:
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

    /// destructor
    virtual ~DDCDiscreteVector() 
    {
        Base::releaseObjectsInStdVector(_local_v);
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

