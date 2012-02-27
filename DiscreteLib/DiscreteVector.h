
#pragma once

#include <vector>


namespace DiscreteLib
{

/**
 * \brief Interface of Vector containers in discrete systems
 */
class IDiscreteVector 
{
public:
    virtual ~IDiscreteVector() {};
    //virtual double size() const = 0;
    //virtual double dot(IDiscreteVector &vec) = 0;
    //virtual double norm1() = 0;
    //virtual double norm2() = 0;
    //virtual double norm_max() = 0;

};

/**
 * \brief Vector container for single memory
 */
template<typename T>
class DiscreteVector : public IDiscreteVector
{
public:
    DiscreteVector(size_t n) : _data(n) {};
    virtual ~DiscreteVector() {};

    virtual size_t size() const {return _data.size();};
    virtual double dot(const DiscreteVector &vec) {return .0;};
    virtual double norm1() {return .0;};
    virtual double norm2() {return .0;};
    virtual double norm_max() {return .0;};

    virtual T& operator[] (size_t idx) 
    {
        return _data[idx];
    }
    virtual const T& operator[] (size_t idx) const
    {
        return _data[idx];
    }

    typename std::vector<T>::iterator begin() 
    {
        return _data.begin();
    }
    typename std::vector<T>::iterator end() 
    {
        return _data.end();
    }

protected:
    std::vector<T> _data;
};

template<typename T>
class DecomposedLocalVector : public DiscreteVector<T>
{     
 public:
    DecomposedLocalVector(size_t local_id, size_t n_global, size_t i_begin, size_t i_end) : DiscreteVector(i_end - i_begin) 
    {
        _local_id = local_id;
        _global_n = n_global;
        _i_start = i_begin;
        _i_end = i_end;
    };

    DecomposedLocalVector(size_t local_id, size_t n_global, size_t i_begin, size_t i_end, std::vector<size_t> &ghost_id) : DiscreteVector(i_end - i_begin + ghost_id.size()), _ghost_id(ghost_id)
    {
        _local_id = local_id;
        _global_n = n_global;
        _i_start = i_begin;
        _i_end = i_end;
    };

    virtual ~DecomposedLocalVector() {};

    T& operator[] (size_t i) 
    {
        return _data[i];
    }

    const T& operator[] (size_t i) const
    {
        return _data[i];
    }
    T& global (size_t i) 
    {
        return _data[access(i)];
    }

   const T& global (size_t i) const
   {
         return _data[access(i)];
    }

    size_t getRangeBegin() const {return _i_start;};
    size_t getRangeEnd() const {return _i_end;};

    size_t getNumberOfGhostElements() const {return _ghost_id.size();};
    size_t getGhostElementId(size_t i) const {return _ghost_id[i];};

private:
    size_t _global_n;
    size_t _i_start;
    size_t _i_end;
    size_t _local_id;
    std::vector<size_t> _ghost_id;

    inline size_t access(size_t global_idx) const 
    {
        if (global_idx < _i_end) {
            assert(_i_start<=global_idx && global_idx<_i_end);
            return global_idx - _i_start;
        } else {
            assert(_ghost_id.end()!=std::find(_ghost_id.begin(), _ghost_id.end(), global_idx));
            return _i_end + (std::find(_ghost_id.begin(), _ghost_id.end(), global_idx) - _ghost_id.begin());
        }
    };
};

/**
 * \brief Vector container for distributed memory
 */
template<typename T>
class DecomposedMasterVector : public IDiscreteVector
{
public:
    DecomposedMasterVector(size_t n_global)
    {
        _global_n = n_global;
    };

    virtual ~DecomposedMasterVector() {};

    void decompose(size_t n)
    {
        _local_v.resize(n);
        _local_v_begin.resize(n);
    }

    DecomposedLocalVector<T>* createLocal(size_t local_id, size_t i_begin, size_t i_end, std::vector<size_t> &ghost_id)
    {

    }

    DecomposedLocalVector<T>* createLocal(size_t local_id, size_t i_begin, size_t i_end)
    {
        DecomposedLocalVector<T> *v = new DecomposedLocalVector<T>(local_id, _global_n, i_begin, i_end);
        _local_v[local_id] = v;
        _local_v_begin[local_id] = i_begin;
        return v;
    }

    void scatter()
    {
        for (size_t i=0; i<_local_v.size(); i++) {
            DecomposedLocalVector<T> *v = _local_v[i];
            for (size_t j=0; j<v->getNumberOfGhostElements(); i++) {
                size_t id = v->getGhostElementId(j);
                v->global(id) = (*this)[id];
            }
        }
    }

    size_t size() const
    {
        return _global_n;
    }

    T& operator[] (size_t idx) 
    {
        return local(idx).global(idx);
    }
    const T& operator[] (size_t idx) const
    {
        return local(idx).global(idx);
    }

private:
    size_t _global_n;
    std::vector<size_t> _local_v_begin;
    std::vector<DecomposedLocalVector<T>*> _local_v;

    DecomposedLocalVector<T>& local(size_t idx)
    {
        return *_local_v[localId(idx)];
    };

    const DecomposedLocalVector<T>& local(size_t idx) const
    {
        return *_local_v[localId(idx)];
    };

    size_t localId(size_t idx) const
    {
        size_t pos = 0;
        for (size_t i=0; i<_local_v_begin.size(); ++i) {
            pos = i;
            if (idx < _local_v_begin[i]) {
                pos = i - 1;
                break;
            }
        }
        return pos;
    }

};


} // end
