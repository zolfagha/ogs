
#pragma once

#include <vector>
#include <cstddef>

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/Core/IDiscreteVector.h"

namespace DiscreteLib
{

/**
 * \brief Vector container for single memory
 */
template<typename T>
class DiscreteVector : public IDiscreteVector<T>
{
protected:
    friend class DiscreteSystem;
    DiscreteVector() {};
    explicit DiscreteVector(size_t n) : _data(n) {};

public:
    virtual ~DiscreteVector() {};

    virtual void resize(size_t n) {_data.resize(n);};
    virtual size_t size() const {return _data.size();};

    virtual DiscreteVector<T>& operator= (T v)
    {
        for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] = v;
        return *this;
    }

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
    virtual size_t getRangeBegin() const
    {
        return 0;
    }
    virtual size_t getRangeEnd() const
    {
        return _data.size();
    }

protected:
    std::vector<T> _data;
};

} // end
