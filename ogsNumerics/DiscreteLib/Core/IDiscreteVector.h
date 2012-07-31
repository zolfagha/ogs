
#pragma once

#include <vector>
#include <cstddef>

#include "BaseLib/CodingTools.h"
#include "IDiscreteResource.h"

namespace DiscreteLib
{


template<typename T>
class IDiscreteVector : public IDiscreteResource
{
public:
    virtual ~IDiscreteVector() {};

    //virtual void resize(size_t n) = 0;
    virtual size_t size() const = 0;
    virtual double dot(const IDiscreteVector<T> &/*vec*/) {return .0;};
    virtual double norm1() {return .0;};
    virtual double norm2() {return .0;};
    virtual double norm_max() {return .0;};

    virtual T& operator[] (size_t idx) = 0;
    virtual const T& operator[] (size_t idx) const = 0;

    virtual IDiscreteVector<T>& operator= (const IDiscreteVector<T> &src)
    {
        for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] = src[i];
        return *this;
    }
    virtual void operator+= (const IDiscreteVector<T>& v)
    {
        for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] += v[i];
    }
    virtual void operator-= (const IDiscreteVector<T>& v)
    {
        for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] -= v[i];
    }

    virtual IDiscreteVector<T>& operator= (T v)
    {
        for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] = v;
        return *this;
    }

    virtual void addSubvector(const std::vector<size_t> &pos, T* local_v)
    {
        for (size_t i=0; i<pos.size(); ++i) {
            if (pos[i]==BaseLib::index_npos) continue;
            (*this)[pos[i]] += local_v[i];
        }
    }

    //virtual typename std::vector<T>::iterator begin() = 0;
    //virtual typename std::vector<T>::iterator end() = 0;
    virtual size_t getRangeBegin() const = 0;
    virtual size_t getRangeEnd() const = 0;
};

} // end
