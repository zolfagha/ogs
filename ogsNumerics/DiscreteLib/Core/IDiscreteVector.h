/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IDiscreteVector.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <cstddef>

#include "BaseLib/CodingTools.h"
#include "IDiscreteObject.h"

namespace DiscreteLib
{

/**
 * \brief Interface of vector class
 * 
 */
template<typename T>
class IDiscreteVector : public IDiscreteObject
{
public:
    typedef IDiscreteVector<T> MyVectorType;
    
    virtual ~IDiscreteVector() {};

    /// return vector size
    virtual size_t size() const = 0;
    /// return a start index of the active data range 
    virtual size_t getRangeBegin() const = 0;
    /// return an end index of the active data range 
    virtual size_t getRangeEnd() const = 0;

    /// access data
    virtual T& operator[] (size_t idx) = 0;
    virtual const T& operator[] (size_t idx) const = 0;
    /// vector operation: set data
    virtual MyVectorType& operator= (const MyVectorType &src)
    {
        for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] = src[i];
        return *this;
    }
    /// vector operation: add 
    virtual void operator+= (const MyVectorType& v)
    {
        for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] += v[i];
    }
    /// vector operation: subtract
    virtual void operator-= (const MyVectorType& v)
    {
        for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] -= v[i];
    }
    /// set all values in this vector
    virtual MyVectorType& operator= (T v)
    {
        for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] = v;
        return *this;
    }
    /// add values to given entries
    virtual void addSubvector(const std::vector<size_t> &pos, T* local_v)
    {
        for (size_t i=0; i<pos.size(); ++i) {
            if (pos[i]==BaseLib::index_npos) continue;
            (*this)[pos[i]] += local_v[i];
        }
    }

};

} // end
