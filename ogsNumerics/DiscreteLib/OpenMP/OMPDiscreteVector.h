/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OMPDiscreteVector.h
 *
 * Created on 2012-08-20 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <cstddef>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "BaseLib/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteVectorAssembler.h"

namespace DiscreteLib
{

/**
 * \brief Vector container for single memory
 */
template<typename T>
class OMPDiscreteVector : public IDiscreteVector<T>
{
protected:
    explicit OMPDiscreteVector(MeshLib::IMesh* msh) : _data(0), _msh(msh), _sys(0) {};
    explicit OMPDiscreteVector(size_t n, IDiscreteSystem* sys) : _data(n), _msh(sys->getMesh()), _sys(sys) {};
    OMPDiscreteVector(MeshLib::IMesh* msh, size_t n) : _data(n), _msh(msh), _sys(0) {};

public:
    static OMPDiscreteVector<T>* createInstance(IDiscreteSystem &sys, size_t n)
    {
        OMPDiscreteVector<T>* vec = new OMPDiscreteVector<T>(n, &sys);
        sys.addVector(vec);
        return vec;
    }
    virtual ~OMPDiscreteVector() {};

    virtual OMPDiscreteVector<T>* clone() const
    {
        OMPDiscreteVector<T>* vec = OMPDiscreteVector<T>::createInstance(*_sys, _data.size());
        *vec = (*this);
        return vec;
    }

    virtual void resize(size_t n) {_data.resize(n);};
    virtual size_t size() const {return _data.size();};

    virtual OMPDiscreteVector<T>& operator= (T v)
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

    /// construct
    virtual void construct(IDiscreteVectorAssembler<T>& assembler)
    {
        assembler.assembly(*_msh, *this);
    }

protected:
    std::vector<T> _data;
    MeshLib::IMesh* _msh;
    IDiscreteSystem* _sys;
};

} // end
