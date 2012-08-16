
#pragma once

#include <vector>
#include <cstddef>

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
class DiscreteVector : public IDiscreteVector<T>
{
protected:
    explicit DiscreteVector(MeshLib::IMesh* msh) : _data(0), _msh(msh) {};
    explicit DiscreteVector(size_t n) : _data(n), _msh(0) {};
    DiscreteVector(MeshLib::IMesh* msh, size_t n) : _data(n), _msh(msh) {};

public:
    static DiscreteVector<T>* createInstance(IDiscreteSystem &sys, size_t n)
    {
        DiscreteVector<T>* vec = new DiscreteVector<T>(n);
        sys.addVector(vec);
        return vec;
    }
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

    /// construct
    virtual void construct(IDiscreteVectorAssembler<T>& assembler)
    {
        assembler.assembly(*_msh, *this);
    }

protected:
    std::vector<T> _data;
    MeshLib::IMesh* _msh;
};

} // end
