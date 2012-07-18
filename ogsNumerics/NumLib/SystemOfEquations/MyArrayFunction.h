
#pragma once

#include "NumLib/Function/IFunction.h"

namespace NumLib
{

template <typename Tarray>
class MyArrayFunction : public IFunction
{
public:
    MyArrayFunction(const Tarray &array)
    {
        _array = new Tarray(array);
    };
    virtual ~MyArrayFunction()
    {
        delete _array;
    };

    Tarray* getArray() { return _array;};
    Tarray* getArray() const { return _array;};

    virtual MyArrayFunction* clone() const
    {
        return new MyArrayFunction(*_array);
    }

private:
    Tarray* _array;
};

}
