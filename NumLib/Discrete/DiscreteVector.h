
#pragma once

#include <vector>

namespace NumLib
{

class IDiscreteVector {};

template<typename T>
class DiscreteVector : public IDiscreteVector
{
public:
    DiscreteVector(size_t n) : _data(n) {};

private:
    std::vector<T> _data;
};

template<typename T>
class DecomposedDiscreteVector : public DiscreteVector<T>
{
public:
    DecomposedDiscreteVector(size_t n) : DiscreteVector(n) {};
};


} // end
