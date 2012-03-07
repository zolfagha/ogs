
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
    virtual size_t size() const = 0;
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
    DiscreteVector() {};
    DiscreteVector(size_t n) : _data(n) {};
    virtual ~DiscreteVector() {};

    virtual void resize(size_t n) {return _data.resize(n);};
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

} // end
