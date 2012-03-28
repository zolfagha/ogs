
#pragma once

#include <vector>


namespace DiscreteLib
{

/**
 * \brief Interface of Vector containers in discrete systems
 */
class IDiscreteVectorBase
{
public:
    virtual ~IDiscreteVectorBase() {};
};

template<typename T>
class IDiscreteVector : public IDiscreteVectorBase
{
public:
    virtual ~IDiscreteVector() {};

    //virtual void resize(size_t n) = 0;
    virtual size_t size() const = 0;
    virtual double dot(const IDiscreteVector<T> &vec) {return .0;};
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
    virtual void operator-= (const IDiscreteVector<T>& v)
    {
    	for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
    		(*this)[i] -= v[i];
    }

    //virtual typename std::vector<T>::iterator begin() = 0;
    //virtual typename std::vector<T>::iterator end() = 0;
    virtual size_t getRangeBegin() const = 0;
    virtual size_t getRangeEnd() const = 0;
};

/**
 * \brief Vector container for single memory
 */
template<typename T>
class DiscreteVector : public IDiscreteVector<T>
{
public:
    DiscreteVector() {};
    DiscreteVector(size_t n) : _data(n) {};
    virtual ~DiscreteVector() {};

    virtual void resize(size_t n) {_data.resize(n);};
    virtual size_t size() const {return _data.size();};
    virtual double dot(const IDiscreteVector<T> &vec) {return .0;};
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
