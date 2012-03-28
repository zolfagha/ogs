
#pragma once

#include <cassert>
#include <cmath>
#include <iostream>

namespace MathLib
{

template<typename T, unsigned N>
class TemplateVector
{
public:
    TemplateVector() {};
    TemplateVector(T v)
    {
        for (size_t i = 0; i < N; i++)
            _data[i] = v;
    }
    TemplateVector(T v[N])
    {
        for (size_t i=0; i<N; i++)
            _data[i] = v[i];
    }
    //TemplateVector(T &v1, T &v2);
    //TemplateVector(T &v1, T &v2, T&v3);

    TemplateVector<T,N>& operator= (T a);
    void operator-= (const TemplateVector<T,N> &v);

    const T* getRaw() const {return _data;};
    T* getRawRef() {return _data;};
private:
    T _data[N];
};

template<typename T, unsigned N> TemplateVector<T,N>& TemplateVector<T,N>::operator= (T a)
{
    for (size_t i = 0; i < N; i++)
        _data[i] = a;

    return *this;
}

//template<>
//TemplateVector<double,2>::TemplateVector(double &v1, double &v2)
//{
//    getRawRef()[0] = v1;
//    getRawRef()[1] = v2;
//}
//
//template<>
//TemplateVector<double,3>::TemplateVector(double &v1, double &v2, double &v3)
//{
//    getRawRef()[0] = v1;
//    getRawRef()[1] = v2;
//    getRawRef()[2] = v3;
//}

template<typename T>
class TemplateVector<T,2>
{
public:
    TemplateVector() 
    {
        for (size_t i=0; i<2; i++)
            _data[i] = .0;
    };
    TemplateVector(T v)
    {
        for (size_t i=0; i<2; i++)
            _data[i] = v;
    }
    TemplateVector(T v1, T v2)
    {
        getRawRef()[0] = v1;
        getRawRef()[1] = v2;
    }
    TemplateVector(const TemplateVector&v)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] = v._data[i];
    }

    const T* getRaw() const {return _data;};
    T* getRawRef() {return _data;};

    T& operator[] (size_t idx) {
        assert (idx <= 2);
        return _data[idx];
    }
    const T& operator[] (size_t idx) const {
        assert (idx <= 2);
        return _data[idx];
    }

    TemplateVector<T,2>& operator= (T a)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] = a;

        return *this;
    }
    TemplateVector<T,2>& operator+= (T a)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] += a;

        return *this;
    }
    TemplateVector<T,2>& operator+= (const TemplateVector<T,2> &v)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] += v._data[i];

        return *this;
    }
    TemplateVector<T,2>& operator-= (const TemplateVector<T,2> &v)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] -= v._data[i];

        return *this;
    }
    TemplateVector<T,2>& operator/= (T a)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] /= a;

        return *this;
    }
    TemplateVector<T,2> operator* (T a)
    {
        TemplateVector<T,2> v;
        for (size_t i = 0; i < 2; i++)
            v._data[i] = _data[i]*a;

        return v;
    }
    TemplateVector<T,2> operator/ (T a)
    {
        TemplateVector<T,2> v;
        for (size_t i = 0; i < 2; i++)
            v._data[i] = _data[i]/a;

        return v;
    }

    TemplateVector<T,2> operator- (const TemplateVector<T,2> &ref) const
    {
        TemplateVector<T,2> v;
        for (size_t i = 0; i < 2; i++)
            v[i] = _data[i] - ref._data[i];

        return v;
    }

    bool operator<(const TemplateVector<T,2> &ref) const
    {
    	return this->magnitude() < ref.magnitude();
    }

    T magnitude() const
    {
    	T m = 0;
        for (size_t i = 0; i < 2; i++)
            m += _data[i]*_data[i];
        return std::sqrt(m);
    }

private:
    T _data[2];
};

template<typename T>
class TemplateVector<T,3>
{
public:
    TemplateVector(T v1, T v2, T v3)
    {
        getRawRef()[0] = v1;
        getRawRef()[1] = v2;
        getRawRef()[2] = v3;
    }
    const T* getRaw() const {return _data;};
    T* getRawRef() {return _data;};
private:
    T _data[3];
};

typedef TemplateVector<double, 2> Vector2D;
typedef TemplateVector<double, 3> Vector3D;
}

inline std::ostream& operator<<(std::ostream& output, const MathLib::TemplateVector<double,2>& p)
{
    output << "[";
    for (size_t i=0; i<2; i++)
        output << p[i] << " ";
    output << "]";
    return output;  // for multiple << operators.
}

namespace std
{
inline MathLib::TemplateVector<double,2> abs(const MathLib::TemplateVector<double,2> &v)
{
	MathLib::TemplateVector<double,2> r;
	for (size_t i=0; i<2; i++)
		r[i] = std::abs(v[i]);
	return r;
}

inline double max(double arg0, const MathLib::TemplateVector<double,2> &arg1)
{
	double r = arg0;
	for (size_t i=0; i<2; i++)
		r = std::max(r, arg1[i]);
	return r;
}
}
