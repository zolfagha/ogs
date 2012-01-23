
#pragma once

namespace MathLib
{
template<typename T, unsigned N>
class TemplateVector
{
public:
    TemplateVector() {};
    TemplateVector(T v[N])
    {
        for (size_t i=0; i<N; i++)
            _data[u] = v[i];
    }

    const T* getRaw() const {return _data;};
    T* getRawRef() {return _data;};
private:
    T _data[N];
};

template<typename T>
class TemplateVector<T,2>
{
public:
    TemplateVector() {};
    TemplateVector(T v1, T v2)
    {
        getRawRef()[0] = v1;
        getRawRef()[1] = v2;
    }
    const T* getRaw() const {return _data;};
    T* getRawRef() {return _data;};
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
