
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
    T* getRaw() {return _data;};
private:
    T _data[N];
};

template<typename T>
class TemplateVector2D : public TemplateVector<T, 2>
{
public:
    TemplateVector2D() {};
    TemplateVector2D(T v1, T v2)
    {
        getRaw()[0] = v1;
        getRaw()[1] = v2;
    }
};

typedef TemplateVector2D<double> Vector2D;
typedef TemplateVector<double, 3> Vector3D;
}
