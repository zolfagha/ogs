
#pragma  once

#include <cstddef>
#include <map>

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);   \
    TypeName &operator=(const TypeName&)

#define RETURN_ENUM_IF_SAME_STRING(TypeName,str) \
	if (str.compare(#TypeName)==0) return TypeName;

namespace BaseLib
{

const size_t index_npos = -1;

template <typename T>
static void releaseObjectsInStdVector(T &object) 
{
    if (object.size()>0) {
        const size_t vec_size(object.size());
        for (size_t i=0; i<vec_size; i++)
            if (object[i]!=0) delete object[i];
        object.clear();
    }
};

template <typename T>
static void releaseObjectsInStdMap(T &object) 
{
    if (object.size()>0) {
        for (typename T::iterator itr=object.begin(); itr!=object.end(); itr++)
            if (itr->second!=0) delete itr->second;
        object.clear();
    }
};

template <typename T>
static void releaseObjectsInStdQueue(T &container) 
{
    while (!container.empty()) {
        if (container.front() !=0) delete container.front();
        container.pop();
    }
};

template <typename T>
void releaseArrayObject(T* &obj)
{
    if (obj) {
        delete [] obj;
        obj = 0;
    }
}

template <typename T>
void releaseObject(T* &obj)
{
    if (obj) {
        delete obj;
        obj = 0;
    }
}

template <typename T1, typename T2>
void releaseObject(T1* &obj1, T2* &obj2)
{
    releaseObject(obj1);
    releaseObject(obj2);
}

template <typename T1, typename T2, typename T3>
void releaseObject(T1* &obj1, T2* &obj2, T3* &obj3)
{
    releaseObject(obj1, obj2);
    releaseObject(obj3);
}

template <typename T1, typename T2, typename T3, typename T4>
void releaseObject(T1* &obj1, T2* &obj2, T3* &obj3, T4* &obj4)
{
    releaseObject(obj1, obj2);
    releaseObject(obj3, obj4);
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
void releaseObject(T1* &obj1, T2* &obj2, T3* &obj3, T4* &obj4, T5* &obj5)
{
    releaseObject(obj1, obj2);
    releaseObject(obj3, obj4);
    releaseObject(obj5);
}

template <typename T>
void zeroObject(T* &obj)
{
    obj = 0;
}

template <typename T1, typename T2>
void zeroObject(T1* &obj1, T2* &obj2)
{
    zeroObject(obj1);
    zeroObject(obj2);
}

template <typename T1, typename T2, typename T3>
void zeroObject(T1* &obj1, T2* &obj2, T3* &obj3)
{
    zeroObject(obj1, obj2);
    zeroObject(obj3);
}

template <typename T1, typename T2, typename T3, typename T4>
void zeroObject(T1* &obj1, T2* &obj2, T3* &obj3, T4* &obj4)
{
    zeroObject(obj1, obj2, obj3);
    zeroObject(obj4);
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
void zeroObject(T1* &obj1, T2* &obj2, T3* &obj3, T4* &obj4, T5* &obj5)
{
    zeroObject(obj1, obj2, obj3, obj4);
    zeroObject(obj5);
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void zeroObject(T1* &obj1, T2* &obj2, T3* &obj3, T4* &obj4, T5* &obj5, T6* &obj6)
{
    zeroObject(obj1, obj2, obj3, obj4, obj5);
    zeroObject(obj6);
}

}


