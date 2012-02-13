
#pragma  once

#include <map>

namespace Base
{

template <typename T>
static void destroyStdVectorWithPointers(T &object) {
    if (object.size()>0) {
        const size_t vec_size(object.size());
        for (size_t i=0; i<vec_size; i++)
            if (object[i]!=0) delete object[i];
        object.clear();
    }
};

template <typename T>
static void destroyStdMapWithPointers(T &object) {
    if (object.size()>0) {
        for (T::iterator itr=object.begin(); itr!=object.end(); itr++)
            if (itr->second!=0) delete itr->second;
        object.clear();
    }
};

}

