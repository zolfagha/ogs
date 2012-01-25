
#pragma  once

#include <map>

template <typename T>
static void destroyStdVectorWithPointers(T &object) {
    if (object.size()>0) {
        const size_t vec_size(object.size());
        for (size_t i=0; i<vec_size; i++)
            delete object[i];
        object.clear();
    }
};

template <typename T>
static void destroyStdMapWithPointers(T &object) {
    if (object.size()>0) {
        for (T::iterator itr=object.begin(); itr!=object.end(); itr++)
            delete itr->second;
        object.clear();
    }
};



