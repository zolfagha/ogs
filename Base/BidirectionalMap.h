
#pragma once

#include <map>

namespace Base
{

/**
 * \brief A bidrectional map class between key A and B
 *
 */
template <typename T1, typename T2>
class BidirectionalMap
{
public:
    typedef std::map<T1, T2> MapA;
    typedef std::map<T2, T1> MapB;

    void insert(const T1 &v1, const T2 &v2)
    {
        _map1[v1] = v2;
        _map2[v2] = v1;
    }

    size_t countInA(const T1 &v) const { return _map1.count(v);};

    size_t countInB(const T2 &v) const { return _map2.count(v);};

    size_t size() const { return _map1.size(); };

    void clear() { _map1.clear(); };

    const T2& mapAtoB(const T1 &v) const 
    {
        MapA::const_iterator itr = _map1.find(v);
        //if (itr!=_map1.end()) 
        return itr->second;
    };

    const T1& mapBtoA(const T2 &v) const 
    {
        MapB::const_iterator itr = _map2.find(v);
        //if (itr!=_map2.end()) 
        return itr->second;
    };

private:
    MapA _map1;
    MapB _map2;
};

}
