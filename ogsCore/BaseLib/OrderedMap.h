
#pragma once

#include <map>
#include <vector>

namespace BaseLib
{

template <typename T_KEY, typename T_VAL>
class OrderedMap
{
public:
    typedef typename std::map<T_KEY, T_VAL>::iterator iterator;
    typedef typename std::map<T_KEY, T_VAL>::const_iterator const_iterator;

    void insert(const T_KEY &key, const T_VAL &val)
    {
        if (_map.count(key)==0)
            _vec.push_back(key);
        _map[key] = val;
    }

    size_t count(const T_KEY &key) const
    {
        return _map.count(key);
    }

    size_t size() const {return _vec.size();};

    iterator begin() {return _map.begin();};
    const_iterator begin() const {return _map.begin();};
    iterator end() {return _map.end();};
    const_iterator end() const {return _map.end();};

    void clear()
    {
        _map.clear();
        _vec.clear();
    }

    iterator operator[] (size_t idx)
    {
        return _map.find(_vec[idx]);
    }

private:
    std::map<T_KEY, T_VAL> _map;
    std::vector<T_KEY> _vec;
};

}
