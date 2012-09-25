/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Options.h
 *
 * Created on 2012-02-02 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <map>
#include <vector>
#include <ostream>
#include <sstream>
#include <limits>

#include "StringTools.h"

namespace BaseLib
{

/**
 * \brief Abstract class of nodes in Options tree
 */
class OptionNode 
{
public:
    virtual ~OptionNode() {};
    virtual bool isValue() const = 0;
    virtual bool isString() const {return false;};
    virtual std::string getText() const = 0;
    virtual void printout(std::ostream &os, size_t depth=0) const = 0;
};

/**
 * \brief Leaf in Options tree
 */
template <class T>
class OptionLeaf : public OptionNode
{
private:
    T _value;
public:
    OptionLeaf(const T &str) : _value(str) {};
    virtual ~OptionLeaf() {};
    bool isValue() const {return true;};
    virtual inline bool isString() const {return false;};
    const T& getValue() const {return _value;};
    virtual inline std::string getText() const
    {
        std::stringstream ss;
        ss << _value;
        return ss.str();
    };
    virtual inline void printout(std::ostream &/*os*/, size_t /*depth*/) const {};
};

template <>
inline bool OptionLeaf<std::string>::isString() const {return true;};

template <>
inline std::string OptionLeaf<std::string>::getText() const {return _value;};

template <>
inline void OptionLeaf<std::string>::printout(std::ostream &os, size_t depth) const
{
    for (size_t i=0; i<depth; i++)
        os << "\t";
    os << "value: " << _value << std::endl;
};

template <>
inline std::string OptionLeaf<double>::getText() const
{
    std::stringstream ss;
    ss.precision(std::numeric_limits<double>::digits10);
    ss << std::scientific << _value;
    return ss.str();
};

/**
 * \brief Options represents a collection of key-value but with hierarchical data structure.
 */
class Options : public OptionNode
{
private:
    typedef std::multimap<std::string, OptionNode*> DictionaryType;
    typedef std::pair<std::string, OptionNode*> PairType;
    typedef std::pair<DictionaryType::iterator, DictionaryType::iterator> pair_of_iterator;
    typedef std::pair<DictionaryType::const_iterator, DictionaryType::const_iterator> const_pair_of_iterator;
    DictionaryType _dictionary;
    Options* _parent;
    const std::string _dummy;
    mutable const_pair_of_iterator _subgroup_range;
    mutable DictionaryType::const_iterator _subgroup_itr;
    mutable const_pair_of_iterator _leaf_range;
    mutable DictionaryType::const_iterator _leaf_itr;

    // return dummy
    template<typename T>
    inline T getDummy() const
    {
        return 0;
    };

public:
    typedef DictionaryType::iterator iterator;
    typedef DictionaryType::const_iterator const_iterator;

    ///
    Options() : _dummy("")
    {
        _parent = 0;
    }

    ///
    explicit Options(Options* parent) : _dummy("")
    {
        _parent = parent;
    }

    ///
    virtual ~Options()
    {
        for (DictionaryType::iterator itr=_dictionary.begin(); itr!=_dictionary.end(); itr++)
            delete itr->second;
    }

    /// check if this node is leaf or not
    bool isValue() const {return false;};

    /// get parent option group
    Options* getParent() {return _parent;};

    ///
    const_iterator begin() const { return _dictionary.begin(); };

    ///
    const_iterator end() const { return _dictionary.end(); };

    /// add new option group 
    Options* addSubGroup(const std::string &key)
    {
        Options *opt = new Options(this);
        _dictionary.insert(PairType(key, opt));
        return opt;
    }

    /// return if a subgroup with the name exists
    bool hasSubGroup(const std::string &key) const
    {
        DictionaryType::const_iterator itr = _dictionary.find(key);
        return (itr!=_dictionary.end() && !itr->second->isValue());
    }

    /// get option with the given key if it exists.
    Options* getSubGroup(const std::string &key)
    {
       DictionaryType::const_iterator itr = _dictionary.find(key);
       if (itr==_dictionary.end() || itr->second->isValue())
           return 0;
       else
           return static_cast<Options*>(itr->second);
    }

    /// get option with the given key if it exists.
    const Options* getSubGroup(const std::string &key) const
    {
       DictionaryType::const_iterator itr = _dictionary.find(key);
       if (itr==_dictionary.end() || itr->second->isValue())
           return 0;
       else
           return static_cast<Options*>(itr->second);
    }

    /// get option with the given key if it exists.
    const Options* getFirstSubGroup(const std::string &key) const
    {
        _subgroup_range = _dictionary.equal_range(key);
        _subgroup_itr = _subgroup_range.first;
        if (_subgroup_itr == _subgroup_range.second || _subgroup_itr->second->isValue())
           return 0;
        else
           return static_cast<Options*>(_subgroup_itr->second);
    }

    /// get option with the given key if it exists.
    const Options* getNextSubGroup() const
    {
        if (_subgroup_itr == _subgroup_range.second)
            return 0;
        ++_subgroup_itr;

        if (_subgroup_itr == _subgroup_range.second || _subgroup_itr->second->isValue())
           return 0;
        else
           return static_cast<Options*>(_subgroup_itr->second);
    }

    /// check if there is a value with the given key
    bool hasOption(const std::string &key) const
    {
        DictionaryType::const_iterator itr = _dictionary.find(key);
        return (itr!=_dictionary.end() && itr->second->isValue());
    }

    /// get value as string
    const std::string getOption(const std::string &key) const
    {
        DictionaryType::const_iterator itr = _dictionary.find(key);
        if (itr==_dictionary.end() || !itr->second->isValue()) {
            return _dummy;
        } else {
            return itr->second->getText();
        }
    }

    /// get raw data
    template<typename T>
    const T getOptionAsNum(const std::string &key) const
    {
        DictionaryType::const_iterator itr = _dictionary.find(key);
        if (itr==_dictionary.end() || !itr->second->isValue()) {
            return getDummy<T>();
        } else {
           if (itr->second->isString()) {
                return BaseLib::str2number<T>(static_cast<OptionLeaf<std::string>*>(itr->second)->getValue());
            } else {
                return static_cast<OptionLeaf<T>*>(itr->second)->getValue();
            }
        }
        //        return Base::str2number<T>(getOption(key));
    }

//    /// get the option value
//    template<typename T>
//    const std::vector<T>* getOptionAsArray(const std::string &key) const
//    {
//        DictionaryType::const_iterator itr = _dictionary.find(key);
//        if (itr==_dictionary.end() || !itr->second->isValue())
//            return 0;
//        else
//            return &static_cast<OptionLeaf<std::vector<T> >*>(itr->second)->getValue();
//    }

    template<typename T>
    const T getFirstOption(const std::string &key) const
    {
        _leaf_range = _dictionary.equal_range(key);
        _leaf_itr = _leaf_range.first;
        if (_leaf_itr == _leaf_range.second || !_leaf_itr->second->isValue())
            return getDummy<T>();
        else
            return static_cast<OptionLeaf<T>*>(_leaf_itr->second)->getValue();
    }

    template<typename T>
    const T getNextOption() const
    {
        if (_leaf_itr == _leaf_range.second)
            return getDummy<T>();
        ++_leaf_itr;

        if (_leaf_itr == _leaf_range.second || !_leaf_itr->second->isValue())
            return getDummy<T>();
        else
            return static_cast<OptionLeaf<T>*>(_leaf_itr->second)->getValue();
    }

//    /// get the option value
//    const std::string getOption(const std::string &key) const
//    {
//        return getOption<std::string>(key);
//    }

    /// add new option
    void addOption(const std::string &key, const std::string &v)
    {
        _dictionary.insert(PairType(key, new OptionLeaf<std::string>(v)));
    }

    /// add option as number
    template<typename T>
    void addOptionAsNum(const std::string &key, T v)
    {
        _dictionary.insert(PairType(key, new OptionLeaf<T>(v)));
    }

//    template <class T>
//    void addOptionAsArray(const std::string &key, const std::vector<T> &v)
//    {
//        _dictionary.insert(PairType(key, new OptionLeaf<std::vector<T> >(v)));
//    }

    virtual inline std::string getText() const {return "";};

    virtual void printout (std::ostream &os, size_t depth=0) const
    {
        DictionaryType::const_iterator itr = _dictionary.begin();
        for (; itr!=_dictionary.end(); ++itr) {
            for (size_t i=0; i<depth; i++)
                os << "\t";
            os << "tag: " << itr->first << std::endl;
            //os << ", is_value: " << (itr->second->isValue() ? "yes" : "no") << std::endl;
            itr->second->printout(os, depth+1);
        }

    }
};

//template<>
//inline std::string Options::getDummy<std::string>() const
//{
//    return _dummy;
//};

}
