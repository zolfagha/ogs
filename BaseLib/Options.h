
#pragma once

#include <string>
#include <map>

#include "StringTools.h"

namespace BaseLib
{

/**
 * \brief Abstract class of nodes in Options tree
 */
class IOptionNode 
{
public:
	virtual ~IOptionNode() {};
    virtual bool isValue() const = 0;
};

/**
 * \brief Leaf in Options tree
 */
template <class T>
class OptionLeaf : public IOptionNode
{
private:
    T _value;
//    std::string _value;
public:
    OptionLeaf(const T &str) : _value(str) {};
    bool isValue() const {return true;};
    const T& getValue() const {return _value;};
};

/**
 * \brief Options represents a collection of key-value but with hierarchical data structure.
 */
class Options : public IOptionNode
{
private:
    typedef std::map<std::string, IOptionNode*> Dictionary;
    Dictionary _dictionary;
    Options* _parent;
    const std::string _dummy;

public:
    typedef Dictionary::const_iterator const_iterator;

    Options() : _dummy("")
    {
        _parent = 0;
    }

    Options(Options* parent) : _dummy("")
    {
        _parent = parent;
    }

    ~Options()
    {
        for (Dictionary::iterator itr=_dictionary.begin(); itr!=_dictionary.end(); itr++)
            delete itr->second;
    }

    /// check if this node is leaf or not
    bool isValue() const {return false;};

    /// get parent option group
    Options* getParent() {return _parent;};

    /// get option with the given key if it exists.
    const Options* getSubGroup(const std::string &key) const
    {
       Dictionary::const_iterator itr = _dictionary.find(key);
       if (itr==_dictionary.end() || itr->second->isValue())
           return 0;
       else
           return static_cast<Options*>(itr->second);
    }

    const_iterator begin() const { return _dictionary.begin(); };
    const_iterator end() const { return _dictionary.end(); };

    /// add new option group 
    Options* addSubGroup(const std::string &key)
    {
        Options *opt = new Options(this);
        _dictionary[key] = opt;
        return opt;
    }

    bool hasSubGroup(const std::string &key) const
    {
        Dictionary::const_iterator itr = _dictionary.find(key);
        return (itr!=_dictionary.end() && !itr->second->isValue());
    }

    /// check if there is a value with the given key
    bool hasOption(const std::string &key) const
    {
        Dictionary::const_iterator itr = _dictionary.find(key);
        return (itr!=_dictionary.end() && itr->second->isValue());
    }

//    /// get option with the given key if exists.
//    const std::string &getOption(const std::string &key) const
//    {
//        Dictionary::const_iterator itr = _dictionary.find(key);
//        if (itr==_dictionary.end() || !itr->second->isValue())
//            return _dummy;
//        else
//            return static_cast<OptionLeaf*>(itr->second)->getValue();
//    }

    template<typename T>
    inline T getDummy() const
    {
    	return 0;
    };

    /// get value as number
    template<typename T>
    const T getOption(const std::string &key) const
    {
        Dictionary::const_iterator itr = _dictionary.find(key);
        if (itr==_dictionary.end() || !itr->second->isValue())
            return getDummy<T>();
        else
            return static_cast<OptionLeaf<T>*>(itr->second)->getValue();
//        return Base::str2number<T>(getOption(key));
    }

    template<typename T>
    const std::vector<T>* getOptionAsArray(const std::string &key) const
    {
        Dictionary::const_iterator itr = _dictionary.find(key);
        if (itr==_dictionary.end() || !itr->second->isValue())
            return 0;
        else
            return &static_cast<OptionLeaf<std::vector<T> >*>(itr->second)->getValue();
    }

    const std::string getOption(const std::string &key) const
    {
    	return getOption<std::string>(key);
    }

    /// add new option
    void addOption(const std::string &key, const std::string &v)
    {
        _dictionary[key] = new OptionLeaf<std::string>(v);
    }

    /// add option as number
    template<typename T>
    void addOptionAsNum(const std::string &key, T v)
    {
        _dictionary[key] = new OptionLeaf<T>(v);
    }

    template <class T>
    void addOptionAsArray(const std::string &key, const std::vector<T> &v)
    {
        _dictionary[key] = new OptionLeaf<std::vector<T> >(v);
    }


};

template<>
inline std::string Options::getDummy<std::string>() const
{
	return _dummy;
};



}
