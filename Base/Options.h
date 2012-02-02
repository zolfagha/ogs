
#pragma once

#include <string>
#include <map>

#include "StringTools.h"

namespace Base
{

/**
 * \brief Abstract class of nodes in Options tree
 */
class IOptionNode 
{
public:
    virtual bool isValue() const = 0;
};

/**
 * \brief Leaf in Options tree
 */
class OptionLeaf : public IOptionNode
{
private:
    std::string _value;
public:
    OptionLeaf(const std::string &str) : _value(str) {};
    bool isValue() const {return true;};
    const std::string &getValue() const {return _value;};
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

public:
    Options()
    {
        _parent = 0;
    }

    Options(Options* parent)
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

    /// add new option group 
    Options* addSubGroup(const std::string &key)
    {
        Options *opt = new Options(this);
        _dictionary[key] = opt;
        return opt;
    }

    /// check if there is a value with the given key
    bool hasOption(const std::string &key) const
    {
        Dictionary::const_iterator itr = _dictionary.find(key);
        return (itr!=_dictionary.end() && itr->second->isValue());
    }

    /// get option with the given key if exists.
    const std::string &getOption(const std::string &key) const
    {
        Dictionary::const_iterator itr = _dictionary.find(key);
        if (itr==_dictionary.end() || !itr->second->isValue())
            return "";
        else 
            return static_cast<OptionLeaf*>(itr->second)->getValue();
    }

    /// get value as number
    template<typename T>
    T getOptionAsNum(const std::string &key) const
    {
        return str2number<T>(getOption(key));
    }

    /// add new option
    void addOption(const std::string &key, const std::string &v)
    {
        _dictionary[key] = new OptionLeaf(v);
    }

    /// add option as number
    template<typename T>
    void addOptionAsNum(const std::string &key, T v)
    {
        addOption(key, number2str<T>(v));
    }
};


}
