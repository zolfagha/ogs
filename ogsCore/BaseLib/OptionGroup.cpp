/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OptionGroup.cpp
 *
 * Created on 2012-07-25 by Norihiro Watanabe
 */

#include "OptionGroup.h"

namespace BaseLib
{

///
OptionGroup::~OptionGroup()
{
    for (DictionaryType::iterator itr=_dictionary.begin(); itr!=_dictionary.end(); itr++)
        delete itr->second;
}


/// add new option group
OptionGroup* OptionGroup::addSubGroup(const std::string &key)
{
    OptionGroup *opt = new OptionGroup(this);
    _dictionary.insert(PairType(key, opt));
    return opt;
}

/// return if a subgroup with the name exists
bool OptionGroup::hasSubGroup(const std::string &key) const
{
    DictionaryType::const_iterator itr = _dictionary.find(key);
    return (itr!=_dictionary.end() && !itr->second->isValue());
}

/// get option with the given key if it exists.
OptionGroup* OptionGroup::getSubGroup(const std::string &key)
{
   DictionaryType::const_iterator itr = _dictionary.find(key);
   if (itr==_dictionary.end() || itr->second->isValue())
       return NULL;
   else
       return static_cast<OptionGroup*>(itr->second);
}

/// get option with the given key if it exists.
const OptionGroup* OptionGroup::getSubGroup(const std::string &key) const
{
   DictionaryType::const_iterator itr = _dictionary.find(key);
   if (itr==_dictionary.end() || itr->second->isValue())
       return 0;
   else
       return static_cast<OptionGroup*>(itr->second);
}

/// get option with the given key if it exists.
const OptionGroup* OptionGroup::getFirstSubGroup(const std::string &key) const
{
    _subgroup_range = _dictionary.equal_range(key);
    _subgroup_itr = _subgroup_range.first;
    if (_subgroup_itr == _subgroup_range.second || _subgroup_itr->second->isValue())
       return 0;
    else
       return static_cast<OptionGroup*>(_subgroup_itr->second);
}

/// get option with the given key if it exists.
const OptionGroup* OptionGroup::getNextSubGroup() const
{
    if (_subgroup_itr == _subgroup_range.second)
        return 0;
    ++_subgroup_itr;

    if (_subgroup_itr == _subgroup_range.second || _subgroup_itr->second->isValue())
       return 0;
    else
       return static_cast<OptionGroup*>(_subgroup_itr->second);
}

/// check if there is a value with the given key
bool OptionGroup::hasOption(const std::string &key) const
{
    DictionaryType::const_iterator itr = _dictionary.find(key);
    return (itr!=_dictionary.end() && itr->second->isValue());
}

/// get value as string
const std::string OptionGroup::getOption(const std::string &key) const
{
    DictionaryType::const_iterator itr = _dictionary.find(key);
    if (itr==_dictionary.end() || !itr->second->isValue()) {
        return getDummy<std::string>();
    } else {
        return itr->second->getText();
    }
}

/// add new option
void OptionGroup::addOption(const std::string &key, const std::string &v)
{
    _dictionary.insert(PairType(key, new OptionLeaf<std::string>(v)));
}

void OptionGroup::printout (std::ostream &os, size_t depth) const
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

}
