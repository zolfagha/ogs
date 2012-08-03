/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Iterator.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IMesh.h"

namespace MeshLib
{

//template<class T>
//class ActiveElementIterator : public std::iterator<std::forward_iterator_tag, T>
//{
//private:
//
//public:
//    ActiveElementIterator() {};
//
//    T operator*() const { return m_value; }
//    ActiveElementIterator &operator++() { gen_rand(); ++m_counter; return *this; }
//    bool operator!=(const ActiveElementIterator &x) const { return m_counter != x.m_counter; }
//
//};

}
