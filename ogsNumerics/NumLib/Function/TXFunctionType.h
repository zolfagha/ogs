/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionType.h
 *
 * Created on 2012-11-14 by Norihiro Watanabe
 */

#pragma once

namespace NumLib
{

struct TXFunctionType
{
    enum type
    {
        CONSTANT,
        ANALYTICAL,
        GEOSPACE,
        INVALID
    };
};

inline TXFunctionType::type convertStringToTXFunctionType(const std::string &str)
{
    if (str.compare("CONSTANT")==0) {
        return TXFunctionType::CONSTANT;
    } else if (str.compare("GEOSPACE")==0) {
        return TXFunctionType::GEOSPACE;
    } else if (str.compare("ANALYTICAL")==0) {
        return TXFunctionType::ANALYTICAL;
    } else {
        //error
        return TXFunctionType::INVALID;
    }
}

inline std::string convertTXFunctionTypeToString(const TXFunctionType::type type)
{
    switch (type) {
    case TXFunctionType::CONSTANT:
        return "CONSTANT";
    case TXFunctionType::GEOSPACE:
        return "GEOSPACE";
    case TXFunctionType::ANALYTICAL:
        return "ANALYTICAL";
    default:
        return "";
    }
}

} //end NumLib
