
#pragma once

#include <vector>
#include "MeshLib/Core/IElement.h"

namespace NumLib
{

class DoF
{
};

class DoFManager
{
public:
    void getMap( const MeshLib::IElement * e, std::vector<size_t>& vec ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }


};
}
