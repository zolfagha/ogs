
#pragma once

#include <vector>
#include "Base/MemoryTools.h"
#include "GeoLib/Point.h"
#include "MeshLib/Core/IElement.h"

namespace MeshLib
{

class IElementMapping
{
public:
    virtual GeoLib::Point* getNodePoint(size_t node_id) = 0;
};

class EleMapInvariant : public IElementMapping
{
public:
    EleMapInvariant(IElement* e) 
    {
        _e = e;
    };
    virtual GeoLib::Point* getNodePoint(size_t node_id) 
    {
        return (GeoLib::Point*)_e->getNodeLocalCoordinates(node_id);
    }
private:
    IElement *_e;
};

class EleMapLocalCoordinates : public IElementMapping 
{
public:
    EleMapLocalCoordinates(IElement* e) 
    {
        //move

        //rotate

    };
    virtual ~EleMapLocalCoordinates()
    {
        destroyStdVectorWithPointers(_point_vec);
    }

    virtual GeoLib::Point* getNodePoint(size_t node_id) {
        return _point_vec[node_id];
    }
private:
    std::vector<GeoLib::Point*> _point_vec;
};

}
