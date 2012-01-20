
#pragma once

namespace MeshLib
{
template<typename Tpos> 
class INode
{
public:
	INode () {};
	virtual ~INode (){};

	virtual size_t getNodeID(size_t id) = 0;
	virtual void setNodeID(size_t id) = 0;
    virtual const Tpos* getData() const = 0;
};

} // end namespace

