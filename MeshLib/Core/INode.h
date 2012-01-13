
#pragma once

namespace MeshLib
{
class INode
{
public:
	INode () {};
	virtual ~INode (){};

	virtual size_t getNodeID(size_t id) = 0;
	virtual void setNodeID(size_t id) = 0;
};

} // end namespace

