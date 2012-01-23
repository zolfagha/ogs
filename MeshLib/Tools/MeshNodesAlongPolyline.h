/*
 * MeshNodesAlongPolyline.h
 *
 *  Created on: Aug 9, 2010
 *      Author: TF
 */

#ifndef MESHNODESALONGPOLYLINE_H_
#define MESHNODESALONGPOLYLINE_H_

// GEOLIB
#include "GeoLib/Polyline.h"

// Mesh
#include "MeshLib/Core/UnstructuredMesh.h"


namespace MeshLib
{
/**
 * This class computes the ids of the mesh nodes along a polyline.
 *
 * The mesh nodes are sorted as follow:
 * [ ... ids of sorted linear nodes ... | ... ids of unsorted higher order nodes ]
 */
class MeshNodesAlongPolyline
{
public:
	MeshNodesAlongPolyline(GeoLib::Polyline const* const ply, UnstructuredMesh2d const* mesh);
	const std::vector<size_t>& getNodeIDs () const;
	const GeoLib::Polyline* getPolyline () const;
	size_t getNumberOfLinearNodes () const;
	std::vector<double> const & getDistOfProjNodeFromPlyStart() const;

private:
	const GeoLib::Polyline* _ply;
	const UnstructuredMesh2d* _mesh;
	size_t _linear_nodes;
	std::vector<size_t> _msh_node_ids;
	std::vector<double> _dist_of_proj_node_from_ply_start;
};
}

#endif /* MESHNODESALONGPOLYLINE_H_ */
