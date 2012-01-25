
#pragma once

#include <memory>

#include "MeshLib/Core/UnstructuredMesh.h"
#include "MathLib/Vector.h"
#include "MeshLib/Core/Element.h"

namespace MeshLib
{

class MeshGenerator
{
public:
    static std::auto_ptr<MeshLib::UnstructuredMesh1d> generateLineMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z) 
    {
        std::auto_ptr<MeshLib::UnstructuredMesh1d> msh(new MeshLib::UnstructuredMesh1d());

        size_t n_eles = subdivision;
        size_t n_nodes = subdivision+1;
        const double unit_length = length / subdivision;
        const size_t n_nodes_per_axis = subdivision+1;

        //nodes
        size_t node_id(0);
        for (size_t i_z=0; i_z<n_nodes_per_axis; i_z++) {
            const double x = unit_length*i_z;
            msh->setNode(node_id++, GeoLib::Point(x+origin_x, origin_y, origin_z));
        }

        //elements
        size_t ele_id(0);
        for (size_t i_z=0; i_z<subdivision; i_z++) {
        }

        return msh;
    }

    static std::auto_ptr<MeshLib::UnstructuredMesh2d> generateRegularQuadMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z) 
    {

        std::auto_ptr<MeshLib::UnstructuredMesh2d> msh(new MeshLib::UnstructuredMesh2d());

        size_t n_eles = static_cast<size_t>(pow(static_cast<double>(subdivision), 2));
        size_t n_nodes = static_cast<size_t>(pow(static_cast<double>(subdivision+1), 2));
        const double unit_length = length / subdivision;
        const size_t n_nodes_per_axis = subdivision+1;

        //nodes
        size_t node_id(0);
        const double z = origin_z;
        for (size_t j_y=0; j_y<n_nodes_per_axis; j_y++) {
            const double y = unit_length*j_y + origin_y;
            for (size_t k_x=0; k_x<n_nodes_per_axis; k_x++) {
                const double x = unit_length*k_x + origin_x;
                msh->setNode(node_id++, GeoLib::Point(x, y, z));
            }
        }

        //elements
        size_t ele_id(0);
        for (size_t j=0; j<subdivision; j++) {
            const size_t offset_y1 = j*n_nodes_per_axis;
            const size_t offset_y2 = (j+1)*n_nodes_per_axis;
            for (size_t k=0; k<subdivision; k++) {
                Quadrirateral *e = new Quadrirateral();
                e->setNodeID(0, offset_y1+k);
                e->setNodeID(1, offset_y1+k+1);
                e->setNodeID(2, offset_y2+k+1);
                e->setNodeID(3, offset_y2+k);
                msh->addElement(e);
            }
        }

        return msh;

    };

};

}
