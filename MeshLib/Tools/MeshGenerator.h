
#pragma once

#include "Mesh.h"

namespace MeshLib
{

class MeshGenerator
{
public:
    static void generateRegularMesh(const int dim, const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z, UnstructuredMesh &mesh) {

        size_t n_eles = static_cast<size_t>(pow(static_cast<double>(subdivision), dim));
        size_t n_nodes = static_cast<size_t>(pow(static_cast<double>(subdivision+1), dim));
        const double unit_length = length / subdivision;
        const size_t n_nodes_per_axis = subdivision+1;

        //nodes
        size_t node_id(0);
        for (size_t i=0; i<n_nodes_per_axis; i++) {
            for (size_t j=0; j<n_nodes_per_axis; j++) {
                for (size_t k=0; k<n_nodes_per_axis; k++) {
                    mesh.setNode(node_id++, unit_length*k, unit_length*j, unit_length*i);
                }
                if (dim<2) break;
            }
            if (dim<3) break;
        }

        //elements
        size_t ele_id(0);
        for (size_t i=0; i<subdivision; i++) {
            for (size_t j=0; j<subdivision; j++) {
                const size_t offset_y1 = j*n_nodes_per_axis;
                const size_t offset_y2 = (j+1)*n_nodes_per_axis;
                for (size_t k=0; k<subdivision; k++) {
                    Quadrirateral *e = new Quadrirateral();
                    e->setNodeID(0, offset_y1+k);
                    e->setNodeID(1, offset_y1+k+1);
                    e->setNodeID(2, offset_y2+k+1);
                    e->setNodeID(3, offset_y2+k);
                    mesh.addElement(e);
                }
                if (dim<2) break;
            }
            if (dim<3) break;
        }

    };

};

}
