

#include "DirichletBC2FEM.h"

#include "GeoLib/Core/Point.h"
#include "MeshLib/Tools/Tools.h"


namespace FemLib
{

DirichletBC2FEM::DirichletBC2FEM(const MeshLib::IMesh &msh, const GeoLib::GeoObject &geo, const NumLib::ITXFunction &bc_func, std::vector<size_t> &vec_nodes, std::vector<double> &vec_values)
{
    // pickup nodes on geometry
    MeshLib::findNodesOnGeometry(&msh, &geo, &vec_nodes);
    const size_t n_bc_nodes = vec_nodes.size();

    if (n_bc_nodes>0) {
        // set values
        vec_values.resize(n_bc_nodes);
        for (size_t i=0; i<n_bc_nodes; i++) {
            const GeoLib::Point* x = msh.getNodeCoordinatesRef(vec_nodes[i]);
            bc_func.eval(x, vec_values[i]);
        }
    } else {
        std::cout << "Dirichlet B.C. was not found." << std::endl;
    }
}

} //end
