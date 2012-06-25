

#pragma once

#include <vector>

#include "GeoLib/Core/GeoObject.h"
#include "MeshLib/Core/IMesh.h"
#include "NumLib/Function/TXFunction.h"

namespace FemLib
{

/**
 * \brief A class constructing DirichletBC data for FEM
 */
class DirichletBC2FEM
{
public:
    ///
	DirichletBC2FEM(const MeshLib::IMesh &msh, const GeoLib::GeoObject &geo, const NumLib::ITXFunction &bc_func, std::vector<size_t> &vec_nodes, std::vector<double> &vec_values);
};


}
