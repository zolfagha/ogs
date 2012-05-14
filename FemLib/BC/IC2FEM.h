

#pragma once

#include <vector>

#include "MeshLib/Core/IMesh.h"
#include "NumLib/Function/TXFunction.h"

namespace FemLib
{

/**
 * \brief A class constructing I.C. data for FEM
 */
class IC2FEM
{
public:
    ///
	IC2FEM(const MeshLib::IMesh &msh, const NumLib::ITXFunction &ic_func, std::vector<size_t> &vec_nodes, std::vector<double> &vec_values);
};


}
