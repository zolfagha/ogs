
#pragma once

#include <vector>

#include "GeoLib/Core/Point.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/DirichletBC2FEM.h"

namespace SolutionLib
{

/**
 * \brief A class contains DirichletBC data for FEM
 */
class FemDirichletBC
{
public:
    ///
    FemDirichletBC(MeshLib::IMesh* msh, GeoLib::GeoObject* geo, NumLib::ITXFunction* bc_func)
    {
        _msh = msh;
        _geo = geo;
        _bc_func = bc_func;
        _is_transient = !bc_func->isTemporallyConst();
        _do_setup = true;
    }

    ///
    virtual ~FemDirichletBC()
    {
    }

    /// setup B.C.
    void setup()
    {
        if (!_do_setup) return;
        if (!_is_transient) _do_setup = false;

        FemLib::DirichletBC2FEM convert(*_msh, *_geo, *_bc_func, _vec_nodes, _vec_values);

        if (!_is_transient)
            _do_setup = false;
    }

    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};
    std::vector<double>& getListOfBCValues() {return _vec_values;};

private:
    MeshLib::IMesh* _msh;
    GeoLib::GeoObject* _geo;
    NumLib::ITXFunction* _bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<double> _vec_values;
    bool _is_transient;
    bool _do_setup;
};


}
