
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MathLib/Function/Function.h"
#include "GeoLib/Core/Point.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "FemLib/Function/FemFunction.h"
#include "IFemBC.h"


namespace FemLib
{



/**
 * \brief BC1 nodal value constructor?
 *
 */
class IFemDirichletBC
{
public:
    ///
    virtual ~IFemDirichletBC() {};

    /// setup B.C.
    virtual void setup() = 0;

    ///
    virtual std::vector<size_t>& getListOfBCNodes() = 0;
    virtual std::vector<double>& getListOfBCValues() = 0;
};

/**
 * DirichletBC class
 */
class FemDirichletBC : public IFemDirichletBC //: IFemBC, public MathLib::SpatialFunction<Tval>
{
public:
    ///
    FemDirichletBC(MeshLib::IMesh *msh, GeoLib::GeoObject *geo, NumLib::ITXFunctionScalar *bc_func)
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

        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(_msh, _geo, &_vec_nodes);
        // set values
        _vec_values.resize(_vec_nodes.size());
        for (size_t i=0; i<_vec_nodes.size(); i++) {
            const GeoLib::Point* x = _msh->getNodeCoordinatesRef(_vec_nodes[i]);
            _bc_func->eval(x, _vec_values[i]);
        }
        if (!_is_transient)
            _do_setup = false;
    }

    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};
    std::vector<double>& getListOfBCValues() {return _vec_values;};

private:
    MeshLib::IMesh *_msh;
    GeoLib::GeoObject *_geo;
    NumLib::ITXFunctionScalar *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<double> _vec_values;
    bool _is_transient;
    bool _do_setup;
};


}
