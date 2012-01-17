
#pragma once

#include "MeshLib/Core/IMesh.h"

namespace FemLib
{

class IPoint;

class IFunction
{
public:
    virtual double getValue(IPoint &pt) = 0;
};

class FEMInterpolation
{
private:
    MeshLib::IMesh* _msh;
    //shape function
};

class FEMNodalFunction : public IFunction
{
public:
    double getValue(IPoint &pt) {
        return .0;
    };

    double getValue(int node_id) {
        return _nodal_values[node_id];
    }

private:
    double* _nodal_values;
    FEMInterpolation* _fe;
};

class FEMElementalFunction : public IFunction
{

};

}
