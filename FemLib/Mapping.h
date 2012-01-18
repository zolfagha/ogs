
#pragma once

#include "MeshLib/Core/IMesh.h"

namespace FemLib
{

class FemIsoparametricMapping
{
public:
    void configure(MeshLib::IElement* ele);
    double* mapNatural2Physical(double* natural_pt);
    double* mapPhysical2Natural(double* physical_pt);

    void getShapeFunctions(double* natural_pt, double* shapes);
    void getJacobians(double* natural_pt, double* jacob, double &jacob_det);
};

class FemAxisymmetric : public FemIsoparametricMapping
{

};

class FemLowerDimension : public FemIsoparametricMapping
{

};

}
