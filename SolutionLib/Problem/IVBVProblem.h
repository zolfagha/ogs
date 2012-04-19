
#pragma once

#include "MathLib/Function/Function.h"
#include "GeoLib/Core/GeoObject.h"

#include "IProblem.h"

namespace SolutionLib
{

/**
 * \brief Initial value boundary value problems
 */
class IVBVProblem : IProblem
{
public:
	virtual ~IVBVProblem() {};

	/// get the number of variables
    virtual size_t getNumberOfVariables() const = 0;

    /// set initial condition
    virtual void setIC(int var_type, MathLib::SpatialFunctionScalar& ic) = 0;

    /// get initial condition
    virtual MathLib::SpatialFunctionScalar* getIC(int var_type) const = 0;

    /// add a Dirichlet boundary condition
    virtual void addDirichletBC(int var_type,  GeoLib::GeoObject &geo, bool is_transient, MathLib::SpatialFunctionScalar& bc1) = 0;

    /// get the number of Dirichlet BCs
    virtual size_t getNumberOfDirichletBC(int var_type) const = 0;

    /// get the Dirichlet BC
    virtual MathLib::SpatialFunctionScalar* getDirichletBC(int var_type, int bc_id) const = 0;

    /// add a Neumann boundary condition
    virtual void addNeumannBC(int var_type,  GeoLib::GeoObject &geo, bool is_transient, MathLib::SpatialFunctionScalar& bc2) = 0;

    /// get the number of Neumann BCs
    virtual size_t getNumberOfNeumannBC(int var_type) const = 0;

    /// get the Neumann BC
    virtual MathLib::SpatialFunctionScalar* getNeumannBC(int var_type, int bc_id) const = 0;
};

}
