
#pragma once

#include "GeoLib/Core/GeoObject.h"
#include "NumLib/TXFunction/TXFunction.h"

#include "IProblem.h"

namespace SolutionLib
{

struct BoundaryConditionType
{
	enum type
	{
		Dirichlet=1,
		Neumann=2,
		Robin=3
	};
};

/**
 * \brief Initial value - boundary value (IVBV) problems
 *
 * IVBV problem has the following elements
 * - variables
 * - equations
 * - IC
 * - BC
 * - domain (space & time)?
 */
class IVBVProblem : IProblem
{
public:
	///
	virtual ~IVBVProblem() {};

	/// get the number of variables
    virtual size_t getNumberOfVariables() const = 0;


    /// set initial condition
    virtual void setIC(size_t var_type, NumLib::ITXFunction &ic) = 0;

    /// get initial condition
    virtual NumLib::ITXFunction* getIC(size_t var_type) const = 0;


    /// add a Dirichlet boundary condition
    virtual void addDirichletBC(size_t var_type,  GeoLib::GeoObject &geo, NumLib::ITXFunction &bc1) = 0;

    /// get the number of Dirichlet BCs
    virtual size_t getNumberOfDirichletBC(size_t var_type) const = 0;

    /// get the Dirichlet BC
    virtual NumLib::ITXFunction* getDirichletBC(size_t var_type, size_t bc_id) const = 0;


    /// add a Neumann boundary condition
    virtual void addNeumannBC(size_t var_type,  GeoLib::GeoObject &geo, NumLib::ITXFunction &bc2) = 0;

    /// get the number of Neumann BCs
    virtual size_t getNumberOfNeumannBC(size_t var_type) const = 0;

    /// get the Neumann BC
    virtual NumLib::ITXFunction* getNeumannBC(size_t var_type, size_t bc_id) const = 0;
};

}
