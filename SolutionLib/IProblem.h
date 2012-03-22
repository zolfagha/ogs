
#pragma once

#include <vector>
#include <map>

#include "Base/CodingTools.h"
#include "MathLib/Function/Function.h"
#include "MeshLib/Core/IMesh.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"

namespace SolutionLib
{

class IProblem {};


typedef MathLib::IFunction<GeoLib::Point, double> SpatialFunction;

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
    virtual void setIC(int var_type, SpatialFunction& ic) = 0;
    /// get initial condition
    virtual SpatialFunction* getIC(int var_type) const = 0;
    /// add a Dirichlet boundary condition
    virtual void addDirichletBC(int var_type,  GeoLib::GeoObject &geo, bool is_transient, SpatialFunction& bc1) = 0;
    /// get the number of Dirichlet BCs
    virtual size_t getNumberOfDirichletBC(int var_type) const = 0;
    /// get the Dirichlet BC
    virtual SpatialFunction* getDirichletBC(int var_type, int bc_id) const = 0;
    /// add a Neumann boundary condition
    virtual void addNeumannBC(int var_type,  GeoLib::GeoObject &geo, bool is_transient, SpatialFunction& bc2) = 0;
    /// get the number of Neumann BCs
    virtual size_t getNumberOfNeumannBC(int var_type) const = 0;
    /// get the Neumann BC
    virtual SpatialFunction* getNeumannBC(int var_type, int bc_id) const = 0;
};

/**
 * \brief Mesh based discrete IVBV problems
 */
class AbstractMeshBasedDiscreteIVBVProblem : public IVBVProblem
{
public:
    AbstractMeshBasedDiscreteIVBVProblem(MeshLib::IMesh &msh) : _msh(&msh), _tim(0) {};
    virtual ~AbstractMeshBasedDiscreteIVBVProblem()
    {
        Base::releaseObject(_tim);
    }
    /// set a mesh
    void setMesh(MeshLib::IMesh &msh);
    /// get the mesh
    MeshLib::IMesh* getMesh() const {return _msh;};
    /// set a time stepping function
    void setTimeSteppingFunction(NumLib::ITimeStepFunction &f)
    {
        _tim = f.clone();
    }
    /// get the time stepping function
    NumLib::ITimeStepFunction* getTimeSteppingFunction() const {return _tim;};
private:
    MeshLib::IMesh* _msh;
    NumLib::ITimeStepFunction* _tim;
};

}
