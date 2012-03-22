
#pragma once

#include "Base/CodingTools.h"
#include "Base/BidirectionalMap.h"

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "DiscreteVector.h"
#include "DiscreteLinearEquation.h"
#include "DiscreteLib/Utils/DiscreteDataContainer.h"
//#include "MeshBasedDiscreteLinearEquation.h"

namespace DiscreteLib
{

class DofEquationIdTable;

/**
 *\brief Interface class for all kinds of discrete systems
 */
class IDiscreteSystem
{
public:
    virtual ~IDiscreteSystem() {};
};

/**
 * \brief Discrete system based on a single mesh
 */
class DiscreteSystem : public IDiscreteSystem
{
public:
    /// @param msh a mesh to represent discrete systems by nodes or elements or edges
    DiscreteSystem(MeshLib::IMesh& msh) : _msh(&msh) {};

    ///
    virtual ~DiscreteSystem() {};

    /// get this mesh
    MeshLib::IMesh* getMesh() const { return _msh; };

    /// create a new linear equation
    /// @tparam T_LINEAR_EQUATION
    /// @tparam T_LINEAR_SOLVER
    /// @tparam T_SPARSITY_BUILDER
    /// @param linear_solver 		Linear solver
    /// @param dofManager			Equation index table
	template <template <class, class> class T_LINEAR_EQUATION, class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
	IDiscreteLinearEquation* createLinearEquation(T_LINEAR_SOLVER &linear_solver, DofEquationIdTable &dofManager)
	{
		IDiscreteLinearEquation *eq = new T_LINEAR_EQUATION<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_msh, linear_solver, dofManager);
		_data.addLinearEquation(eq);
		return eq;
	}

	/// create a new vector
	/// @param n	vector length
	/// @return vector object
    template<typename T>
    DiscreteVector<T>* createVector(const size_t n)
    {
        DiscreteVector<T>* v = new DiscreteVector<T>(n);
        _data.addVector(v);
        return v;
    };

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteSystem);

protected:
    MeshLib::IMesh* _msh;
    DiscreteDataContainer _data;
};

} //end
