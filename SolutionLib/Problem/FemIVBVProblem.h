
#pragma once

#include <vector>

#include "Base/CodingTools.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"

#include "AbstractMeshBasedIVBVProblem.h"

namespace SolutionLib
{

/**
 * \brief IVBV problems using FEM
 *
 *- Variables
 *- Equations (local assembly)
 *- IC
 *- BC
 */
template <
	class T_ASSEMBLER_RESIDUAL,
	class T_ASSEMBLER_JACOBIAN
	>
class FemIVBVProblem : public AbstractMeshBasedDiscreteIVBVProblem
{
public:
	typedef T_ASSEMBLER_RESIDUAL ReisdualAssemblerType;
	typedef T_ASSEMBLER_JACOBIAN JacobianAssemblerType;

	///
    FemIVBVProblem(	DiscreteLib::DiscreteSystem &dis,
    				MeshLib::IMesh &msh,
    				ReisdualAssemblerType &user_assembly
    				)
        : AbstractMeshBasedDiscreteIVBVProblem(msh),
          _residual_assembler(user_assembly), _jacobian_assembler(0),
          _discrete_system(&dis)
    {
        Base::zeroObject(_map_var, _map_ic);
    }

    FemIVBVProblem(	DiscreteLib::DiscreteSystem &dis,
    				MeshLib::IMesh &msh,
    				ReisdualAssemblerType &user_assembly,
    				JacobianAssemblerType *jacobian_assembler
    				)
        : AbstractMeshBasedDiscreteIVBVProblem(msh),
          _residual_assembler(user_assembly), _jacobian_assembler(jacobian_assembler),
          _discrete_system(&dis)
    {
        Base::zeroObject(_map_var, _map_ic);
    }

    ///
    virtual ~FemIVBVProblem()
    {
        Base::releaseObject(_map_var, _map_ic);
        Base::releaseObjectsInStdVector(_map_bc1);
        Base::releaseObjectsInStdVector(_map_bc2);
    }

    /// get the number of variables
    size_t getNumberOfVariables() const { return 1; }

    /// create FE approximation field
    size_t createField(FemLib::PolynomialOrder::type order)
    {
        if (_map_var==0) {
            _map_var = new FemLib::FemNodalFunctionScalar(*_discrete_system, *getMesh(), order);
        }
        return 0;
    }

    /// get the FE field
    FemLib::FemNodalFunctionScalar* getField(size_t) const
    {
        return _map_var;
    }

    ///
    void setIC(int, MathLib::SpatialFunctionScalar& ic)
    {
        _map_ic = (MathLib::SpatialFunctionScalar*) ic.clone();
    }

    ///
    MathLib::SpatialFunctionScalar* getIC(int) const
    {
        return _map_ic;
    };

    ///
    void addDirichletBC(int, GeoLib::GeoObject &geo, bool is_transient, MathLib::SpatialFunctionScalar& bc1)
    {
        addDirichletBC(*new FemLib::FemDirichletBC<double>(_map_var, &geo, is_transient, &bc1, new FemLib::DiagonalizeMethod()));
    }


    ///
    size_t getNumberOfDirichletBC(int) const {return getNumberOfDirichletBC();};
    size_t getNumberOfDirichletBC() const {return _map_bc1.size();};

    ///
    MathLib::SpatialFunctionScalar* getDirichletBC(int, int bc_id) const
    {
        return _map_bc1[bc_id];
    };

    ///
    FemLib::FemDirichletBC<double>* getFemDirichletBC(int bc_id) const 
    {
        return _map_bc1[bc_id];
    };

    ///
    void addNeumannBC(int, GeoLib::GeoObject &geo, bool is_transient, MathLib::SpatialFunctionScalar& bc2)
    {
        addNeumannBC(*new FemLib::FemNeumannBC<double, double>(_map_var, &geo, is_transient, &bc2));
    }

    ///
    size_t getNumberOfNeumannBC(int) const {return getNumberOfNeumannBC();};
    size_t getNumberOfNeumannBC() const {return _map_bc2.size();};

    ///
    MathLib::SpatialFunctionScalar* getNeumannBC(int, int bc_id) const
    {
        return _map_bc2[bc_id];
    };

    ///
    FemLib::FemNeumannBC<double, double>* getFemNeumannBC(int bc_id) const 
    {
        return _map_bc2[bc_id];
    };

    ReisdualAssemblerType& getResidualAssembler()
    {
        return _residual_assembler;
    }

    JacobianAssemblerType* getJacobianAssembler()
    {
        return _jacobian_assembler;
    }

private:
    void addDirichletBC(FemLib::FemDirichletBC<double>& bc1)
    {
        _map_bc1.push_back(&bc1);
    }

    void addNeumannBC(FemLib::FemNeumannBC<double, double>& bc2)
    {
        _map_bc2.push_back(&bc2);
    }

    DISALLOW_COPY_AND_ASSIGN(FemIVBVProblem);

private:
    ReisdualAssemblerType _residual_assembler;
    JacobianAssemblerType *_jacobian_assembler;
    DiscreteLib::DiscreteSystem* _discrete_system;
    FemLib::FemNodalFunctionScalar* _map_var;
    MathLib::SpatialFunctionScalar* _map_ic;
    std::vector<FemLib::FemDirichletBC<double>*> _map_bc1;
    std::vector<FemLib::FemNeumannBC<double, double>*> _map_bc2;
};


} //end
