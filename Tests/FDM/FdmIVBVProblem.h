
#pragma once

#include <vector>

#include "Base/CodingTools.h"
#include "SolutionLib/Problem/AbstractMeshBasedIVBVProblem.h"
#include "FdmFunction.h"
#include "BoundaryConditions.h"

namespace FdmLib
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
	class T_LOCAL_ASSEMBLER_LINEAR
//	class T_LOCAL_ASSEMBLER_RESIDUAL,
//	class T_LOCAL_ASSEMBLER_JACOBIAN
	>
class FdmIVBVProblem : public SolutionLib::AbstractMeshBasedDiscreteIVBVProblem
{
public:
	typedef T_LOCAL_ASSEMBLER_LINEAR LinearAssemblerType;
//	typedef T_LOCAL_ASSEMBLER_RESIDUAL ResidualAssemblerType;
//	typedef T_LOCAL_ASSEMBLER_JACOBIAN JacobianAssemblerType;

	///
	FdmIVBVProblem(	DiscreteLib::DiscreteSystem &dis,
    				MeshLib::IMesh &msh,
    				LinearAssemblerType *linear_assembly
//    				ResidualAssemblerType *residual_assembly,
//    				JacobianAssemblerType *jacobian_assembly
    				)
        : AbstractMeshBasedDiscreteIVBVProblem(msh),
          	_discrete_system(&dis),
			_linear_assembler(linear_assembly)
//			_residual_assembler(residual_assembly),
//			_jacobian_assembler(jacobian_assembly)
    {
        Base::zeroObject(_map_var, _map_ic);
    }

    ///
    virtual ~FdmIVBVProblem()
    {
        Base::releaseObject(_map_var, _map_ic);
        Base::releaseObjectsInStdVector(_map_bc1);
        Base::releaseObjectsInStdVector(_map_bc2);
    }

    /// get the number of variables
    size_t getNumberOfVariables() const { return 1; }

    /// create FD approximation field
    size_t createField()
    {
        if (_map_var==0) {
            _map_var = new FdmLib::FdmFunctionScalar(*_discrete_system, *getMesh());
        }
        return 0;
    }

    /// get the FE field
    FdmLib::FdmFunctionScalar* getField(size_t) const
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
        addDirichletBC(*new FdmLib::FdmDirichletBC<double>(_map_var, &geo, is_transient, &bc1));
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
    FdmLib::FdmDirichletBC<double>* getFdmDirichletBC(int bc_id) const
    {
        return _map_bc1[bc_id];
    };

    ///
    void addNeumannBC(int, GeoLib::GeoObject &geo, bool is_transient, MathLib::SpatialFunctionScalar& bc2)
    {
        addNeumannBC(*new FdmLib::FdmNeumannBC<double, double>(_map_var, &geo, is_transient, &bc2));
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
    FdmLib::FdmNeumannBC<double, double>* getFdmNeumannBC(int bc_id) const
    {
        return _map_bc2[bc_id];
    };

    ///
    LinearAssemblerType* getLinearAssembler() const { return _linear_assembler; }

//    ///
//    ResidualAssemblerType* getResidualAssembler() const { return _residual_assembler; }
//
//    ///
//    JacobianAssemblerType* getJacobianAssembler() const { return _jacobian_assembler; }

private:
    void addDirichletBC(FdmLib::FdmDirichletBC<double>& bc1)
    {
        _map_bc1.push_back(&bc1);
    }

    void addNeumannBC(FdmLib::FdmNeumannBC<double, double>& bc2)
    {
        _map_bc2.push_back(&bc2);
    }

    DISALLOW_COPY_AND_ASSIGN(FdmIVBVProblem);

private:
    DiscreteLib::DiscreteSystem* _discrete_system;
    FdmLib::TemplateFDMFunction<double>* _map_var;
    MathLib::SpatialFunctionScalar* _map_ic;
    std::vector<FdmLib::FdmDirichletBC<double>*> _map_bc1;
    std::vector<FdmLib::FdmNeumannBC<double, double>*> _map_bc2;
    LinearAssemblerType* _linear_assembler;
//    ResidualAssemblerType* _residual_assembler;
//    JacobianAssemblerType* _jacobian_assembler;
};

} //end
