
#pragma once

#include <vector>

#include "Base/CodingTools.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"

#include "IVBVProblem.h"
#include "MeshBasedProblem.h"
#include "TimeSteppingProblem.h"
#include "LocalAssemblerProblem.h"

namespace SolutionLib
{


class FemVariable
{
public:
	FemVariable(size_t id, const std::string &name) : _id(id), _name(name)
	{

	}

	//----------------------------------------------------------------------
	size_t getID() const {return _id;};
	const std::string& getName() const { return _name;}


	//----------------------------------------------------------------------
    void setIC(NumLib::ITXFunction* ic) { _f_ic = ic; };
    NumLib::ITXFunction* getIC() const { return _f_ic; };


	//----------------------------------------------------------------------
    void addDirichletBC(FemLib::IFemDirichletBC *bc)
    {
        _map_bc1.push_back(bc);
    }
    size_t getNumberOfDirichletBC() const {return _map_bc1.size();};
    FemLib::IFemDirichletBC* getDirichletBC(size_t bc_id) const
    {
        return _map_bc1[bc_id];
    };


	//----------------------------------------------------------------------
    void addNeumannBC(FemLib::IFemNeumannBC& bc2)
    {
        _map_bc2.push_back(&bc2);
    }
    size_t getNumberOfNeumannBC() const {return _map_bc2.size();};
    FemLib::IFemNeumannBC* getNeumannBC(size_t bc_id) const
    {
        return _map_bc2[bc_id];
    };

private:
    size_t _id;
    std::string _name;
    NumLib::ITXFunction* _f_ic;
    std::vector<FemLib::IFemDirichletBC*> _map_bc1;
    std::vector<FemLib::IFemNeumannBC*> _map_bc2;
};


/**
 * \brief IVBV problems using FEM
 *
 *- Variables
 *- Equations (local assembly)
 *- IC
 *- BC
 */
template
	<
	class T_LOCAL_ASSEMBLER_LINEAR,
	class T_LOCAL_ASSEMBLER_RESIDUAL,
	class T_LOCAL_ASSEMBLER_JACOBIAN
	>
class FemIVBVProblem
: //public IVBVProblem,
  public MeshBasedProblem,
  public TimeSteppingProblem,
  public LocalAssemblerProblem<T_LOCAL_ASSEMBLER_LINEAR, T_LOCAL_ASSEMBLER_RESIDUAL, T_LOCAL_ASSEMBLER_JACOBIAN>
{
public:
	typedef T_LOCAL_ASSEMBLER_LINEAR LinearAssemblerType;
	typedef T_LOCAL_ASSEMBLER_RESIDUAL ResidualAssemblerType;
	typedef T_LOCAL_ASSEMBLER_JACOBIAN JacobianAssemblerType;
	typedef LocalAssemblerProblem<T_LOCAL_ASSEMBLER_LINEAR, T_LOCAL_ASSEMBLER_RESIDUAL, T_LOCAL_ASSEMBLER_JACOBIAN> MyLocalAssemblerProblem;

	///
    FemIVBVProblem(	DiscreteLib::DiscreteSystem &dis,
    				LinearAssemblerType *linear_assembly,
    				ResidualAssemblerType *residual_assembly,
    				JacobianAssemblerType *jacobian_assembly
    				)
        : MeshBasedProblem(dis.getMesh()),
        	MyLocalAssemblerProblem(linear_assembly, residual_assembly, jacobian_assembly),
          	_discrete_system(&dis)
    {
    }


    ///
    virtual ~FemIVBVProblem()
    {
        Base::releaseObjectsInStdVector(_variables);
    }

    /// get the number of variables
    size_t getNumberOfVariables() const { return _variables.size(); }

    /// create FE approximation field
    FemVariable* createVariables(const std::string name)
    {
    	_variables.push_back(new FemVariable(_variables.size(), name));
        return _variables.back();
    }

    /// get the FE field
    FemVariable* getVariable(size_t i) const { return _variables[i]; }

private:
    DISALLOW_COPY_AND_ASSIGN(FemIVBVProblem);

private:
    DiscreteLib::DiscreteSystem* _discrete_system;
    std::vector<FemVariable*> _variables;
};

} //end
