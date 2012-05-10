
#pragma once

#include <vector>

#include "Base/CodingTools.h"
#include "IVBVProblem.h"
#include "MeshBasedProblem.h"
#include "TimeSteppingProblem.h"
#include "LocalAssemblerProblem.h"

#include "FemVariable.h"

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
template
	<
	class T_FEM_EQUATION
	>
class FemIVBVProblem
{
public:
	typedef T_FEM_EQUATION EquationType;

	///
    FemIVBVProblem(	DiscreteLib::DiscreteSystem* dis)
        : _discrete_system(&dis), _eqs(0)
    {
    }

    ///
    virtual ~FemIVBVProblem()
    {
        Base::releaseObjectsInStdVector(_variables);
    }

    /// get this discrete system
    DiscreteLib::DiscreteSystem* getDiscreteSystem() {return _discrete_system;};

    /// create FE approximation field
    FemVariable* addVariable(const std::string name)
    {
    	_variables.push_back(new FemVariable(_variables.size(), name));
        return _variables.back();
    }

    /// get a variable
    FemVariable* getVariable(size_t i) const { return _variables[i]; }

    /// get the number of variables
    size_t getNumberOfVariables() const { return _variables.size(); }

    /// set an equation
    void setEquation(EquationType* eqs) {_eqs = eqs;};

    /// get this equation
    EquationType* getEquation() const {return _eqs;};

private:
    DISALLOW_COPY_AND_ASSIGN(FemIVBVProblem);

private:
    DiscreteLib::DiscreteSystem* _discrete_system;
    EquationType* _eqs;
    std::vector<FemVariable*> _variables;
};

} //end
