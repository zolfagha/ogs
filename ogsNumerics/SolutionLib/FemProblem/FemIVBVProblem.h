
#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"
#include "IVBVProblem.h"
#include "MeshBasedProblem.h"
#include "TimeSteppingProblem.h"
#include "LocalAssemblerProblem.h"

#include "AbstractFemIVBVProblem.h"

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
class FemIVBVProblem : public AbstractFemIVBVProblem
{
public:
    typedef T_FEM_EQUATION EquationType;

    ///
    explicit FemIVBVProblem(    DiscreteLib::DiscreteSystem* dis)
        : AbstractFemIVBVProblem(dis), _eqs(0)
    {
    }

    ///
    virtual ~FemIVBVProblem()
    {
    }

    /// set an equation
    void setEquation(EquationType* eqs) {_eqs = eqs;};

    /// get this equation
    EquationType* getEquation() const {return _eqs;};

private:
    DISALLOW_COPY_AND_ASSIGN(FemIVBVProblem);

private:
    EquationType* _eqs;
};

} //end
