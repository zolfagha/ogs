/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTransientResidualAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <valarray>

#include "DiscreteLib/EquationId/DofEquationIdTable.h"
#include "DiscreteLib/Vector/DiscreteVector.h"
#include "DiscreteLib/Assembler/IDiscreteVectorAssembler.h"

#include "IElementWiseTransientResidualLocalAssembler.h"

namespace MeshLib
{
class IMesh;
}

namespace NumLib
{

class TimeStep;


/**
 * \brief Element-based discrete system assembler classes
 */
class ElementWiseTransientResidualAssembler : public DiscreteLib::IDiscreteVectorAssembler<double>
{
public:
    typedef DiscreteLib::IDiscreteVectorAssembler<double>::VectorType GlobalVectorType;

    /// @param time
    /// @param u0
    /// @param u1
    /// @param a
    ElementWiseTransientResidualAssembler(const TimeStep* time, const std::vector<GlobalVectorType*>* u0, const std::vector<GlobalVectorType*>* u1, IElementWiseTransientResidualLocalAssembler* a)
        : _transient_e_assembler(a), _timestep(time), _u0(u0), _u1(u1)
    { };


    /// Conduct the element by element assembly procedure
    ///
    /// @param msh                 Mesh
    /// @param dofManager         Dof map manager
    /// @param list_dofId         List of Dof IDs used in this problem
    /// @param r                 Residual
    virtual void assembly( const MeshLib::IMesh &msh, const DiscreteLib::DofEquationIdTable &dofManager, GlobalVectorType &r);

private:
    IElementWiseTransientResidualLocalAssembler* _transient_e_assembler;
    const TimeStep* _timestep;
    const std::vector<GlobalVectorType*>* _u0;
    const std::vector<GlobalVectorType*>* _u1;
};


}
