/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseLinearEquationAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "DiscreteLib/Core/IDiscreteLinearEquationAssembler.h"
#include "DiscreteLib/Core/IElemenetWiseLinearEquationLocalAssembler.h"

namespace MeshLib 
{
class IMesh;
}

namespace MathLib
{
    class ILinearEquation;
}

namespace DiscreteLib
{

/**
 * \brief Element-based discrete system assembler classes
 */
template <class T_UPDATER>
class OMPElementWiseLinearEquationAssembler : public IDiscreteLinearEquationAssembler
{
public:
    typedef T_UPDATER UpdaterType;

    ///
    explicit OMPElementWiseLinearEquationAssembler(UpdaterType* a) : _e_assembler(a) {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::ILinearEquation &eqs);

private:
    UpdaterType* _e_assembler;
};

template <class T>
void OMPElementWiseLinearEquationAssembler<T>::assembly(MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::ILinearEquation &eqs)
{
    const size_t n_ele = msh.getNumberOfElements();
    for (size_t i=0; i<n_ele; i++) {
        MeshLib::IElement *e = msh.getElemenet(i);
        _e_assembler->update(*e, eqs);
    }
};

}
