/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseVectorAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Core/LocalDataType.h"
#include "DiscreteLib/Core/IDiscreteVectorAssembler.h"
#include "DiscreteLib/Core/IElemenetWiseVectorLocalAssembler.h"

namespace MeshLib 
{
class IMesh;
}

namespace DiscreteLib
{

/**
 * \brief Element-based discrete vector assembler classes
 */
template <class T_VALUE, class T_UPDATER>
class ElementWiseVectorAssembler : public IDiscreteVectorAssembler<T_VALUE>
{
public:
    typedef typename IDiscreteVectorAssembler<T_VALUE>::VectorType GlobalVectorType;
    typedef T_UPDATER UpdaterType;

    ///
    explicit ElementWiseVectorAssembler(UpdaterType* a) : _e_assembler(a) {};
    ///
    virtual ~ElementWiseVectorAssembler() {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param vec Discrete vector
    virtual void assembly(const MeshLib::IMesh &msh, GlobalVectorType &globalVec);
    //void assembly(const MeshLib::IMesh &msh, const DofEquationIdTable &dofManager, GlobalVectorType &globalVec);

private:
    UpdaterType* _e_assembler;
};


template <class T1, class T2>
void ElementWiseVectorAssembler<T1,T2>::assembly(const MeshLib::IMesh &msh, GlobalVectorType &globalVec)
{
    const size_t n_ele = msh.getNumberOfElements();

    for (size_t i=0; i<n_ele; i++) {
        MeshLib::IElement *e = msh.getElemenet(i);
        _e_assembler->update(*e, globalVec);
    }
};

}
