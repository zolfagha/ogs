/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DeformationInPorousMediaLinearEQSLocalAssembler.h
 *
 * Created on 2012-11-15 by Norihiro Watanabe
 */


#pragma once

#include "MeshLib/Core/IElement.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Utils/Tools.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTransientCoupledLinearEQSLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"

namespace THMmf
{

class DeformationInPorousMediaLinearEQSLocalAssembler
: public NumLib::ElementWiseTransientCoupledLinearEQSLocalAssembler
{
public:
    typedef MathLib::LocalVector LocalVectorType;
    typedef MathLib::LocalMatrix LocalMatrixType;

    explicit DeformationInPorousMediaLinearEQSLocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects, size_t n_var, const std::vector<size_t> &vec_order)
    : NumLib::ElementWiseTransientCoupledLinearEQSLocalAssembler(n_var, vec_order), _feObjects(*feObjects)
    {
    };

    virtual ~DeformationInPorousMediaLinearEQSLocalAssembler() {};

protected:
    virtual void assembleComponents(  const NumLib::TimeStep &/*time*/,  
                            const MeshLib::IElement &e, 
                            const std::vector<size_t> &vec_order, 
                            const std::vector<LocalVectorType> &vec_u0, 
                            const std::vector<LocalVectorType> &vec_u1, 
                            std::vector<std::vector<LocalMatrixType> > &vec_K,
                            std::vector<LocalVectorType> &vec_F
                            );

private:
    FemLib::LagrangianFeObjectContainer _feObjects;
};

}
