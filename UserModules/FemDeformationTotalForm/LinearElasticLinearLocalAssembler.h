/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearElasticLinearLocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


#pragma once

#include "MeshLib/Core/IElement.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"

class FemLinearElasticLinearLocalAssembler: public NumLib::IElementWiseTransientLinearEQSLocalAssembler
{
public:
    typedef NumLib::LocalVector LocalVectorType;
    typedef NumLib::LocalMatrix LocalMatrixType;

    explicit FemLinearElasticLinearLocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects)
    : _feObjects(*feObjects)
    {
    };

    virtual ~FemLinearElasticLinearLocalAssembler() {};

    virtual void assembly(const NumLib::TimeStep &/*time*/,  const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &localDofManager, const LocalVectorType &/*local_u_n1*/, const LocalVectorType &/*local_u_n*/, NumLib::LocalEquation &eqs);

private:
    FemLib::LagrangianFeObjectContainer _feObjects;
};

