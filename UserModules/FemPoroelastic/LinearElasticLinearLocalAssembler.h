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
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Utils/Tools.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/TransientCoupledProlbemLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"



class FemPoroelasticLinearLocalAssembler
: public NumLib::TransientCoupledProlbemLocalAssembler
{
public:
    typedef NumLib::LocalVector LocalVectorType;
    typedef NumLib::LocalMatrix LocalMatrixType;

    explicit FemPoroelasticLinearLocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects, size_t n_var, const std::vector<size_t> &vec_order)
    : NumLib::TransientCoupledProlbemLocalAssembler(n_var, vec_order), _feObjects(*feObjects)
    {
    };

    virtual ~FemPoroelasticLinearLocalAssembler() {};

protected:
    //virtual void assembly(const NumLib::TimeStep &/*time*/,  const MeshLib::IElement &e, const LocalVectorType &/*local_u_n1*/, const LocalVectorType &/*local_u_n*/, NumLib::LocalEquation &eqs);
    virtual void assemble(  const NumLib::TimeStep &/*time*/,  
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

