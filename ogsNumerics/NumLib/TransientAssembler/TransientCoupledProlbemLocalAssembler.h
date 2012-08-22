/* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TransientCoupledProlbemLocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "MeshLib/Core/IElement.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "NumLib/DataType.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "IElementWiseTransientLinearEQSLocalAssembler.h"

namespace NumLib
{

class TransientCoupledProlbemLocalAssembler
: public IElementWiseTransientLinearEQSLocalAssembler
{
public:
    typedef LocalVector LocalVectorType;
    typedef LocalMatrix LocalMatrixType;

    TransientCoupledProlbemLocalAssembler(size_t n_var, const std::vector<size_t> &vec_order)
    : _n_var(n_var), _vec_order(vec_order)
    {
    };

    virtual ~TransientCoupledProlbemLocalAssembler() {};

    virtual void assembly(  const NumLib::TimeStep &timestep,
                            const MeshLib::IElement &e,
                            const DiscreteLib::DofEquationIdTable &localDofManager,
                            const LocalVectorType &local_u_n1,
                            const LocalVectorType &local_u_n,
                            NumLib::LocalEquation &eqs
                            );


protected:
    virtual void assemble(  const NumLib::TimeStep &/*time*/,
                            const MeshLib::IElement &e,
                            const std::vector<size_t> &vec_order,
                            const std::vector<LocalVectorType> &vec_u0,
                            const std::vector<LocalVectorType> &vec_u1,
                            std::vector<std::vector<LocalMatrixType> > &vec_K,
                            std::vector<LocalVectorType> &vec_F
                            ) = 0;

private:
    size_t _n_var;
    std::vector<size_t> _vec_order;
};

} //end
