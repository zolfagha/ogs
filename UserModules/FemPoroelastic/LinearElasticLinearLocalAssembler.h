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
#include "NumLib/TimeStepping/TimeStep.h"

class TransientCoupledProlbemLocalAssembler : public NumLib::IElementWiseTransientLinearEQSLocalAssembler
{
public:
    typedef NumLib::LocalVector LocalVectorType;
    typedef NumLib::LocalMatrix LocalMatrixType;

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
                            )
    {
        // parameters need to be passed
        const size_t n_var = _n_var;
        std::vector<size_t> &vec_order = _vec_order;

        //
        std::vector<std::vector<size_t> > vec_local_pos(n_var);
        for (size_t i=0; i<n_var; i++) {
            std::vector<size_t> list_nodeid;
            e.getNodeIDList(vec_order[i], list_nodeid);
            localDofManager.mapEqsID(0, 0, list_nodeid, vec_local_pos[i]);
        }
        
        std::vector<LocalVectorType> vec_u0(n_var), vec_u1(n_var);
        for (size_t i=0; i<n_var; i++) {
            DiscreteLib::getLocalVector(vec_local_pos[i], local_u_n, vec_u0[i]);
            DiscreteLib::getLocalVector(vec_local_pos[i], local_u_n1, vec_u1[i]);
        }
        
        std::vector<std::vector<LocalMatrixType> > vec_K(n_var);
        std::vector<LocalVectorType> vec_F(n_var);
        
        this->assemble(timestep, e, vec_order, vec_u0, vec_u1, vec_K, vec_F);
        
        for (size_t i=0; i<n_var; i++) {
            for (size_t j=0; j<n_var; j++) {
                eqs.addAsub(vec_local_pos[i], vec_local_pos[j], vec_K[i][j]);
            }
            eqs.addRHSsub(vec_local_pos[i], vec_F[i]);
        }
    }
    
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


class FemPoroelasticLinearLocalAssembler
: public TransientCoupledProlbemLocalAssembler
{
public:
    typedef NumLib::LocalVector LocalVectorType;
    typedef NumLib::LocalMatrix LocalMatrixType;

    explicit FemPoroelasticLinearLocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects, size_t n_var, const std::vector<size_t> &vec_order)
    : TransientCoupledProlbemLocalAssembler(n_var, vec_order), _feObjects(*feObjects)
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

