/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearTransportJacobianLocalAssembler.h
 *
 * Created on 2012-09-19 by Haibing Shao
 */

/**
  * This file is same as the MassTransportTimeODELocalAssembler.h
  * The difference is, the compound molecular diffusion coefficient is disabled, 
  */

#ifndef LINEAR_TRANSPORT_JACOBIAN_LOCAL_ASSEMBLER_H
#define LINEAR_TRANSPORT_JACOBIAN_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"

#include "Ogs6FemData.h"

class LinearTransportJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    LinearTransportJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects)
        : _feObjects(*feObjects), _vel(NULL)
    {
    };

    virtual ~LinearTransportJacobianLocalAssembler() {};

    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

    void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const NumLib::LocalVector &/*u1*/, const NumLib::LocalVector &/*u0*/, NumLib::LocalMatrix &localJ)
    {
        FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        double cmp_mol_diffusion = .0;
        // _cmp->molecular_diffusion->eval(0, cmp_mol_diffusion);

        NumLib::LocalMatrix matM(localJ);
        NumLib::LocalMatrix matDiff(localJ);
        NumLib::LocalMatrix matAdv(localJ);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        NumLib::LocalMatrix poro(1,1);
        NumLib::LocalMatrix d_poro(1,1);
        NumLib::ITXFunction::DataType v;

        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

            pm->porosity->eval(real_x, poro);
            d_poro(0,0) = cmp_mol_diffusion * poro(0,0);
            _vel->eval(real_x, v);

            fe->integrateWxN(j, poro, matM);
            fe->integrateDWxDN(j, d_poro, matDiff);
            fe->integrateWxDN(j, v, matAdv);
        }

        double dt = time.getTimeStepSize();
        double theta = 1.0;
        matM *= 1.0 / dt;
        matDiff *= theta;
        matAdv *= theta;
        localJ = matM;
        localJ += matDiff;
        localJ += matAdv;

        //std::cout << "M="; localM.write(std::cout); std::cout << std::endl;
        //std::cout << "L="; matDiff.write(std::cout); std::cout << std::endl;
        //std::cout << "A="; matAdv.write(std::cout); std::cout << std::endl;
    }

private:
    FemLib::LagrangianFeObjectContainer _feObjects;
    NumLib::ITXFunction* _vel;
};

#endif  // end of ifndef
