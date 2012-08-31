/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PressureBasedGWTimeODELocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Fluid.h"
#include "Ogs6FemData.h"

/**
 * \brief Local assembly class for time-ODE GW equation
 */
template <class T>
class PressureBasedGWTimeODELocalAssembler: public T
{
public:
    typedef NumLib::LocalVector LocalVector;
    typedef NumLib::LocalMatrix LocalMatrix;

    explicit PressureBasedGWTimeODELocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects)
    : _feObjects(feObjects)
    {
    };

    virtual ~PressureBasedGWTimeODELocalAssembler() {};

protected:
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const LocalVector &/*u1*/, const LocalVector &/*u0*/, LocalMatrix &localM, LocalMatrix &localK, LocalVector &/*localF*/)
    {
        FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());

        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        MaterialLib::Fluid* fluid = Ogs6FemData::getInstance()->list_fluid[0];
        double mu = .0;
        fluid->dynamic_viscosity->eval(e_pos, mu);
        double rho_f = .0;
        fluid->density->eval(e_pos, rho_f);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            LocalMatrix &Np = *fe->getBasisFunction();
            LocalMatrix &dNp = *fe->getGradBasisFunction();
            fe->getRealCoordinates(real_x);
            double fac = fe->getDetJ() * q->getWeight(j);

            double k;
            pm->permeability->eval(real_x, k);
            double s;
            pm->storage->eval(real_x, s);
            double k_mu;
            k_mu = k / mu;

            // M_pp += Np^T * S * Np
            localM.noalias() += fac * Np.transpose() * s * Np;

            // K_pp += dNp^T * K * dNp
            localK.noalias() += fac * dNp.transpose() * k_mu * dNp;
        }
    }

private:
    FemLib::LagrangianFeObjectContainer _feObjects;
};

