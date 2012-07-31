
#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "Ogs6FemData.h"

/**
 * \brief Local assembly class for time-ODE GW equation
 */
template <class T>
class GroundwaterFlowTimeODELocalAssembler: public T
{
public:
    typedef NumLib::LocalVector LocalVector;
    typedef NumLib::LocalMatrix LocalMatrix;

//    GroundwaterFlowTimeODELocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, MaterialLib::PorousMedia &pm)
//    : _pm(&pm), _feObjects(&feObjects)
//    {
//    };

    explicit GroundwaterFlowTimeODELocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects)
    : _feObjects(&feObjects)
    {
    };

    virtual ~GroundwaterFlowTimeODELocalAssembler() {};

protected:
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const LocalVector &/*u1*/, const LocalVector &/*u0*/, LocalMatrix &localM, LocalMatrix &localK, LocalVector &/*localF*/)
    {
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

            LocalMatrix k;
            pm->hydraulic_conductivity->eval(real_x, k);
            LocalMatrix s;
            pm->storage->eval(real_x, s);

            fe->integrateWxN(j, s, localM);
            fe->integrateDWxDN(j, k, localK);
        }
    }

private:
    //MaterialLib::PorousMedia* _pm;
    FemLib::LagrangianFeObjectContainer* _feObjects;
};

