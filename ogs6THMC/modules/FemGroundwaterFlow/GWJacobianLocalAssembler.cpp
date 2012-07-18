
#include "GWJacobianLocalAssembler.h"

#include "NumLib/Function/TXFunction.h"
#include "Ogs6FemData.h"

namespace Geo
{

void GroundwaterFlowJacobianLocalAssembler::assembly(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const NumLib::LocalVector &/*u1*/, const NumLib::LocalVector &/*u0*/,  NumLib::LocalMatrix &localJ)
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

        NumLib::LocalMatrix k;
        pm->hydraulic_conductivity->eval(real_x, k);

        //fe->integrateWxN(_pm->storage, localM);
        fe->integrateDWxDN(j, k, localJ);
    }
}


} //end
