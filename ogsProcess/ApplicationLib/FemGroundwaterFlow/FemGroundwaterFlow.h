
#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"

namespace Geo
{

template <class T>
class GroundwaterFlowTimeODELocalAssembler: public T
{
private:
	MaterialLib::PorousMedia* _pm;
	FemLib::LagrangianFeObjectContainer* _feObjects;
public:
    typedef NumLib::LocalVector LocalVector;
    typedef NumLib::LocalMatrix LocalMatrix;

	GroundwaterFlowTimeODELocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, MaterialLib::PorousMedia &pm)
	: _pm(&pm), _feObjects(&feObjects)
	{
	};

	virtual ~GroundwaterFlowTimeODELocalAssembler() {};

protected:
	virtual void assembleODE(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const LocalVector &/*u1*/, const LocalVector &/*u0*/, LocalMatrix &/*localM*/, LocalMatrix &localK, LocalVector &/*localF*/)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

        	LocalMatrix k;
        	_pm->hydraulic_conductivity->eval(real_x, k);

    		//fe->integrateWxN(_pm->storage, localM);
    		fe->integrateDWxDN(j, k, localK);
        }
	}
};

class GroundwaterFlowJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
private:
	MaterialLib::PorousMedia* _pm;
	FemLib::LagrangianFeObjectContainer* _feObjects;
public:
	GroundwaterFlowJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, MaterialLib::PorousMedia &pm)
	: _pm(&pm), _feObjects(&feObjects)
	{
	};

	void assembly(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const NumLib::LocalVector &/*u1*/, const NumLib::LocalVector &/*u0*/,  NumLib::LocalMatrix &localJ)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

        	NumLib::LocalMatrix k;
        	_pm->hydraulic_conductivity->eval(real_x, k);

    		//fe->integrateWxN(_pm->storage, localM);
    		fe->integrateDWxDN(j, k, localJ);
        }
	}
};


} //end
