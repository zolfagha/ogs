
#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "Tests/Geo/Material/PorousMedia.h"

namespace Geo
{

template <class T>
class GroundwaterFlowTimeODELocalAssembler: public T
{
private:
	PorousMedia* _pm;
	FemLib::LagrangianFeObjectContainer* _feObjects;
public:
    typedef typename T::LocalVectorType LocalVectorType;
    typedef typename T::LocalMatrixType LocalMatrixType;

	GroundwaterFlowTimeODELocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, PorousMedia &pm)
	: _pm(&pm), _feObjects(&feObjects)
	{
	};

	virtual ~GroundwaterFlowTimeODELocalAssembler() {};

protected:
	virtual void assembleODE(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/, LocalMatrixType &/*localM*/, LocalMatrixType &localK, LocalVectorType &/*localF*/)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

        	double k;
        	_pm->hydraulic_conductivity->eval(real_x, k);

    		//fe->integrateWxN(_pm->storage, localM);
    		fe->integrateDWxDN(j, k, localK);
        }
	}
};

class GroundwaterFlowJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
private:
	PorousMedia* _pm;
	FemLib::LagrangianFeObjectContainer* _feObjects;
public:
	GroundwaterFlowJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, PorousMedia &pm)
	: _pm(&pm), _feObjects(&feObjects)
	{
	};

	void assembly(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/,  LocalMatrixType &localJ)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

        	double k;
        	_pm->hydraulic_conductivity->eval(real_x, k);

    		//fe->integrateWxN(_pm->storage, localM);
    		fe->integrateDWxDN(j, k, localJ);
        }
	}
};


} //end
