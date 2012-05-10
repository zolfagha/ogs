
#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "NumLib/Function/Function.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Material/Compound.h"

namespace Geo
{

template <class T>
class MassTransportTimeODELocalAssembler: public T
{
public:
    typedef typename T::LocalVectorType LocalVectorType;
    typedef typename T::LocalMatrixType LocalMatrixType;

	MassTransportTimeODELocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, PorousMedia &pm, Compound &cmp)
	: _pm(&pm), _cmp(&cmp), _feObjects(&feObjects)
	{
        Base::zeroObject(_vel);
	};

	virtual ~MassTransportTimeODELocalAssembler() {};

    void initialize(const NumLib::SpatialFunctionVector *vel)
    {
        _vel = const_cast<NumLib::SpatialFunctionVector*>(vel);
    }

protected:
	virtual void assembleODE(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/, LocalMatrixType &localM, LocalMatrixType &localK, LocalVectorType &/*localF*/)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        LocalMatrixType matDiff(localK);
        LocalMatrixType matAdv(localK);
        NumLib::SpatialCompositFunction<double, NumLib::SpatialFunctionScalar, NumLib::SpatialFunctionScalar, NumLib::Multiplication> f_diff_poro(*_cmp->molecular_diffusion, *_pm->porosity);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

        	double poro;
        	_pm->porosity->eval(real_x, poro);
        	double d_poro;
        	f_diff_poro.eval(real_x, d_poro);
        	MathLib::Vector v;
        	_vel->eval(real_x, v);

     		fe->integrateWxN(j, poro, localM);
    		fe->integrateDWxDN(j, d_poro, matDiff);
            fe->integrateWxDN(j, v, matAdv);
        }

        localK = matDiff;
        localK += matAdv;

        //std::cout << "M="; localM.write(std::cout); std::cout << std::endl;
        //std::cout << "L="; matDiff.write(std::cout); std::cout << std::endl;
        //std::cout << "A="; matAdv.write(std::cout); std::cout << std::endl;
	}

private:
    PorousMedia* _pm;
    Compound* _cmp;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    NumLib::SpatialFunctionVector* _vel;
};



class MassTransportJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
	MassTransportJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, PorousMedia &pm, Compound &cmp)
	: _pm(&pm), _cmp(&cmp), _feObjects(&feObjects)
	{
        Base::zeroObject(_vel);
	};

    void initialize(const NumLib::SpatialFunctionVector *vel)
    {
        _vel = const_cast<NumLib::SpatialFunctionVector*>(vel);
    }

	void assembly(const NumLib::TimeStep &time, MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/, LocalMatrixType &localJ)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        LocalMatrixType matM(localJ);
        LocalMatrixType matDiff(localJ);
        LocalMatrixType matAdv(localJ);
        NumLib::SpatialCompositFunction<double, NumLib::SpatialFunctionScalar, NumLib::SpatialFunctionScalar, NumLib::Multiplication> f_diff_poro(*_cmp->molecular_diffusion, *_pm->porosity);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

        	double poro;
        	_pm->porosity->eval(real_x, poro);
        	double d_poro;
        	f_diff_poro.eval(real_x, d_poro);
        	MathLib::Vector v;
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
    PorousMedia* _pm;
    Compound* _cmp;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    NumLib::SpatialFunctionVector* _vel;
};
} //end
