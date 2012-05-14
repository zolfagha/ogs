
#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "NumLib/Function/TXFunction.h"
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
    typedef NumLib::LocalVector LocalVectorType;
    typedef NumLib::LocalMatrix LocalMatrixType;

	MassTransportTimeODELocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, PorousMedia &pm, Compound &cmp)
	: _pm(&pm), _cmp(&cmp), _feObjects(&feObjects)
	{
        Base::zeroObject(_vel);
	};

	virtual ~MassTransportTimeODELocalAssembler() {};

    void initialize(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

protected:
	virtual void assembleODE(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/, LocalMatrixType &localM, LocalMatrixType &localK, LocalVectorType &/*localF*/)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);
        const size_t n_dim = e.getDimension();

        LocalMatrixType matDiff(localK);
        LocalMatrixType matAdv(localK);
        NumLib::TXCompositFunction<NumLib::ITXFunction, NumLib::ITXFunction, NumLib::Multiplication> f_diff_poro(*_cmp->molecular_diffusion, *_pm->porosity);

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
        	NumLib::ITXFunction::DataType v;
        	_vel->eval(real_x, v);
            NumLib::ITXFunction::DataType v2 = v.topRows(n_dim).transpose();
        	NumLib::LocalMatrix mat_poro(1,1);
        	mat_poro(0,0) = d_poro;

     		fe->integrateWxN(j, poro, localM);
    		fe->integrateDWxDN(j, mat_poro, matDiff);
            fe->integrateWxDN(j, v2, matAdv);
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
    NumLib::ITXFunction* _vel;
};



class MassTransportJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
	MassTransportJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, PorousMedia &pm, Compound &cmp)
	: _pm(&pm), _cmp(&cmp), _feObjects(&feObjects)
	{
        Base::zeroObject(_vel);
	};

    void initialize(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

	void assembly(const NumLib::TimeStep &time, MeshLib::IElement &e, const NumLib::LocalVector &/*u1*/, const NumLib::LocalVector &/*u0*/, NumLib::LocalMatrix &localJ)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

		NumLib::LocalMatrix matM(localJ);
		NumLib::LocalMatrix matDiff(localJ);
		NumLib::LocalMatrix matAdv(localJ);
        NumLib::TXCompositFunction<NumLib::ITXFunction, NumLib::ITXFunction, NumLib::Multiplication> f_diff_poro(*_cmp->molecular_diffusion, *_pm->porosity);

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
        	NumLib::ITXFunction::DataType v;
        	_vel->eval(real_x, v);
        	NumLib::LocalMatrix mat_poro(1,1);
        	mat_poro(0,0) = d_poro;

     		fe->integrateWxN(j, poro, matM);
    		fe->integrateDWxDN(j, mat_poro, matDiff);
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
    NumLib::ITXFunction* _vel;
};
} //end
