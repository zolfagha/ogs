
#pragma once

#include "MathLib/Function/IFunction.h"
#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Material/Compound.h"

namespace Geo
{

class WeakFormMassTransport: public NumLib::ITimeODEElementAssembler
{
public:
	WeakFormMassTransport(FemLib::LagrangianFeObjectContainer &feObjects, PorousMedia &pm, Compound &cmp)
	: _pm(&pm), _cmp(&cmp), _feObjects(&feObjects)
	{
        Base::zeroObject(_vel);
	};

    void initialize(const MathLib::SpatialFunctionVector *vel)
    {
        _vel = const_cast<MathLib::SpatialFunctionVector*>(vel);
    }

	void assembly(const NumLib::TimeStep &time, MeshLib::IElement &e, const LocalVectorType &u1, const LocalVectorType &u0, LocalMatrixType &localM, LocalMatrixType &localK, LocalVectorType &localF)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        LocalMatrixType matDiff(localK);
        LocalMatrixType matAdv(localK);
        SpatialCompositFunction<double, MathLib::SpatialFunctionScalar, MathLib::SpatialFunctionScalar, MathLib::Multiplication> f_diff_poro(*_cmp->molecular_diffusion, *_pm->porosity);

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
    MathLib::SpatialFunctionVector* _vel;
};

} //end
