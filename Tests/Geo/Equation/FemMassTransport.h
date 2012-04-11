
#pragma once

#include "MathLib/Function/IFunction.h"
#include "FemLib/Core/IFemElement.h"
#include "FemLib/Core/Integration.h"
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

 		fe->integrateWxN(_pm->porosity, localM);

        LocalMatrixType matDiff(localK);
        LocalMatrixType matAdv(localK);
        SpatialCompositFunction<double, MathLib::SpatialFunctionScalar, MathLib::SpatialFunctionScalar, MathLib::Multiplication> f_diff_poro(*_cmp->molecular_diffusion, *_pm->porosity);
		fe->integrateDWxDN(&f_diff_poro, matDiff);
        fe->integrateWxDN(_vel, matAdv);
        localK = matDiff;
        localK += matAdv;
	}

private:
    PorousMedia* _pm;
    Compound* _cmp;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    MathLib::SpatialFunctionVector* _vel;
};

} //end
