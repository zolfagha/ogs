
#pragma once

#include "FemLib/Core/IFemElement.h"
#include "Tests/Geo/Material/PorousMedia.h"

namespace Geo
{

class WeakFormGroundwaterFlow: public NumLib::ITimeODEElementAssembler
{
private:
	PorousMedia* _pm;
	FemLib::LagrangianFeObjectContainer* _feObjects;
public:
	WeakFormGroundwaterFlow(FemLib::LagrangianFeObjectContainer &feObjects, PorousMedia &pm)
	: _pm(&pm), _feObjects(&feObjects)
	{
	};

	//protected:
	void assembly(const NumLib::TimeStep &time, MeshLib::IElement &e, const LocalVectorType &u1, const LocalVectorType &u0, LocalMatrixType &localM, LocalMatrixType &localK, LocalVectorType &localF)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

		//fe->integrateWxN(_pm->storage, localM);
		fe->integrateDWxDN(_pm->hydraulic_conductivity, localK);
		//localF.resize(localF.size(), .0);
	}
};

} //end
