
#pragma once

#include "FemLib/Core/IFemElement.h"

class GWAssembler: public NumLib::ITimeODEElementAssembler
{
private:
	MathLib::TemplateFunction<double*, double>* _matK;
	FemLib::LagrangianFeObjectContainer* _feObjects;
public:
	GWAssembler(FemLib::LagrangianFeObjectContainer &feObjects, MathLib::TemplateFunction<double*, double> &mat)
	: _matK(&mat), _feObjects(&feObjects)
	{
	};

	//protected:
	void assembly(const NumLib::TimeStep &time, MeshLib::IElement &e, const LocalVectorType &u1, const LocalVectorType &u0, LocalMatrixType &localM, LocalMatrixType &localK, LocalVectorType &localF)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

		//localM = .0;
		//localF.resize(localF.size(), .0);
		fe->integrateDWxDN(_matK, localK);
	}
};
