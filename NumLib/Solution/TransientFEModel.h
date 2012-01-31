
#pragma once

#include "NumLib/Core/TimeStep.h"
#include "NumLib/Coupling/TransientSystems.h"

namespace NumLib
{

class TransientFEModel: public ITransientSystem
{
private:
	//Geo *GEO;
	//Mesh *MSH;
	//Mat *MAT;
	//TimeDisc TIM;
	//Fe FE;
	//Ic IC;
	//Bc BC;
	//St ST;
	//Num NUM;
    TimeStepping *TIM;

public:	
	//void initialize(Input*) {
	//	//GEO, MSH, TIM, IC, BC, ST, NUM
	//}
	
	TimeStep suggestNext(TimeStep t_n) {
		return TIM->next();
	}
	bool isAwake(TimeStep t_n1) {
		return (TIM->next()==t_n1);
	}
	bool solveNextStep(TimeStep time) {
		return solveTimeStep(time);
	};
	virtual void setupModel(int*) = 0;
	virtual bool solveTimeStep(TimeStep t_n1) = 0;
};

}
