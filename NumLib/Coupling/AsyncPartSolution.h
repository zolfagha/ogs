
#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

#include "CouplingSolution.h"
#include "TransientSystems.h"
#include "PartitionedAlgorithm.h"

namespace NumLib
{




/**
 * \brief 
 */
class AsyncPartSolution : public ITransientSystem
{
private:
	IPartitionedAlgorithm *algorithm;
	std::vector<ITransientSystem*> _listChildren;
    std::vector<TimeStep> _list_synchronize_time;
	
public:
    void addChildren(ITransientSystem* sys);
    size_t getNumberOfChildren() const;
	
	TimeStep suggestNext(TimeStep time_current);
	
	bool solveNextStep(TimeStep time);
	
	bool isAwake(TimeStep time);
};

}