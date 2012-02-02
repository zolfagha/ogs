
#pragma once

#include <vector>

#include "TransientSystems.h"
#include "PartitionedAlgorithm.h"

namespace NumLib
{
class AsyncPartSolution : public ITransientSystem
{
private:
	PartitionedAlgorithm *algorithm;
	std::vector<ITransientSystem*> listChildren;
	
public:
	void setAlgorithm(PartitionedAlgorithm *algo);

    void addChildren(ITransientSystem* sys);
    size_t getNumberOfChildren() const;
	
	TimeStep suggestNext(TimeStep time_current);
	
	bool solveNextStep(TimeStep time);
	
	bool isAwake(TimeStep time);
};
}