
#include "AsyncPartSolution.h"

#include <algorithm>

namespace NumLib
{
void AsyncPartSolution::setAlgorithm(PartitionedAlgorithm *algo) 
{
	algorithm = algo;
}

void AsyncPartSolution::addChildren(ITransientSystem* sys)
{
    listChildren.push_back(sys);
}

size_t AsyncPartSolution::getNumberOfChildren() const
{
    return listChildren.size();
}

TimeStep AsyncPartSolution::suggestNext(TimeStep time_current) {
    TimeStep t = .0; 
    for (size_t i=0; i<listChildren.size(); i++) {
        ITransientSystem *solution = listChildren[i];
        t = std::min(t, solution->suggestNext(time_current));
    }
	return t;
}
	
bool AsyncPartSolution::solveNextStep(TimeStep time) {
	for (size_t i=0; i<listChildren.size(); i++) {
        ITransientSystem *solution = listChildren[i];
		if (solution->isAwake(time)) {
			if (!solution->solveNextStep(time)) 
				return false;
		}
	}
	return true;
};
	
bool AsyncPartSolution::isAwake(TimeStep time) {
    for (size_t i=0; i<listChildren.size(); i++) {
        ITransientSystem *solution = listChildren[i];
        if (solution->isAwake(time))
            return true;
    }
    return false;
};
}
