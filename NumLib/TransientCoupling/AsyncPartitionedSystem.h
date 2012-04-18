
#pragma once

#include <vector>

#include "MathLib/Coupling/ICoupledProblem.h"
#include "MathLib/Coupling/MonolithicProblem.h"
#include "MathLib/Coupling/PartitionedProblem.h"
#include "MathLib/Coupling/Algorithm/TransientPartitionedAlgorithm.h"
#include "TransientCoupledSystem.h"

namespace NumLib
{

/**
 * \brief Asynchronized partitioned problem
 */
class AsyncPartitionedSystem : public MathLib::AbstractPartitionedProblem<ITransientCoupledSystem>
{
public:
	typedef size_t InternalID;
	typedef size_t ExternalKey;

    AsyncPartitionedSystem() : _algorithm(0)
    {
    }

    virtual ~AsyncPartitionedSystem() {};

    virtual void setAlgorithm(MathLib::ITransientPartitionedAlgorithm &algo)
    {
    	_algorithm = &algo;
    }

	double suggestNext(const TimeStep &time_current);
	
	int solveTimeStep(const TimeStep &time);
	
	bool isAwake(const TimeStep &time);

    void accept(const TimeStep &time);

    void getActiveProblems(const TimeStep &time, std::vector<MathLib::ICoupledSystem*> &list_active_problems);

private:
    MathLib::ITransientPartitionedAlgorithm *_algorithm;
    std::vector<TimeStep> _list_synchronize_time;
	
    MathLib::UnnamedParameterSet _vars_t_n;
};


}
