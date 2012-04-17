
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

//    /// find subproblem
//    int find(const ITransientCoupledSystem &sub) const;
//
//    /// check consistency
//    bool check() const;
//
//
//    /// add parameter without giving reference
//    size_t addInputParameter(const std::string &name);
//
//    /// add parameter and reference
//    /// @param name variable name
//    /// @param sys problem
//    /// @param internal_id parameter id in the sys
//    /// @return parameter id
//    size_t addOutputParameter(const std::string &name, ITransientCoupledSystem & sub_problem, size_t para_id_in_sub_problem);
//
//    /// connect system input and shared variable
//    void connectInput(const std::string &this_para_name, ITransientCoupledSystem &subproblem, size_t subproblem_para_id);

	double suggestNext(const TimeStep &time_current);
	
	int solveTimeStep(const TimeStep &time);
	
	bool isAwake(const TimeStep &time);

    void accept(const TimeStep &time);

    void getActiveProblems(const TimeStep &time, std::vector<MathLib::ICoupledSystem*> &list_active_problems);

private:
    MathLib::ITransientPartitionedAlgorithm *_algorithm;
//    std::vector<ITransientCoupledSystem*> _list_subproblems;
    std::vector<TimeStep> _list_synchronize_time;
	
    MathLib::UnnamedParameterSet _vars_t_n;
    //MathLib::ParameterProblemMappingTable _map;

    //size_t addSubProblem(ITransientCoupledSystem &sub_problem);
};


}
