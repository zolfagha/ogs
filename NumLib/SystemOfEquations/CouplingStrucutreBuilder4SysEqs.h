
#pragma once

#include <vector>
#include <string>

#include "Base/Options.h"

#include "NumLib/Coupling/ICoupledProblem.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/Coupling/PartitionedProblem.h"
#include "NumLib/Coupling/Algorithm/PartitionedAlgorithmFactory.h"
#include "MyCouplingEQS.h"

namespace NumLib
{

template <class T_I, class T_M, class T_P, class T_ALGORITHM>
class TemplateCouplingStrucutreBuilder4SysEqs
{
public:
	template <class T_EQS_FACTORY, class T_CHECK_FACTORY>
	T_I* build(const Base::Options *option, T_EQS_FACTORY eqs_fac, T_CHECK_FACTORY check_fack)
	{
		const Base::Options* op_cpl = option->getSubGroup("coupling");
        if (op_cpl->begin()==op_cpl->end()) return 0;

        std::string str = op_cpl->begin()->first;
        Base::Options* op_sub = static_cast<Base::Options*>(op_cpl->begin()->second);
        if (str.find("M")==0) {
            //const Base::Options* op_sub = op_cpl->getSubGroup("M");
            T_M *sys = buildMonolithicSystem(op_sub, eqs_fac);
            return sys;
        } else if (str.find("P")==0) {
            //const Base::Options* op_sub = op_cpl->getSubGroup("P");
            T_P *sys = buildPartitionedSystem(op_sub, eqs_fac, check_fack);
            return sys;
        }

	    //if (op_cpl->hasSubGroup("M")) {
	    //    const Base::Options* op_sub = op_cpl->getSubGroup("M");
	    //    T_M *sys = buildMonolithicSystem(op_sub, eqs_fac);
	    //    return sys;
	    //} else if (op_cpl->hasSubGroup("P")) {
	    //    const Base::Options* op_sub = op_cpl->getSubGroup("P");
	    //    T_P *sys = buildPartitionedSystem(op_sub, eqs_fac, check_fack);
	    //    return sys;
	    //}
	    return 0;
	}

private:
	template <class T_EQS_FACTORY>
	T_M* buildMonolithicSystem(const Base::Options *option, T_EQS_FACTORY &eqs_fac)
	{
		const std::vector<std::string>* list_var_name = option->getOptionAsArray<std::string>("variable");
		T_M* eqs = eqs_fac.create(*list_var_name);
		return eqs;
	}

	template <class T_EQS_FACTORY, class T_CHECK_FACTORY>
	T_P* buildPartitionedSystem(const Base::Options *option, T_EQS_FACTORY &eqs_fac, T_CHECK_FACTORY check_fac)
	{
		T_P* part = new T_P();
		//alg
		size_t max_itr = option->getOption<size_t>("max_itr");
		double epsilon = option->getOption<double>("epsilon");
		IConvergenceCheck* checker = check_fac.create(option->getOption("convergence"));
		part->setAlgorithm(*T_ALGORITHM::create(option->getOption("algorithm"), checker, max_itr, epsilon));
		//problems
		const Base::Options* op_problems = option->getSubGroup("problems");
		std::vector<T_I*> list_subproblems;
		for (Base::Options::const_iterator itr=op_problems->begin(); itr!=op_problems->end(); ++itr) {
			std::string str = itr->first;
			Base::Options* op_sub = static_cast<Base::Options*>(itr->second);
			T_I* sys = 0;
			if (str.find("M")==0) {
				sys = buildMonolithicSystem(op_sub, eqs_fac);
			} else if (str.compare("P")==0) {
				sys = buildPartitionedSystem(op_sub, eqs_fac, check_fac);
			}
	        if (sys!=0) {
	        	part->addProblem(*sys);
	        }
		}
		eqs_fac.setPartitionedProblem(part);
		return part;
	}
};

typedef class TemplateCouplingStrucutreBuilder4SysEqs<ICoupledSystem,MyCouplingEQS,PartitionedProblem,PartitionedAlgorithmFactory> CouplingStrucutreBuilder4SysEqs;


} //end
