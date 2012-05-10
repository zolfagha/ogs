
#pragma once

#include <vector>
#include <string>

#include "Base/Options.h"

#include "ICoupledProblem.h"
#include "MonolithicProblem.h"
#include "PartitionedProblem.h"
#include "Algorithm/PartitionedAlgorithmFactory.h"


namespace NumLib
{

template <class T_I, class T_M, class T_P, class T_ALGORITHM>
class TemplateCouplingStrucutreBuilder
{
public:
	template <class T_EQS_FACTORY, class T_CHECK_FACTORY>
	T_I* build(const Base::Options *option, T_EQS_FACTORY eqs_fac, T_CHECK_FACTORY check_fack)
	{
		const Base::Options* op_cpl = option->getSubGroup("coupling");
	    if (op_cpl->hasSubGroup("M")) {
	        const Base::Options* op_sub = op_cpl->getSubGroup("M");
	        T_M *sys = buildMonolithicSystem(op_sub, eqs_fac);
	        return sys;
	    } else if (op_cpl->hasSubGroup("P")) {
	        const Base::Options* op_sub = op_cpl->getSubGroup("P");
	        T_P *sys = buildPartitionedSystem(op_sub, eqs_fac, check_fack);
	        return sys;
	    }
	    return 0;
	}

private:
	template <class T_EQS_FACTORY>
	T_M* buildMonolithicSystem(const Base::Options *option, T_EQS_FACTORY &eqs_fac)
	{
		T_M* eqs = eqs_fac.create(option->getOption("name"));
		const std::vector<std::string>* in_names = option->getOptionAsArray<std::string>("in");
		const std::vector<std::string>* out_names = option->getOptionAsArray<std::string>("out");
		if (in_names!=0) {
			for (size_t i=0; i<in_names->size(); i++) {
				eqs->setInputParameterName(i, (*in_names)[i]);
			}
		}
		if (out_names!=0) {
			for (size_t i=0; i<out_names->size(); i++) {
				eqs->setOutputParameterName(i, (*out_names)[i]);
			}
		}
		return eqs;
	}

	template <class T_EQS_FACTORY, class T_CHECK_FACTORY>
	T_P* buildPartitionedSystem(const Base::Options *option, T_EQS_FACTORY &eqs_fac, T_CHECK_FACTORY check_fac)
	{
		T_P* part = new T_P();
		//para
		const std::vector<std::string>* in_names = option->getOptionAsArray<std::string>("in");
		const std::vector<std::string>* out_names = option->getOptionAsArray<std::string>("out");
		if (in_names!=0) {
	        part->resizeInputParameter(in_names->size());
			for (size_t i=0; i<in_names->size(); i++) {
				part->setInputParameterName(i, (*in_names)[i]);
			}
		}
		if (out_names!=0) {
	        part->resizeOutputParameter(out_names->size());
			for (size_t i=0; i<out_names->size(); i++) {
				part->setOutputParameterName(i, (*out_names)[i]);
			}
		}
		//alg
		size_t max_itr = option->getOption<size_t>("max_itr");
		double epsilon = option->getOption<double>("epsilon");
		IConvergenceCheck* checker = check_fac.create(option->getOption("convergence"));
		part->setAlgorithm(*T_ALGORITHM::create(option->getOption("algorithm"), checker, max_itr, epsilon));
//		IPartitionedAlgorithm* alg = T_ALGORITHM::create(option->getOption("algorithm"), checker, max_itr, epsilon);
//		if (alg!=0) {
//		}
		//problems
		const Base::Options* op_problems = option->getSubGroup("problems");
		for (Base::Options::const_iterator itr=op_problems->begin(); itr!=op_problems->end(); ++itr) {
			std::string str = itr->first;
			Base::Options* op_sub = static_cast<Base::Options*>(itr->second);
			T_I* sys = 0;
			if (str.find("M")==0) {
				sys = buildMonolithicSystem(op_sub, eqs_fac);
			} else if (str.compare("P")==0) {
				sys = buildPartitionedSystem(op_sub, eqs_fac, check_fac);
			}
	        if (sys!=0) part->addProblem(*sys);
		}
	    part->connectParameters();
		return part;
	}
};

typedef class TemplateCouplingStrucutreBuilder<ICoupledSystem,TemplateSteadyMonolithicSystem,PartitionedProblem,PartitionedAlgorithmFactory> CouplingStrucutreBuilder;


} //end
