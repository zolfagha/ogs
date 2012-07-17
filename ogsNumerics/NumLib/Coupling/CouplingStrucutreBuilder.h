
#pragma once

#include <vector>
#include <string>

#include "BaseLib/Options.h"

#include "ICoupledProblem.h"
#include "MonolithicProblem.h"
#include "PartitionedProblem.h"
#include "Algorithm/PartitionedAlgorithmFactory.h"


namespace NumLib
{

template <class T_I, class T_M, class T_P, class T_ALGORITHM>
class TemplateCouplingStrucutreBuilder
{
	std::vector<T_M*> _vec_m;
	std::vector<std::string> _vec_m_name;
public:
	TemplateCouplingStrucutreBuilder() {};

	template <class T_EQS_FACTORY, class T_CHECK_FACTORY>
	T_I* build(const BaseLib::Options *option, T_EQS_FACTORY &eqs_fac, T_CHECK_FACTORY &check_fack)
	{
		const BaseLib::Options* op_cpl = option->getSubGroup("coupling");
	    if (op_cpl->hasSubGroup("M")) {
	        const BaseLib::Options* op_sub = op_cpl->getSubGroup("M");
	        T_M *sys = buildMonolithicSystem(op_sub, eqs_fac);
	        return sys;
	    } else if (op_cpl->hasSubGroup("P")) {
	        const BaseLib::Options* op_sub = op_cpl->getSubGroup("P");
	        T_P *sys = buildPartitionedSystem(op_sub, eqs_fac, check_fack);
	        return sys;
	    }
	    return 0;
	}

	std::vector<T_M*>& getListOfMonolithicSystem()
	{
		return _vec_m;
	}

	std::vector<std::string>& getListOfMonolithicSystemName()
	{
		return _vec_m_name;
	}

private:
	template <class T_EQS_FACTORY>
	T_M* buildMonolithicSystem(const BaseLib::Options *option, T_EQS_FACTORY &eqs_fac)
	{
		std::string eqs_name = option->getOption("name");
		T_M* eqs = eqs_fac.create(eqs_name);
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
		_vec_m.push_back(eqs);
		_vec_m_name.push_back(eqs_name);
		return eqs;
	}

	template <class T_EQS_FACTORY, class T_CHECK_FACTORY>
	T_P* buildPartitionedSystem(const BaseLib::Options *option, T_EQS_FACTORY &eqs_fac, T_CHECK_FACTORY check_fac)
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
		const BaseLib::Options* op_problems = option->getSubGroup("problems");
		for (BaseLib::Options::const_iterator itr=op_problems->begin(); itr!=op_problems->end(); ++itr) {
			std::string str = itr->first;
			BaseLib::Options* op_sub = static_cast<BaseLib::Options*>(itr->second);
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
