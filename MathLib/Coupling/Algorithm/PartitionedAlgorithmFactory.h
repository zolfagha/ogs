
#pragma once

#include <string>
#include "BlockGaussSeidelMethod.h"
#include "BlockJacobiMethod.h"
#include "IConvergenceCheck.h"

namespace MathLib
{

class PartitionedAlgorithmFactory
{
public:
//	template<class T>
	static IPartitionedAlgorithm* create(const std::string &name, IConvergenceCheck *checker, size_t max_itr, double epsilon)
	{
		if (name.compare("Jacobi")==0) {
			return new BlockJacobiMethod(*checker, epsilon, max_itr);
		} else if (name.compare("Gauss")==0) {
			return new BlockGaussSeidelMethod(*checker, epsilon, max_itr);
		}
		return 0;
	};
};


} //end
