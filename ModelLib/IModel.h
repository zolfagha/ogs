
#pragma once

#include <vector>

#include "MathLib/LinAlg/Dense/Matrix.h"


namespace ModelLib
{

class IDiscretizationModel
{
public:
    //virtual void localAssembly(size_t i_ele, MathLib::Matrix<double> localA, double *localRHS) = 0;
    //virtual void getKnownXi(std::vector<size_t> vec_dof_if, std::vector<double> vec_values) = 0;
    //virtual void addRHS(std::vector<size_t> vec_dof_if, std::vector<double> vec_values) = 0;
};

}
