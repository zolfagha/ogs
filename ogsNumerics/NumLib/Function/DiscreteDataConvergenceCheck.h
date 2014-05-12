/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteDataConvergenceCheck.h
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>

#include "MathLib/Vector.h"
#include "NumLib/Function/NormOfDiscreteDataFunction.h"
#include "NumLib/Coupling/Algorithm/IConvergenceCheck.h"
#include "NumLib/Function/ITXDiscreteFunction.h"

namespace NumLib
{

struct DiscreteDataType
{
	enum type
	{
		NodalScalar,
		NodalVector,
		IntegrationPointScalar,
		IntegrationPointVector
	};
};

/**
 * \brief Convergence check for iterative calculation of discrete data
 *
 */
class DiscreteDataConvergenceCheck : public NumLib::IConvergenceCheck
{
public:
    ///
    DiscreteDataConvergenceCheck() { };

    ///
    virtual ~DiscreteDataConvergenceCheck() {};

    ///
    virtual bool isConverged(const std::vector<unsigned> &vec_var_id, NumLib::UnnamedParameterSet& vars_prev, NumLib::UnnamedParameterSet& vars_current, double eps, double &v_diff)
    {
        if (_vec_value_type.size()!=vec_var_id.size())
        	return true;
        for (size_t ii=0; ii<vec_var_id.size(); ii++) {
        	size_t i = vec_var_id[ii];
            v_diff = .0;
			if (_vec_value_type[ii]==DiscreteDataType::NodalScalar)
				v_diff = calculateDifference<double>(vars_prev, vars_current, i);
			else if (_vec_value_type[ii]==DiscreteDataType::NodalVector)
				v_diff = calculateDifference<MathLib::LocalVector>(vars_prev, vars_current, i);
			else if (_vec_value_type[ii]==DiscreteDataType::IntegrationPointScalar)
				v_diff = calculateDifference<MathLib::TemplateVectorX<double> >(vars_prev, vars_current, i);
			else if (_vec_value_type[ii]==DiscreteDataType::IntegrationPointVector)
				v_diff = calculateDifference<MathLib::TemplateVectorX<MathLib::LocalVector> >(vars_prev, vars_current, i);

            if (v_diff>eps) {
                return false;
            }
        }
        return true;
    }

    ///
    template <class T>
    double calculateDifference(NumLib::UnnamedParameterSet &vars_prev, NumLib::UnnamedParameterSet &vars_current, size_t i)
    {
        const NumLib::ITXDiscreteFunction<T>* f_fem_prev = vars_prev.get<NumLib::ITXDiscreteFunction<T> >(i);
        const NumLib::ITXDiscreteFunction<T>* f_fem_cur = vars_current.get<NumLib::ITXDiscreteFunction<T> >(i);
        NumLib::NormOfDiscreteDataFunction<T> _norm;
        return _norm(*f_fem_prev, *f_fem_cur);
    }

    void addValueType(DiscreteDataType::type t) { _vec_value_type.push_back(t);}

private:
    std::vector<DiscreteDataType::type> _vec_value_type;
};

}
