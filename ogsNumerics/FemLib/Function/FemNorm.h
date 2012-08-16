/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemNorm.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "DiscreteLib/Serial/DiscreteVector.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "FemLib/Function/FemIntegrationPointFunction.h"


namespace FemLib
{

template <typename T_DIS_SYS, typename Tvalue>
class NormOfFemNodalFunction
{
public:
    typedef typename FEMIntegrationPointFunctionVector<T_DIS_SYS>::type MyIntegrationPointFunctionVector;

    explicit NormOfFemNodalFunction(T_DIS_SYS* dis) : _dis(dis), _v_diff(0)
    {
    }

    ~NormOfFemNodalFunction()
    {
        if (_v_diff!=0)
            _dis->deleteVector(_v_diff);
    }

    ///
    double operator()(const TemplateFEMNodalFunction<T_DIS_SYS, Tvalue> &ref1, const TemplateFEMNodalFunction<T_DIS_SYS, Tvalue> &ref2)
    {
        const DiscreteLib::IDiscreteVector<Tvalue>* vec1 = ref1.getNodalValues();
        const DiscreteLib::IDiscreteVector<Tvalue>* vec2 = ref2.getNodalValues();
        const size_t n = vec1->size();
        if (n!=vec2->size()) {
            std::cout << "***Warning in NormOfFemNodalFunction::norm(): size of two vectors is not same." << std::endl;
            return .0;
        }

        if (_v_diff==0) {
            _v_diff = _dis->template createVector<Tvalue >(n);
        } else if (_v_diff->size() != n) {
            _dis->deleteVector(_v_diff);
            _v_diff = _dis->template createVector<Tvalue >(n);
        }

        *_v_diff = *vec1;
        *_v_diff -= *vec2;

        return MathLib::norm_max(*_v_diff, n);
    }

    double operator()(const MyIntegrationPointFunctionVector &ref1, const MyIntegrationPointFunctionVector &ref2) const
    {
        const typename MyIntegrationPointFunctionVector::MyDiscreteVector* vec1 = ref1.getElementValues();
        const typename MyIntegrationPointFunctionVector::MyDiscreteVector* vec2 = ref2.getElementValues();
        const size_t n = vec1->size();
        if (n!=vec2->size()) {
            std::cout << "***Warning in NormOfFemNodalFunction::norm(): size of two vectors is not same." << std::endl;
            return .0;
        }

        double mnorm = .0;
        for (size_t i=0; i<n; ++i) {
            const typename MyIntegrationPointFunctionVector::IntegrationPointVectorType &val1 = (*vec1)[i];
            const typename MyIntegrationPointFunctionVector::IntegrationPointVectorType &val2 = (*vec2)[i];
            const size_t n_gp = val1.size();
            if (n_gp!=val2.size()) {
                std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is not same." << std::endl;
                return .0;
            } else if (n_gp==0) {
                std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is zero." << std::endl;
                return .0;
            }
            typename MyIntegrationPointFunctionVector::IntegrationPointVectorType val_diff = val1 - val2;

//            val_diff = std::abs(val_diff);
            double val_diff_max = .0; // val_diff.max
            for (size_t j=0; j<val_diff.size(); j++) {
                val_diff_max = std::max(val_diff_max, val_diff[j].array().abs().maxCoeff());
            }
            mnorm = std::max(mnorm, val_diff_max);
        }

        return mnorm;
    }

private:
    T_DIS_SYS* _dis;
    DiscreteLib::IDiscreteVector<Tvalue>* _v_diff;
};


} //end
