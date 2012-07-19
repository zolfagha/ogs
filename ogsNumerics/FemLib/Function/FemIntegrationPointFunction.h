
#pragma once

#include <valarray>
#include <cmath>

#include "MathLib/Vector.h"
#include "FemLib/Core/DataType.h"
#include "DiscreteLib/Vector/DiscreteVector.h"


namespace FemLib
{

/**
 * \brief Template class for FEM integration point-based functions
 */
template<typename Tvalue>
class TemplateFEMIntegrationPointFunction : public NumLib::ITXFunction
{
public:
    typedef MathLib::TemplateVectorX<Tvalue> IntegrationPointVectorType;
    //typedef std::valarray<Tvalue> IntegrationPointVectorType;
    typedef DiscreteLib::DiscreteVector<IntegrationPointVectorType > DiscreteVectorType;

    TemplateFEMIntegrationPointFunction(DiscreteLib::DiscreteSystem &dis)
    {
        initialize(dis);
    };

    TemplateFEMIntegrationPointFunction(const TemplateFEMIntegrationPointFunction &src)
    {
        initialize(*src._discrete_system);
        (*this->_values) = (*src._values);
    };

    TemplateFEMIntegrationPointFunction<Tvalue>* clone() const
    {
        TemplateFEMIntegrationPointFunction<Tvalue> *obj = new TemplateFEMIntegrationPointFunction<Tvalue>(*this);
        return obj;
    };

    const MeshLib::IMesh* getMesh() const
    {
        return this->_discrete_system->getMesh();
    }

    DiscreteLib::DiscreteSystem* getDiscreteSystem() const {return _discrete_system;};

    void eval(const NumLib::TXPosition /*pt*/, NumLib::ITXFunction::DataType &v) const
    {
        v = (*_values)[0][0]; //TODO
    };

    void setIntegrationPointValue( size_t i_e, size_t ip, Tvalue &q )
    {
        assert(ip<(*_values)[i_e].size());
        (*_values)[i_e][ip] = q;
    }

    void setNumberOfIntegationPoints(size_t i_e, size_t n)
    {
        (*_values)[i_e].resize(n);
    }

    const IntegrationPointVectorType& getIntegrationPointValues(size_t i_e) const
    {
        return (*_values)[i_e];
    }

//    double norm_diff(const TemplateFEMIntegrationPointFunction<Tvalue> &ref) const
//    {
//        const size_t n = _values->size();
//        if (n!=ref._values->size()) {
//            std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is not same." << std::endl;
//            return .0;
//        }
//
//        double mnorm = .0;
//        for (size_t i=0; i<n; ++i) {
//            const IntegrationPointVectorType &val1 = (*_values)[i];
//            const IntegrationPointVectorType &val2 = (*ref._values)[i];
//            const size_t n_gp = val1.size();
//            if (n_gp!=val2.size()) {
//                std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is not same." << std::endl;
//                return .0;
//            } else if (n_gp==0) {
//                std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is zero." << std::endl;
//                return .0;
//            }
//            IntegrationPointVectorType val_diff = val1 - val2;
//
////            val_diff = std::abs(val_diff);
//            double val_diff_max = .0; // val_diff.max
//            for (size_t j=0; j<val_diff.size(); j++) {
//                val_diff_max = std::max(val_diff_max, val_diff[j].array().abs().maxCoeff());
//            }
//            mnorm = std::max(mnorm, val_diff_max);
//        }
//
//        return mnorm;
//    }

    const DiscreteVectorType* getNodalValues() const
    {
        return _values;
    }

    void printout() const
    {
        std::cout << "integration_pt_values = ";
        for (size_t i=_values->getRangeBegin(); i<_values->getRangeEnd(); ++i) {
            const IntegrationPointVectorType &val1 = (*_values)[i];
            std::cout << "(";
            for (size_t j=0; j<val1.size(); ++j) std::cout << val1[j] << " ";
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

private:
    DiscreteLib::DiscreteSystem* _discrete_system;
    DiscreteVectorType* _values;

    void initialize(DiscreteLib::DiscreteSystem &dis)
    {
        _discrete_system = &dis;
        _values = _discrete_system->createVector<DiscreteVectorType>(dis.getMesh()->getNumberOfElements());
    }
};

//typedef TemplateFEMIntegrationPointFunction<double> FEMIntegrationPointFunctionScalar;
//typedef TemplateFEMIntegrationPointFunction<MathLib::Vector2D> FEMIntegrationPointFunctionVector2d;
typedef TemplateFEMIntegrationPointFunction<LocalVector> FEMIntegrationPointFunctionVector;

} //end
