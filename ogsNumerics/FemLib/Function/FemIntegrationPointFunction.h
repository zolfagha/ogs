
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
    typedef DiscreteLib::DiscreteVector<IntegrationPointVectorType > DiscreteVectorType;

    ///
    explicit TemplateFEMIntegrationPointFunction(DiscreteLib::DiscreteSystem &dis)
    {
        initialize(dis);
    };

    ///
    explicit TemplateFEMIntegrationPointFunction(const TemplateFEMIntegrationPointFunction &src)
    {
        initialize(*src._discrete_system);
        (*this->_values) = (*src._values);
    };

    ///
    TemplateFEMIntegrationPointFunction<Tvalue>* clone() const
    {
        TemplateFEMIntegrationPointFunction<Tvalue> *obj = new TemplateFEMIntegrationPointFunction<Tvalue>(*this);
        return obj;
    };

    ///
    const MeshLib::IMesh* getMesh() const
    {
        return this->_discrete_system->getMesh();
    }

    ///
    DiscreteLib::DiscreteSystem* getDiscreteSystem() const {return _discrete_system;};

    ///
    virtual void eval(const NumLib::TXPosition x, NumLib::ITXFunction::DataType &v) const
    {
        switch (x.getIdObjectType()) {
        case NumLib::TXPosition::IntegrationPoint:
            {
                size_t ele_id = x.getId(0);
                size_t gp_id = x.getId(1);
                IntegrationPointVectorType &gp_values = (*_values)[ele_id];
                if (gp_values.size()==0) return;
                v = gp_values[gp_id]; 
            }
            break;
        case NumLib::TXPosition::Element:
            {
                // calculate mean value
                size_t ele_id = x.getId();
                IntegrationPointVectorType &gp_values = (*_values)[ele_id];
                if (gp_values.size()==0) return;
                Tvalue val = gp_values[0]; 
                val *= .0;
                for (size_t i=0; i<gp_values.size(); i++)
                    val += gp_values[i];
                val /= gp_values.size();
                v = val;
            }
            break;
        default:
            break;
        }
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

template <> 
void TemplateFEMIntegrationPointFunction<double>::eval(const NumLib::TXPosition x,  NumLib::ITXFunction::DataType &v) const
{
    size_t ele_id = x.getId();
    IntegrationPointVectorType &gp_values = (*_values)[ele_id];
    if (gp_values.size()==0) return;
    double val = gp_values[0]; 
    val *= .0;
    for (size_t i=0; i<gp_values.size(); i++)
        val += gp_values[i];
    val /= gp_values.size();

    v.resize(1,1);
    v(0,0) = val;
    NumLib::ITXFunction::DataType mat(1,1);
};

typedef TemplateFEMIntegrationPointFunction<LocalVector> FEMIntegrationPointFunctionVector;

} //end
