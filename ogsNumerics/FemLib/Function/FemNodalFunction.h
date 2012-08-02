
#pragma once

#include <vector>
#include <algorithm>

#include "BaseLib/CodingTools.h"

#include "MathLib/Vector.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/Vector/DiscreteVector.h"

#include "NumLib/Function/TXFunction.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/FemElementObjectContainer.h"


#include "FemLib/Core/PolynomialOrder.h"


namespace FemLib
{

/**
 * \brief Template class for FEM node-based functions (assuming Lagrangian elements)
 *
 * This class represents the following
 * u^h(x) = N(x)*u_i
 *
 * @tparam Tvalue Nodal value type, e.g. double, vector
 */
template<typename Tvalue>
class TemplateFEMNodalFunction : public NumLib::ITXFunction
{
public:
    /// @param dis         Discrete system
    /// @param order     Polynomial order
    /// @param v0        initial value
    TemplateFEMNodalFunction(DiscreteLib::DiscreteSystem &dis, PolynomialOrder::type order, Tvalue v0)
    {
        initialize(dis, *dis.getMesh(), order);
        resetNodalValues(v0);
    }

    /// @param dis         Discrete system
    /// @param order Polynomial order
    TemplateFEMNodalFunction(DiscreteLib::DiscreteSystem &dis, PolynomialOrder::type order)
    {
        initialize(dis, *dis.getMesh(), order);
    }

    /// @param org source object for copying
    TemplateFEMNodalFunction(const TemplateFEMNodalFunction<Tvalue> &org)
    {
        assign(org);
    }

    ///
    virtual ~TemplateFEMNodalFunction()
    {
        BaseLib::releaseObject(_feObjects);
    }

    ///
    TemplateFEMNodalFunction &operator=(const TemplateFEMNodalFunction &org)
    {
        assign(org);
        return *this;
    }

    /// make a clone of this object
    /// @return MathLib::IFunction*
    TemplateFEMNodalFunction<Tvalue>* clone() const
    {
        TemplateFEMNodalFunction<Tvalue> *obj = new TemplateFEMNodalFunction<Tvalue>(*this);
        return obj;
    };

    /// get the mesh
    const MeshLib::IMesh* getMesh() const { return _msh; }

    DiscreteLib::DiscreteSystem* getDiscreteSystem() const {return _discrete_system;};

    ///
    size_t getNumberOfNodes() const {return _nodal_values->size();};


    /// evaluate this function at the given point
    //virtual void eval(const NumLib::TXPosition x, Tvalue &v) const
    //{
    //    v = (*_nodal_values)[x.getId()];
    //};

    virtual void eval(const NumLib::TXPosition x, NumLib::ITXFunction::DataType &v) const
    {
        switch (x.getIdObjectType()) {
        case NumLib::TXPosition::Node:
            {
                v = (*_nodal_values)[x.getId()];
            }
            break;
        default:
            break;
        }
    };

    /// get nodal value
    Tvalue& getValue(size_t node_id)
    {
        return (*_nodal_values)[node_id];
    }

    /// get nodal values
    DiscreteLib::IDiscreteVector<Tvalue>* getNodalValues()
    {
        return _nodal_values;
    }

    const DiscreteLib::IDiscreteVector<Tvalue>* getNodalValues() const
    {
        return _nodal_values;
    }

    /// set nodal values
    void setNodalValues( Tvalue* x, size_t i_start, size_t n )
    {
        for (size_t i=0; i<n; ++i)
            (*_nodal_values)[i+i_start] = x[i];
        //std::copy(x, x+getNumberOfNodes(), _nodal_values->begin());
    }

    /// set nodal values
    void setNodalValues( const DiscreteLib::IDiscreteVector<Tvalue> &x )
    {
        *_nodal_values = x;
    }

    /// reset nodal values with the given value
    void resetNodalValues (Tvalue &v)
    {
        //std::fill(_nodal_values->begin(), _nodal_values->end(), v);
        *_nodal_values = v;
    }

    /// get Finite element object container
    LagrangianFeObjectContainer* getFeObjectContainer()
    {
        return _feObjects;
    }

//    ///
//    double norm_diff(const TemplateFEMNodalFunction<Tvalue> &ref) const
//    {
//        const size_t n = _nodal_values->size();
//        if (n!=ref._nodal_values->size()) {
//            std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is not same." << std::endl;
//            return .0;
//        }
//
//        const DiscreteLib::DiscreteVector<double>* vec_prev = ref.getNodalValues();
//        const DiscreteLib::DiscreteVector<double>* vec_cur = this->getNodalValues();
//        DiscreteLib::DiscreteVector<double> vec_diff(vec_prev->size());
//        vec_diff = *vec_cur;
//        vec_diff -= *vec_prev;
//        return MathLib::norm_max(vec_diff, vec_diff.size());
//    }

    /// printout internal data for debugging
    void printout() const
    {
        std::cout << "nodal_values = ";
        for (size_t i=_nodal_values->getRangeBegin(); i<_nodal_values->getRangeEnd(); ++i)
            std::cout << (*_nodal_values)[i] << " ";
        std::cout << std::endl;
    }
private:
    /// initialize
    void initialize(DiscreteLib::DiscreteSystem &dis, MeshLib::IMesh &msh, PolynomialOrder::type order)
    {
        _discrete_system = &dis;
        _msh = &msh;
        _order = order;
        size_t nnodes = msh.getNumberOfNodes();
        _nodal_values = dis.createVector<DiscreteLib::DiscreteVector<Tvalue> >(nnodes);
        _feObjects = new LagrangianFeObjectContainer(msh);
    }

    /// Assign this object from the given object
    void assign(const TemplateFEMNodalFunction<Tvalue> &org)
    {
        initialize(*org._discrete_system, *org._msh, org._order);
        //std::copy(org._nodal_values->begin(), org._nodal_values->end(), _nodal_values->begin());
        for (size_t i=org._nodal_values->getRangeBegin(); i<org._nodal_values->getRangeEnd(); ++i)
            (*_nodal_values)[i] = (*org._nodal_values)[i];
    }

private:
    DiscreteLib::DiscreteSystem* _discrete_system;
    DiscreteLib::IDiscreteVector<Tvalue>* _nodal_values;
    MeshLib::IMesh* _msh;
    PolynomialOrder::type _order;
    LagrangianFeObjectContainer* _feObjects;
};

/// evaluate this function at the given point
template <> 
void TemplateFEMNodalFunction<double>::eval(const NumLib::TXPosition x,  NumLib::ITXFunction::DataType &v) const
{
    NumLib::ITXFunction::DataType val(1,1);
    val(0,0) = (*_nodal_values)[x.getId()];
    v = val;
};

typedef TemplateFEMNodalFunction<double> FemNodalFunctionScalar;
typedef TemplateFEMNodalFunction<NumLib::LocalVector> FemNodalFunctionVector2d;

} //end
