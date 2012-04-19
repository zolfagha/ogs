
#pragma once

#include <vector>
#include <algorithm>

#include "Base/CodingTools.h"

#include "MathLib/Vector.h"
#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/Core/DiscreteVector.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/FemElementObjectContainer.h"



namespace FemLib
{


/// Polynomial order
struct PolynomialOrder
{
    enum type {
        Linear = 1,
        Quadratic = 2,
        INVALID = -1
    };
};


/**
 * \brief Template class for FEM node-based functions (assuming Lagrangian elements)
 *
 * @tparam Tvalue Nodal value type, e.g. double, vector
 */
template<typename Tvalue>
class TemplateFEMNodalFunction : public MathLib::SpatialFunction<Tvalue>
{
public:
    /// @param msh 		Mesh
    /// @param order 	Polynomial order
    TemplateFEMNodalFunction(DiscreteLib::DiscreteSystem &dis, MeshLib::IMesh &msh, PolynomialOrder::type order, Tvalue v0)
    {
        initialize(dis, msh, order);
        resetNodalValues(v0);
    }

    /// @param msh Mesh
    /// @param order Polynomial order
    TemplateFEMNodalFunction(DiscreteLib::DiscreteSystem &dis, MeshLib::IMesh &msh, PolynomialOrder::type order)
    {
        initialize(dis, msh, order);
    }

    /// @param org source object for copying
    TemplateFEMNodalFunction(const TemplateFEMNodalFunction<Tvalue> &org)
    {
        assign(org);
    }

    ///
    virtual ~TemplateFEMNodalFunction()
    {
        Base::releaseObject(_feObjects);
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

    /// get the spatial dimension
    size_t getDimension() const { return _msh->getDimension(); }

    /// get the mesh
    const MeshLib::IMesh* getMesh() const { return _msh; }

    /// get the number of nodes
    size_t getNumberOfNodes() const { return _msh->getNumberOfNodes(); }


    /// evaluate this function at the given point
    void eval(const MathLib::SpatialPosition &pt, Tvalue &v)
    {
        throw "eval() is not implemented yet.";
        v = (*_nodal_values)[0];
    };

    /// get nodal value
    Tvalue& getValue(int node_id)
    {
        return (*_nodal_values)[node_id];
    }

//    /// get an array of nodal values
//    Tvalue* getNodalValues()
//    {
//        if (_nodal_values->size()>0)
//            return &(*_nodal_values)[0];
//        else
//            return 0;
//    }

    /// get nodal values
    DiscreteLib::DiscreteVector<Tvalue>* getNodalValues()
    {
        return _nodal_values;
    }

    const DiscreteLib::DiscreteVector<Tvalue>* getNodalValues() const
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
    void setNodalValues( const DiscreteLib::DiscreteVector<Tvalue> &x )
    {
    	*_nodal_values = x;
    }

    /// reset nodal values with the given value
    void resetNodalValues (Tvalue &v)
    {
        std::fill(_nodal_values->begin(), _nodal_values->end(), v);
    }

    LagrangianFeObjectContainer* getFeObjectContainer()
    {
        return _feObjects;
    }

    ///// get finite element object
    //IFiniteElement* getFiniteElement(MeshLib::IElement &e)
    //{
    //    _feObjects->setPolynomialOrder(_order);
    //    IFiniteElement* fe = _feObjects->getFeObject(e);
    //    fe->configure(e);
    //    fe->getIntegrationMethod()->initialize(e, 2);

    //    return fe;
    //}

    double norm_diff(const TemplateFEMNodalFunction<Tvalue> &ref) const
    {
    	const size_t n = _nodal_values->size();
    	if (n!=ref._nodal_values->size()) {
    		std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is not same." << std::endl;
    		return .0;
    	}

		const DiscreteLib::DiscreteVector<double>* vec_prev = ref.getNodalValues();
		const DiscreteLib::DiscreteVector<double>* vec_cur = this->getNodalValues();
		DiscreteLib::DiscreteVector<double> vec_diff(vec_prev->size());
		vec_diff = *vec_cur;
		vec_diff -= *vec_prev;
		return MathLib::norm_max(vec_diff, vec_diff.size());
    }

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
        _nodal_values = dis.createVector<Tvalue>(nnodes);
        _feObjects = new LagrangianFeObjectContainer(msh);
    }

    /// Assign this object from the given object
    void assign(const TemplateFEMNodalFunction<Tvalue> &org)
    {
        initialize(*org._discrete_system, *org._msh, org._order);
        std::copy(org._nodal_values->begin(), org._nodal_values->end(), _nodal_values->begin());
    }

private:
    DiscreteLib::DiscreteSystem* _discrete_system;
    DiscreteLib::DiscreteVector<Tvalue>* _nodal_values;
    MeshLib::IMesh* _msh;
    PolynomialOrder::type _order;
    LagrangianFeObjectContainer* _feObjects;
};

typedef TemplateFEMNodalFunction<double> FemNodalFunctionScalar;
typedef TemplateFEMNodalFunction<MathLib::Vector> FemNodalFunctionVector2d;

} //end
