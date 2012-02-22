
#pragma once

#include <vector>
#include <algorithm>

#include "Base/CodingTools.h"

#include "MathLib/Vector.h"
#include "MathLib/Function/Function.h"
#include "MeshLib/Core/IMesh.h"


#include "FemLib/Core/IFemElement.h"
#include "FemLib/FemElementObjectContainer.h"

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
class TemplateFEMNodalFunction : public MathLib::IFunction<Tvalue, GeoLib::Point>
{
public:
    /// @param msh Mesh
    /// @param order Polynomial order
    TemplateFEMNodalFunction(MeshLib::IMesh &msh, PolynomialOrder::type order, Tvalue v0) 
    {
        initialize(msh, order);
        resetNodalValues(v0);
    }

    /// @param msh Mesh
    /// @param order Polynomial order
    TemplateFEMNodalFunction(MeshLib::IMesh &msh, PolynomialOrder::type order) 
    {
        initialize(msh, order);
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
    MathLib::IFunction<Tvalue, GeoLib::Point>* clone() const 
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
    Tvalue eval(const GeoLib::Point &pt) 
    {
        throw std::exception("eval() is not implemented yet.");
        return _nodal_values[0];
    };

    /// get nodal value
    Tvalue& getValue(int node_id)
    {
        return _nodal_values[node_id];
    }

    /// get an array of nodal values
    Tvalue* getNodalValues() 
    {
        if (_nodal_values.size()>0)
            return &_nodal_values[0];
        else
            return 0;
    }

    std::vector<Tvalue>* getNodalValuesAsStdVec() 
    {
        return &_nodal_values;
    }

    /// set nodal values
    void setNodalValues( Tvalue* x ) 
    {
        std::copy(x, x+getNumberOfNodes(), _nodal_values.begin());
    }

    /// reset nodal values with the given value
    void resetNodalValues (Tvalue &v)
    {
        std::fill(_nodal_values.begin(), _nodal_values.end(), v);
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


private:
    std::vector<Tvalue> _nodal_values;
    MeshLib::IMesh* _msh;
    PolynomialOrder::type _order;
    LagrangianFeObjectContainer* _feObjects;

    /// initialize 
    void initialize(MeshLib::IMesh &msh, PolynomialOrder::type order)
    {
        _msh = &msh;
        _order = order;
        size_t nnodes = msh.getNumberOfNodes();
        _nodal_values.resize(nnodes);
        _feObjects = new LagrangianFeObjectContainer(msh);
    }

    /// assigne this object from the given object
    void assign(const TemplateFEMNodalFunction<Tvalue> &org)
    {
        initialize(*org._msh, org._order);
        std::copy(org._nodal_values.begin(), org._nodal_values.end(), _nodal_values.begin());
    }

};

typedef TemplateFEMNodalFunction<double> FemNodalFunctionScalar;
typedef TemplateFEMNodalFunction<MathLib::Vector2D> FemNodalFunctionVector2d;
typedef TemplateFEMNodalFunction<MathLib::Vector3D> FEMNodalFunctionVector3d;


/**
 * \brief Template class for FEM element-based functions
 */
template<typename Tvalue>
class TemplateFEMElementalFunction : public MathLib::IFunction<Tvalue,GeoLib::Point>
{
public:
    Tvalue eval(const GeoLib::Point &pt) {
        return _ele_values[0];
    };
    Tvalue& getValue(size_t id) const {
        return _ele_values[id];
    };
private:
    Tvalue* _ele_values;
};

/**
 * \brief Template class for FEM integration point-based functions
 */
template<typename Tvalue>
class TemplateFEMIntegrationPointFunction : public MathLib::IFunction<Tvalue,GeoLib::Point>
{
public:
    TemplateFEMIntegrationPointFunction(MeshLib::IMesh* msh) {
        _msh = msh;
        _values.resize(msh->getNumberOfElements());
    };

    MathLib::IFunction<Tvalue, GeoLib::Point>* clone() const {return 0;};

    const MeshLib::IMesh* getMesh() const {
        return _msh;
    }

    Tvalue eval(const GeoLib::Point &pt) {
        throw std::exception("The method or operation is not implemented.");
    };

    void setIntegrationPointValue( size_t i_e, size_t ip, Tvalue &q ) 
    {
        assert(ip<_values[i_e].size());
        _values[i_e][ip] = q;
    }

    void setNumberOfIntegationPoints(size_t i_e, size_t n) {
        _values[i_e].resize(n);
    }

    const std::vector<Tvalue>& getIntegrationPointValues(size_t i_e) const 
    {
        return _values[i_e];
    }

private:
    MeshLib::IMesh* _msh;
    std::vector<std::vector<Tvalue>> _values;

};

typedef TemplateFEMIntegrationPointFunction<double> FEMIntegrationPointFunctionScalar2d;
typedef TemplateFEMIntegrationPointFunction<MathLib::Vector2D> FEMIntegrationPointFunctionVector2d;

}



