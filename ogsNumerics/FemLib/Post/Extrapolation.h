/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Extrapolation.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cmath>
#include <limits>

#include "MeshLib/Topology/Topology.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/Core/Integration/Integration.h"

namespace FemLib
{
/**
 * \brief Extrapolation methods for integration point values
 */
template<typename Tvalue>
class IFemExtrapolation
{
public:
    virtual ~IFemExtrapolation() {};
    virtual void extrapolate(TemplateFEMIntegrationPointFunction<Tvalue> &ele, TemplateFEMNodalFunction<Tvalue> &nod) = 0;
};

/**
 * \brief
 */
template<typename Tvalue>
class FemExtrapolationAverage : public IFemExtrapolation<Tvalue>
{
public:

    void extrapolate(TemplateFEMIntegrationPointFunction<Tvalue> &ele_var, TemplateFEMNodalFunction<Tvalue> &nod_var)
    {
        const MeshLib::IMesh* msh = ele_var.getMesh();
        LagrangianFeObjectContainer* feObjects = nod_var.getFeObjectContainer();
        MeshLib::TopologySequentialNodes2Elements node2eles(*msh);
        Tvalue v0 = ele_var.getIntegrationPointValues(0)[0];
        v0 *= .0;
        DiscreteLib::IDiscreteVector<Tvalue>* node_vec = nod_var.getNodalValues();
        (*node_vec) = v0;

        for (size_t i=0; i<msh->getNumberOfElements(); i++) {
            MeshLib::IElement* e = msh->getElemenet(i);
            const typename TemplateFEMIntegrationPointFunction<Tvalue>::IntegrationPointVectorType &gp_values = ele_var.getIntegrationPointValues(i);
            std::vector<Tvalue> vec_gp_values(&gp_values[0], &gp_values[0]+gp_values.size());
            const size_t e_nnodes = e->getNumberOfNodes();
            std::vector<Tvalue> nodal_values(e_nnodes);
            IFiniteElement *fe = feObjects->getFeObject(*e);
            fe->extrapolate(vec_gp_values, nodal_values);

            for (size_t j=0; j<e_nnodes; j++) {
                const size_t nod_id = e->getNodeID(j);
                size_t n_conn_eles = node2eles.getConnectedElements(nod_id).size();
                (*node_vec)[nod_id] += nodal_values[j] / static_cast<double>(n_conn_eles);
            }
        }
    }
};


struct FEMExtrapolationMethod
{
    enum type {
        Average,
        INVALID
    };
};

template<typename Tvalue>
class FEMExtrapolationFactory
{
public:
    static IFemExtrapolation<Tvalue>* create(FEMExtrapolationMethod::type tp)
    {
        switch (tp) {
            case FEMExtrapolationMethod::Average:
                return new FemExtrapolationAverage<Tvalue>();
        }
        return 0;
    };
};


}
