
#pragma once

#include <cmath>
#include <limits>

#include "MeshLib/Topology/Topology.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/Core/Integration.h"

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

template<typename Tvalue>
class FemExtrapolationAverage : public IFemExtrapolation<Tvalue>
{
public:
    void extrapolate(TemplateFEMIntegrationPointFunction<Tvalue> &ele_var, TemplateFEMNodalFunction<Tvalue> &nod_var)
    {
        const MeshLib::IMesh* msh = ele_var.getMesh();
        LagrangianFeObjectContainer* feObjects = nod_var.getFeObjectContainer();
        MeshLib::TopologySequentialNodes2Elements node2eles(*msh);
        std::vector<Tvalue> vec_v(msh->getNumberOfNodes());

        for (size_t i=0; i<msh->getNumberOfElements(); i++) {
            MeshLib::IElement* e = msh->getElemenet(i);
            const std::valarray<Tvalue> &gp_values = ele_var.getIntegrationPointValues(i);
            std::vector<Tvalue> vec_gp_values(&gp_values[0], &gp_values[0]+gp_values.size());
            const size_t e_nnodes = e->getNumberOfNodes();
            std::vector<Tvalue> nodal_values(e_nnodes);
            IFiniteElement *fe = feObjects->getFeObject(*e);
            fe->extrapolate(vec_gp_values, nodal_values);

            for (size_t j=0; j<e_nnodes; j++) {
                const size_t nod_id = e->getNodeID(j);
                size_t n_conn_eles = node2eles.getConnectedElements(nod_id).size();
                vec_v[nod_id] += nodal_values[j] / static_cast<double>(n_conn_eles);
            }
        }

        nod_var.setNodalValues(&vec_v[0], 0, vec_v.size());
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
