
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
    virtual void extrapolate(TemplateFEMIntegrationPointFunction<Tvalue> &ele, TemplateFEMNodalFunction<Tvalue> &nod) = 0;
};

/**
 * \brief Extrapolate sampling point values by taking average
 */
template<typename Tvalue>
class FEMExtrapolationAverage : public IFemExtrapolation<typename Tvalue>
{
public:
    void extrapolate(TemplateFEMIntegrationPointFunction<Tvalue> &ele_var, TemplateFEMNodalFunction<Tvalue> &nod_var)
    {
        const MeshLib::IMesh* msh = ele_var.getMesh();
        MeshLib::TopologyNode2Elements node2eles(&msh);
        std::vector<double> vec_v(msh->getNumberOfNodes());

        for (size_t i=0; i<msh->getNumberOfElements(); i++) {
            const MeshLib::IElement* e = msh->getElemenet(i);
            const std::vector<Tvalue> &gp_values = ele_var.getIntegrationPointValues(i);
            Tvalue ele_avg = .0;
            for (size_t j=0; j<gp_values.size(); j++)
                ele_avg += gp_values[j];
            ele_avg /= gp_values.size();

            for (size_t j=0; j<e->getNumberOfNodes(); j++) {
                size_t n_conn_eles = node2eles.getConnectedElements(e->getNodeID(j)).size();
                vec_v[e->getNodeID(j)] += ele_avg / static_cast<double>(n_conn_eles);
            }
        }

        nod_var.setNodalValues(&vec_v[0]);
    }
};

/**
 * \brief Linear extrapolation method
 */
template<typename Tvalue>
class FEMExtrapolationLinear : public IFemExtrapolation<typename Tvalue>
{
private:
    double calcXi_p(MeshLib::IElement* e, FemIntegrationGaussBase* gauss)
    {
        double Xi_p = 0.0;
        if (e->getElementType() == MeshLib::ElementType::QUAD || e->getElementType() == MeshLib::ElementType::HEXAHEDRON) {
            double r = .0;
            const size_t nGauss = gauss->getSamplingLevel();
            for (size_t gp=0; gp<nGauss; gp++) {
                r = MathLib::GaussLegendre::getPoint(nGauss, gp)
                if (fabs(r) > Xi_p)
                    Xi_p = fabs(r);
            }
            r = 1.0 / Xi_p;
            Xi_p = r;
        }

        return Xi_p;
    }

    int getLocalIndex(MeshLib::IElement* e, FemIntegrationGaussBase* gauss, size_t igp)
    {
        int LoIndex = -1;
        double r,s,t;
        size_t nGauss = gauss->getSamplingLevel();
        size_t gp_r = igp/nGauss;
        size_t gp_s = igp%nGauss;
        size_t gp_t = 0;
        const double MKleinsteZahl = std::numeric_limits<double>::epsilon();

        switch (e->getElementType())
        {
        case MeshLib::ElementType::QUAD:
            r = MathLib::GaussLegendre::getPoint(nGauss, gp_r);
            s = MathLib::GaussLegendre::getPoint(nGauss, gp_s);
            if (r > 0.0 && s > 0.0)
                LoIndex = 0;
            else if (r < 0.0 && s > 0.0)
                LoIndex = 1;
            else if (r < 0.0 && s < 0.0)
                LoIndex = 2;
            else if (r > 0.0 && s < 0.0)
                LoIndex = 3;
            else if (fabs(r) < MKleinsteZahl && s > 0.0)
                LoIndex = 4;
            else if (r < 0.0 && fabs(s) < MKleinsteZahl)
                LoIndex = 5;
            else if (fabs(r) < MKleinsteZahl && s < 0.0)
                LoIndex = 6;
            else if (r > 0.0 && fabs(s) < MKleinsteZahl)
                LoIndex = 7;
            else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
                LoIndex = 8;
            break;
        case MeshLib::ElementType::HEXAHEDRON:
            r = MathLib::GaussLegendre::getPoint(nGauss, gp_r);
            s = MathLib::GaussLegendre::getPoint(nGauss, gp_s);
            t = MathLib::GaussLegendre::getPoint(nGauss, gp_t);

            if (t > 0.0)
            {
                if (r > 0.0 && s > 0.0)
                    LoIndex = 0;
                else if (r < 0.0 && s > 0.0)
                    LoIndex = 1;
                else if (r < 0.0 && s < 0.0)
                    LoIndex = 2;
                else if (r > 0.0 && s < 0.0)
                    LoIndex = 3;
                else if (fabs(r) < MKleinsteZahl && s > 0.0)
                    LoIndex = 8;
                else if (r < 0.0 && fabs(s) < MKleinsteZahl)
                    LoIndex = 9;
                else if (fabs(r) < MKleinsteZahl && s < 0.0)
                    LoIndex = 10;
                else if (r > 0.0 && fabs(s) < MKleinsteZahl)
                    LoIndex = 11;
                else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
                    return -1;
            }
            else if (fabs(t) < MKleinsteZahl)
            {
                if (fabs(r) < MKleinsteZahl || fabs(s) < MKleinsteZahl)
                    return -1;
                if (r > 0.0 && s > 0.0)
                    LoIndex = 16;
                else if (r < 0.0 && s > 0.0)
                    LoIndex = 17;
                else if (r < 0.0 && s < 0.0)
                    LoIndex = 18;
                else if (r > 0.0 && s < 0.0)
                    LoIndex = 19;
            }
            if (t < 0.0)
            {
                if (r > 0.0 && s > 0.0)
                    LoIndex = 4;
                else if (r < 0.0 && s > 0.0)
                    LoIndex = 5;
                else if (r < 0.0 && s < 0.0)
                    LoIndex = 6;
                else if (r > 0.0 && s < 0.0)
                    LoIndex = 7;
                else if (fabs(r) < MKleinsteZahl && s > 0.0)
                    LoIndex = 12;
                else if (r < 0.0 && fabs(s) < MKleinsteZahl)
                    LoIndex = 13;
                else if (fabs(r) < MKleinsteZahl && s < 0.0)
                    LoIndex = 14;
                else if (r > 0.0 && fabs(s) < MKleinsteZahl)
                    LoIndex = 15;
                else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
                    return -1;
            }
            break;
        default:
            std::cerr << "CElement::GetLocalIndex invalid mesh element type given"
                << std::endl;
        }
        return LoIndex;
    }

    void getExtropoGaussPoints(const MeshLib::IElement *e, const int i, double Xi_p, double* u)
    {
        int j = 0;
        MeshLib::ElementType::type ElementType = e->getElementType();
        //
        switch (ElementType) {
        case MeshLib::ElementType::LINE:
            break;
        case MeshLib::ElementType::TRIANGLE:
            switch (i) {
            case 0:
                unit[0] = -0.1666666666667;
                unit[1] = -0.1666666666667;
                break;
            case 1:
                unit[0] = 1.6666666666667;
                unit[1] = -0.1666666666667;
                break;
            case 2:
                unit[0] = -0.1666666666667;
                unit[1] = 1.6666666666667;
                break;
            }
            break;
        case MeshLib::ElementType::QUAD:
            switch (i) {
            case 0:
                unit[0] = Xi_p;
                unit[1] = Xi_p;
                break;
            case 1:
                unit[0] = -Xi_p;
                unit[1] = Xi_p;
                break;
            case 2:
                unit[0] = -Xi_p;
                unit[1] = -Xi_p;
                break;
            case 3:
                unit[0] = Xi_p;
                unit[1] = -Xi_p;
                break;
            }
            break;
        case MeshLib::ElementType::HEXAHEDRON:
            if (i < 4) {
                j = i;
                unit[2] = Xi_p;
            } else  {
                j = i - 4;
                unit[2] = -Xi_p;
            }
            switch (j) {
            case 0:
                unit[0] = Xi_p;
                unit[1] = Xi_p;
                break;
            case 1:
                unit[0] = -Xi_p;
                unit[1] = Xi_p;
                break;
            case 2:
                unit[0] = -Xi_p;
                unit[1] = -Xi_p;
                break;
            case 3:
                unit[0] = Xi_p;
                unit[1] = -Xi_p;
                break;
            }
            break;
        case MeshLib::ElementType::TETRAHEDRON:
            switch (i) {
            case 0:
                unit[0] = -0.166666666666667;
                unit[1] = -0.166666666666667;
                unit[2] = -0.166666666666667;
                break;
            case 1:
                unit[0] = 1.5;
                unit[1] = -0.166666666666667;
                unit[2] = -0.166666666666667;
                break;
            case 2:
                unit[0] = -0.166666666666667;
                unit[1] = 1.5;
                unit[2] = -0.166666666666667;
                break;
            case 3:
                unit[0] = -0.166666666666667;
                unit[1] = -0.166666666666667;
                unit[2] = 1.5;
            }
            break;
        default:
            unit[0] = unit[1] = unit[2] = 0.;
            break;
        }
    }

public:
    void extrapolate(TemplateFEMIntegrationPointFunction<Tvalue> &ele_var, TemplateFEMNodalFunction<Tvalue> &nod_var)
    {
        MeshLib::IMesh* msh = ele_var.getMesh();
        MeshLib::TopologyNode2Elements node2eles(msh);
        std::vector<double> vec_v(msh->getNumberOfNodes());

        for (size_t i=0; i<msh->getNumberOfElements(); i++) {
            const MeshLib::IElement* e = msh->getElemenet(i);
            const std::vector<Tvalue> &gp_values = ele_var.getIntegrationPointValues(i);

            IFiniteElement *fe = nod_var.getFiniteElement(e);
            FemIntegrationGaussBase* gauss = static_cast<FemIntegrationGaussBase*>(fe->getIntegrationMethod());

            std::vector<Tvalue> reordered_gp_values(gp_values.size());
            for (size_t j=0; j<gp_values.size(); j++) {
                size_t nod_id = j;
                if (e->getElementType()==MeshLib::ElementType::QUAD || e->getElementType()==MeshLib::ElementType::HEAD) {
                    nod_id = getLocalIndex(e, j);
                }
                reordered_gp_values[nod_id] = gp_values[j];
            }

            double Xi_p = .0;
            if (e->getElementType()==MeshLib::ElementType::QUAD || e->getElementType()==MeshLib::ElementType::HEAD) {
                Xi_p = calcXi_p(e, gauss);
            }

            size_t i_s = 0;
            size_t i_e = e->getNumberOfNodes();
            size_t ish = 0;
            if (e->getElementType()==MeshLib::ElementType::TETRAHEDRON) {
                i_s = 1;
                i_e = e->getNumberOfNodes() + 1;
                ish = 1;
            }

            double x[3];
            for (size_t j=0; j<e->getNumberOfNodes(); j++) {
                getExtropoGaussPoints(e, i, Xi_p, x);
                fe->computeBasisFunctions(x);
                Matrix<double> *N = fe->getBasisFunction();
                //ComputeShapefct(1); // Linear interpolation function
                double EV = .0;
                for(size_t k=i_s; k<i_e; k++)
                    EV += reordered_gp_values[k] * (*N)(0,k - ish);

                size_t n_conn_eles = node2eles.getConnectedElements(e->getNodeID(j)).size();
                vec_v[e->getNodeID(j)] += EV / static_cast<double>(n_conn_eles);
            }
        }

        nod_var.setNodalValues(&vec_v[0]);
    }
};


struct FEMExtrapolationMethod
{
    enum type {
        Linear = 1,
        Average = 2,
        LeastSquare = 3,
        INVALID = -1
    };
};

template<typename Tvalue>
class FEMExtrapolationFactory
{
public:
    static IFemExtrapolation<Tvalue>* create(FEMExtrapolationMethod::type tp)
    {
        switch (tp) {
            case FEMExtrapolationMethod::Linear:
                return new FEMExtrapolationLinear<Tvalue>();
            case FEMExtrapolationMethod::Average:
                return new FEMExtrapolationAverage<Tvalue>();
        }
        return 0;
    };
};


}
