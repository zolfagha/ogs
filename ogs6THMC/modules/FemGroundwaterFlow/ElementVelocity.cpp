
#include "ElementVelocity.h"

#include "MathLib/Vector.h"

namespace Geo
{

void FunctionElementVelocity::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();

    size_t msh_id = option.getOption<size_t>("MeshID");
    DiscreteLib::DiscreteSystem* dis = femData->list_dis_sys[msh_id];
    MeshLib::IMesh* msh = dis->getMesh();

    _dis = dis;
    _vel = new FemLib::FEMIntegrationPointFunctionVector(*dis);
    //this->setOutput(Velocity, _vel);
}


int FunctionElementVelocity::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    const MeshLib::IMesh *msh = _dis->getMesh();
    FemLib::FemNodalFunctionScalar *head = (FemLib::FemNodalFunctionScalar*)getInput(Head);
    FemLib::FEMIntegrationPointFunctionVector *vel = _vel;;


    FemLib::LagrangianFeObjectContainer* feObjects = head->getFeObjectContainer();
    //calculate vel (vel=f(h))
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        MeshLib::IElement* e = msh->getElemenet(i_e);
        size_t mat_id = e->getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

        FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
        NumLib::LocalVector local_h(e->getNumberOfNodes());
        for (size_t j=0; j<e->getNumberOfNodes(); j++)
            local_h[j] = head->getValue(e->getNodeID(j));
        // for each integration points
        FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
        double r[2] = {};
        const size_t n_gp = integral->getNumberOfSamplingPoints();
        vel->setNumberOfIntegationPoints(i_e, n_gp);
        NumLib::LocalVector xi(e->getNumberOfNodes());
        NumLib::LocalVector yi(e->getNumberOfNodes());
        for (size_t i=0; i<e->getNumberOfNodes(); i++) {
            const GeoLib::Point* pt = msh->getNodeCoordinatesRef(e->getNodeID(i));
            xi[i] = (*pt)[0];
            yi[i] = (*pt)[1];
        }
        for (size_t ip=0; ip<n_gp; ip++) {
            NumLib::LocalVector q(2);
            q[0] = .0;
            q[1] = .0;
            integral->getSamplingPoint(ip, r);
            fe->computeBasisFunctions(r);
            const NumLib::LocalMatrix* dN = fe->getGradBasisFunction();
            NumLib::LocalMatrix* N = fe->getBasisFunction();
            std::vector<double> xx(3, .0);
            NumLib::LocalVector tmp_v;
            tmp_v = (*N) * xi;
            xx[0] = tmp_v[0];
            tmp_v = (*N) * yi;
            xx[1] = tmp_v[0];
            NumLib::TXPosition pos(&xx[0]);

            NumLib::LocalMatrix k;
            pm->hydraulic_conductivity->eval(pos, k);
            if (k.rows()==1) {
                q.noalias() = (*dN) * local_h * (-1.0) * k(0,0);
            } else {
                q.noalias() = (*dN) * k * local_h * (-1.0);
            }
            vel->setIntegrationPointValue(i_e, ip, q);
        }
    }
    setOutput(Velocity, vel);
    return 0;
}
}

