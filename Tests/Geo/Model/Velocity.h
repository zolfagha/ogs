
#pragma once

#include "Base/CodingTools.h"
#include "MathLib/Function/Function.h"
#include "MathLib/Vector.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TimeStepping/TimeStep.h"

#include "Tests/Geo/Material/PorousMedia.h"

using namespace NumLib;

namespace Geo
{

class FunctionVelocity
	: public NumLib::TemplateTransientMonolithicSystem
{
    enum In { Head=0 };
    enum Out { Velocity=0 };
public:

    FunctionVelocity() 
    {
        TemplateTransientMonolithicSystem::resizeInputParameter(1);
        TemplateTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteLib::DiscreteSystem &dis, PorousMedia &pm)
    {
    	_dis = &dis;
    	_K = pm.hydraulic_conductivity;
        _vel = new FemLib::FEMIntegrationPointFunctionVector2d(dis, *dis.getMesh());
        //this->setOutput(Velocity, _vel);
    }

    int solveTimeStep(const TimeStep &time)
    {
        const MeshLib::IMesh *msh = _dis->getMesh();
        FemLib::FemNodalFunctionScalar *head = (FemLib::FemNodalFunctionScalar*)getInput(Head);
        FemLib::FEMIntegrationPointFunctionVector2d *vel = _vel;;

        FemLib::LagrangianFeObjectContainer* feObjects = head->getFeObjectContainer();
        //calculate vel (vel=f(h))
        for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
            MeshLib::IElement* e = msh->getElemenet(i_e);
            FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
            std::vector<double> local_h(e->getNumberOfNodes());
            for (size_t j=0; j<e->getNumberOfNodes(); j++)
                local_h[j] = head->getValue(e->getNodeID(j));
            // for each integration points
            FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
            double r[2] = {};
            const size_t n_gp = integral->getNumberOfSamplingPoints();
            vel->setNumberOfIntegationPoints(i_e, n_gp);
            std::vector<double> xi(e->getNumberOfNodes());
            std::vector<double> yi(e->getNumberOfNodes());
            for (size_t i=0; i<e->getNumberOfNodes(); i++) {
                const GeoLib::Point* pt = msh->getNodeCoordinatesRef(e->getNodeID(i));
                xi[i] = (*pt)[0];
                yi[i] = (*pt)[1];
            }
            for (size_t ip=0; ip<n_gp; ip++) {
                MathLib::Vector q(2);
                q[0] = .0;
                q[1] = .0;
                integral->getSamplingPoint(ip, r);
                fe->computeBasisFunctions(r);
                const MathLib::Matrix<double> *dN = fe->getGradBasisFunction();
                MathLib::Matrix<double>*N = fe->getBasisFunction();
                std::vector<double> xx(2);
                N->axpy(1.0, &xi[0], .0, &xx[0]);
                N->axpy(1.0, &yi[0], .0, &xx[1]);

                double k;
                _K->eval(&xx[0], k);
                dN->axpy(-k, &local_h[0], .0, &q[0]); // q = - K * dN * local_h;
                vel->setIntegrationPointValue(i_e, ip, q);
            }
        }
        setOutput(Velocity, vel);
        return 0;
    }

    double suggestNext(const TimeStep &time_current)
    {
        return .0;
    }

    bool isAwake(const TimeStep &time)
    {
        return true;
    }

    void accept(const TimeStep &time)
    {
        //std::cout << "Velocity=" << std::endl;
        //_vel->printout();
    };

private:
    DiscreteLib::DiscreteSystem* _dis;
    FemLib::FEMIntegrationPointFunctionVector2d* _vel;
    MathLib::SpatialFunctionScalar *_K;

    DISALLOW_COPY_AND_ASSIGN(FunctionVelocity);
};

}
