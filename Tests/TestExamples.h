
#pragma once

#include <vector>

#include "Base/BidirectionalMap.h"

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/Assembler/IElemenetWiseLinearEquationLocalAssembler.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"
#include "FemLib/Post/Extrapolation.h"


class GWFemTest
{
public:
    GeoLib::Rectangle *rec;
    DiscreteLib::DiscreteSystem *dis;
    MeshLib::IMesh *msh;
    MathLib::SpatialFunctionScalar *_K;
    FemLib::FemNodalFunctionScalar *head;
    FemLib::FEMIntegrationPointFunctionVector2d *vel;
    std::vector<FemLib::FemDirichletBC<double>*> vec_bc1;
    std::vector<FemLib::FemNeumannBC<double, double>*> vec_bc2;

    void define(MeshLib::IMesh *msh, MathLib::SpatialFunctionScalar *K=0)
    {
        //#Define a problem
        //geometry
        rec = new GeoLib::Rectangle(GeoLib::Point(0.0, 0.0, 0.0),  GeoLib::Point(2.0, 2.0, 0.0));
        GeoLib::Polyline* poly_left = rec->getLeft();
        GeoLib::Polyline* poly_right = rec->getRight();
        //mesh
        this->msh = msh;
        dis = new DiscreteLib::DiscreteSystem(*msh);
        //discretization
        head = new FemLib::FemNodalFunctionScalar(*dis, *msh, FemLib::PolynomialOrder::Linear);
        vel = new FemLib::FEMIntegrationPointFunctionVector2d(*dis, *msh);
        //bc
        vec_bc1.push_back(new FemLib::FemDirichletBC<double>(head, poly_right, false, new MathLib::SpatialFunctionConstant<double>(.0), new FemLib::DiagonalizeMethod()));
        vec_bc2.push_back(new FemLib::FemNeumannBC<double, double>(head, poly_left, false, new MathLib::SpatialFunctionConstant<double>(-1e-5)));
        // mat
        _K = (K!=0) ? K : new MathLib::SpatialFunctionConstant<double>(1.e-11);
    }

    static void calculateHead(GWFemTest &gw)
    {
        const MeshLib::IMesh *msh = gw.msh;
        for (size_t i=0; i<gw.vec_bc1.size(); i++) gw.vec_bc1[i]->setup();
        for (size_t i=0; i<gw.vec_bc2.size(); i++) gw.vec_bc2[i]->setup();
        // global EQS
        const size_t n_dof = msh->getNumberOfNodes();
        MathLib::DenseLinearEquations eqs;
        eqs.create(n_dof);

        MathLib::DenseLinearEquations::MatrixType* globalA = eqs.getA();
        double* globalRHS = eqs.getRHS();

        //assembly
        FemLib::LagrangianFeObjectContainer* feObjects = gw.head->getFeObjectContainer();
        MathLib::Matrix<double> localK;
        std::vector<size_t> e_node_id_list;
        double gp_x[3], real_x[3];
        for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
            MeshLib::IElement *e = msh->getElemenet(i_e);
            FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
            const size_t &n_dof = fe->getNumberOfVariables();
            localK.resize(n_dof, n_dof);
            localK = .0;
            FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
            for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            	q->getSamplingPoint(j, gp_x);
            	fe->computeBasisFunctions(gp_x);
                fe->getRealCoordinates(real_x);
            	double k;
            	gw._K->eval(real_x, k);
            	fe->integrateDWxDN(j, k, localK);
            }
            e->getNodeIDList(e_node_id_list);
            globalA->add(e_node_id_list, localK);
        }

        //outputLinearEQS(*globalA, globalRHS);

        //apply BC
        for (size_t i=0; i<gw.vec_bc2.size(); i++) gw.vec_bc2[i]->apply(globalRHS);
        //outputLinearEQS(globalA, globalRHS);
        for (size_t i=0; i<gw.vec_bc1.size(); i++) gw.vec_bc1[i]->apply(eqs);
        //outputLinearEQS(*globalA, globalRHS);

        //solve
        eqs.solve();

        //update head
        gw.head->setNodalValues(eqs.getX(), 0, eqs.getDimension());
    }

    static void calculateVelocity(GWFemTest &gw)
    {
        const MeshLib::IMesh *msh = gw.msh;
        FemLib::LagrangianFeObjectContainer* feObjects = gw.head->getFeObjectContainer();
        //calculate vel (vel=f(h))
        for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
            MeshLib::IElement* e = msh->getElemenet(i_e);
            FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
            std::vector<double> local_h(e->getNumberOfNodes());
            for (size_t j=0; j<e->getNumberOfNodes(); j++)
                local_h[j] = gw.head->getValue(e->getNodeID(j));
            // for each integration points
            FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
            double r[2] = {};
            const size_t n_gp = integral->getNumberOfSamplingPoints();
            gw.vel->setNumberOfIntegationPoints(i_e, n_gp);
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
                fe->getRealCoordinates(&xx[0]);

                double k;
                gw._K->eval(&xx[0], k);
                dN->axpy(-k, &local_h[0], .0, &q[0]); //TODO  q = - K * dN * local_h;
                gw.vel->setIntegrationPointValue(i_e, ip, q);
            }
        }
    }
};



struct DiscreteExample1
{
    std::vector<size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const size_t dim_eqs = 9;
    std::vector<double> exH;

    DiscreteExample1()
    {
        size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        list_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        list_dirichlet_bc_value.resize(6);
        fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
        fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);
        exH.resize(9);
        for (size_t i=0; i<9; i++) {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
    }

    void setLocalDirichletBC(const Base::BidirectionalMap<size_t, size_t> &map_global2localNodeId, std::vector<size_t> &local_dirichlet_bc_id, std::vector<double> &local_dirichlet_bc_value)
    {
        for (size_t i=0; i<list_dirichlet_bc_id.size(); i++) {
            if (map_global2localNodeId.countInA(list_dirichlet_bc_id[i])>0) {
                size_t local_id = map_global2localNodeId.mapAtoB(list_dirichlet_bc_id[i]);
                local_dirichlet_bc_id.push_back(local_id);
                local_dirichlet_bc_value.push_back(list_dirichlet_bc_value[i]);
            }
        }
    }


    class TestElementAssembler : public DiscreteLib::IElemenetWiseLinearEquationLocalAssembler
    {
        MathLib::Matrix<double> _m;
    public:
        TestElementAssembler()
        {
            _m.resize(4,4);
            _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0; 
            _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
            _m(2,2) = 4.0; _m(2,3) = -1.0;
            _m(3,3) = 4.0;
            for (size_t i=0; i<4; i++)
                for (size_t j=0; j<i; j++) _m(i,j) = _m(j,i);
            _m *= 1.e-11/6.0;
        }
        void assembly(MeshLib::IElement &e, MathLib::DenseLinearEquations &eqs)
        {
            (*eqs.getA()) = _m;
        }
    };
};
