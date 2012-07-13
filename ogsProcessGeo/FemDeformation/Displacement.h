
#pragma once

#include "MathLib/LinAlg/VectorNorms.h"
#include "GeoLib/Shape/Rectangle.h"
#include "NumLib/Function/Function.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "FemLib/Function/FemFunction.h"
#include "SolutionLib/FemProblem/FemEquation.h"
#include "SolutionLib/FemProblem/FemIVBVProblem.h"
#include "SolutionLib/Solution/SingleStepFEM.h"

#include "Tests/Geo/Equation/FemLinearElasticity.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Material/Solid.h"

using namespace GeoLib;
using namespace MathLib;
using namespace FemLib;
using namespace NumLib;
using namespace MeshLib;
using namespace SolutionLib;
using namespace DiscreteLib;

namespace Geo
{

typedef TemplateFemEquation<
		Geo::FemLinearElasticLinearLocalAssembler,
		Geo::FemLinearElasticResidualLocalAssembler,
		Geo::FemLinearElasticJacobianLocalAssembler
		>
		FemLinearElasticEquation;

typedef FemIVBVProblem< FemLinearElasticEquation > FemLinearElasticProblem;



template <
	class T_LINEAR_SOLVER
	>
class FunctionDisplacement : public TemplateTransientMonolithicSystem
{
    enum Out { u_x=0, u_y=1, Strain=2, Stress=3 };
public:
    typedef SingleStepFEM
    		<
    			FemLinearElasticProblem,
    			T_LINEAR_SOLVER
    		> MySolution;

    FunctionDisplacement()
    {
        TemplateTransientMonolithicSystem::resizeOutputParameter(4);
    };

    void define(DiscreteSystem* dis, FemLinearElasticProblem* problem_u, Geo::PorousMedia* pm, BaseLib::Options &option)
    {
    	_dis = dis;
    	_pm = pm;
        //solution algorithm
        _sol_u = new MySolution(dis, problem_u);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        typename MySolution::LinearSolverType* linear_solver = _sol_u->getLinearEquationSolver();
        linear_solver->setOption(option);
        _sol_u->getNonlinearSolver()->setOption(option);
        this->setOutput(u_x, problem_u->getVariable(0)->getIC());
        this->setOutput(u_y, problem_u->getVariable(1)->getIC());

        _strain = new FemLib::FEMIntegrationPointFunctionVector(*dis);
        _stress = new FemLib::FEMIntegrationPointFunctionVector(*dis);
    }

    void calculateStressStrain()
    {
        const MeshLib::IMesh *msh = _dis->getMesh();
        FemLib::FemNodalFunctionScalar *ux = _sol_u->getCurrentSolution(0);
        FemLib::FemNodalFunctionScalar *uy = _sol_u->getCurrentSolution(1);
        FemLib::FEMIntegrationPointFunctionVector* strain = _strain;
        FemLib::FEMIntegrationPointFunctionVector* stress = _stress;
        FemLib::LagrangianFeObjectContainer* feObjects = ux->getFeObjectContainer();

		const size_t dim = msh->getDimension();
		const size_t n_strain_components = (dim==2 ? 4 : 6);
		Solid *solidphase = _pm->solidphase;

        // set D
		NumLib::LocalMatrix matD = NumLib::LocalMatrix::Zero(n_strain_components, n_strain_components);
        //matD *= .0;
        double nv = solidphase->poisson_ratio;
        double E = solidphase->Youngs_modulus;
        double Lambda, G, K;
        calculateLameConstant(nv, E, Lambda, G, K);
        setElasticConsitutiveTensor(dim, Lambda, G, matD);

        //calculate strain, stress
        for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++)
        {
        	// element setup
            MeshLib::IElement* e = msh->getElemenet(i_e);
    		const size_t nnodes = e->getNumberOfNodes();
            FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
            FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
            const size_t n_gp = integral->getNumberOfSamplingPoints();
            strain->setNumberOfIntegationPoints(i_e, n_gp);
            stress->setNumberOfIntegationPoints(i_e, n_gp);

            // local u
            NumLib::LocalVector local_u(dim*nnodes);
            for (size_t j=0; j<nnodes; j++) {
            	const size_t node_id = e->getNodeID(j);
				local_u[j*dim] = ux->getValue(node_id);
				local_u[j*dim+1] = uy->getValue(node_id);
            }

            // for each integration points
            NumLib::LocalMatrix matB(n_strain_components, nnodes*dim);
    		NumLib::LocalMatrix matN(dim, nnodes*dim);
            double r[3] = {};
            double x[3] = {};
            for (size_t ip=0; ip<n_gp; ip++) {
                integral->getSamplingPoint(ip, r);
                fe->computeBasisFunctions(r);
                NumLib::LocalMatrix &N = *fe->getBasisFunction();
                const NumLib::LocalMatrix &dN = *fe->getGradBasisFunction();
                fe->getRealCoordinates(x);

                // set N,B
                setNu_Matrix(dim, nnodes, N, matN);
                setB_Matrix(dim, nnodes, dN, matB);

                // strain
                NumLib::LocalVector gp_strain(n_strain_components);
                gp_strain.noalias() = matB * local_u;
                strain->setIntegrationPointValue(i_e, ip, gp_strain);

                // stress
                NumLib::LocalVector gp_stress(n_strain_components);
                gp_stress = matD * gp_strain;
                stress->setIntegrationPointValue(i_e, ip, gp_stress);

//                std::cout << "strain=\n" << gp_strain << std::endl;
//                std::cout << "D=\n" << matD << std::endl;
//                std::cout << "stress=\n" << gp_stress << std::endl;
            }
        }

        setOutput(Strain, strain);
        setOutput(Stress, stress);
    }

    int solveTimeStep(const TimeStep &time)
    {
        _sol_u->solveTimeStep(time);
        setOutput(u_x, _sol_u->getCurrentSolution(0));
        setOutput(u_y, _sol_u->getCurrentSolution(1));

        calculateStressStrain();

        return 0;
    }

    double suggestNext(const TimeStep &time_current) { return _sol_u->suggestNext(time_current); }

    bool isAwake(const TimeStep &time) { return _sol_u->isAwake(time);  }

    void accept(const TimeStep &time)
    {
        _sol_u->accept(time);

        //std::cout << "Head=" << std::endl;
        //_solHead->getCurrentSolution(0)->printout();
    };

private:
    //FemLinearElasticProblem* _problem;
    MySolution* _sol_u;
    DiscreteSystem* _dis;
    Geo::PorousMedia* _pm;
    LagrangianFeObjectContainer* _feObjects;
    FemLib::FEMIntegrationPointFunctionVector* _strain;
    FemLib::FEMIntegrationPointFunctionVector* _stress;

    DISALLOW_COPY_AND_ASSIGN(FunctionDisplacement);
};


} //end


