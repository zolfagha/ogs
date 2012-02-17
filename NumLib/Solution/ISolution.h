
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/Sparse/CRSMatrix.h"
#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"

#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/Discrete/DoF.h"
#include "NumLib/Discrete/DiscreteSystem.h"

namespace NumLib
{

class ISolutionAlgorithm
{};


template<class T_USER_ASSEMBLY>
class TimeEulerSpFEMLinearSolution : public ISolutionAlgorithm
{
private:
    TimeODEElementAssembler<T_USER_ASSEMBLY> _element_ode_assembler;
    DiscreteSystem *_discrete_system;
	DofMapManager _dofManager;
	//IC,BC,ST
	FemLib::FemNodalFunctionScalar* _u0;

	//
	bool isAccepted;

public:
    TimeEulerSpFEMLinearSolution() {};

	void initialize(FemLib::FemNodalFunctionScalar* u0) 
    {
		//# Assume: A mesh does not change during simulation #
		const MeshLib::IMesh *mesh = u0->getMesh();
		//_dofManager = u0->getDOFManager();
        this->_u0 = u0;
        GeoLib::GeoObject *obj1;
        MathLib::IFunction<double, GeoLib::Point> *f;
	}
	
	double suggestNextTimeStep(TimeStep t_n) 
    {
		return .0;
	}
	
	bool isTimeStepAccepted() 
    {
		return isAccepted;
	}
	
	FemLib::FemNodalFunctionScalar* solve(const TimeStep &t_n1, FemLib::FemNodalFunctionScalar &u_n) 
    {
		FemLib::FemNodalFunctionScalar *u_n1 = new FemLib::FemNodalFunctionScalar(u_n); //copy mesh, etc..
		
		// collect data
		double delta_t = t_n1.getTimeStep();
        double theta = 1.0;
		
		// initialization
		isAccepted = false;
		
		// pre
		doPreAssembly();
		
		// assembly
        const MeshLib::IMesh *msh = u_n1->getMesh();
        IDiscreteLinearEquation* discreteEqs = _discrete_system->getLinearEquation();
        discreteEqs->construct(ElementBasedAssembler(TimeEulerElementAssembler<T_USER_ASSEMBLY>(t_n1, theta, _element_ode_assembler)));

		// BC/ST
		//if (!_st->isConstant())
		//	_st.NodesValueList = doSetST();
		//_linearEQS->addRHS(_st.NodesValueList, dofmap);
		//if (!_bc->isConstant())
		//	_bc.NodesValueList = doSetDirectBC();
		//_linearEQS->setKnownX(_bc.NodesValueList, dofmap);
		
		// solve
		discreteEqs->solve();

		//u_n1->setNodalValues(_linearEQS->getX());
        discreteEqs->getX();

		// 
		doRightAfterSolvingPrimaryVariable();
		
		// check 
		isAccepted = checkTimestep();
		
		if (isAccepted) {
			// post
			doPostTimeStep();
		} else {
			//_tim->setNextTimeStep(t_n + delta_t/10.0);
		}
		
		return u_n1;
	}

    bool checkTimestep();

    void doPreAssembly();

    void doPostTimeStep();

    void doRightAfterSolvingPrimaryVariable();

};

}
