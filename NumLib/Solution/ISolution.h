//
//#pragma once
//
//#include "MathLib/LinAlg/Sparse/CRSMatrix.h"
//#include "MeshLib/Core/IMesh.h"
//
//#include "FemLib/Function/FemFunction.h"
//#include "FemLib/BC/FemDirichletBC.h"
//#include "FemLib/BC/FemNeumannBC.h"
//
//#include "NumLib/Core/TimeStep.h"
//#include "NumLib/Core/LinearEquation.h"
//#include "NumLib/Core/LinearSolver.h"
//#include "NumLib/Core/DoF.h"
//
//typedef MathLib::CRSMatrix<double,size_t> SparseMatrix;
//typedef std::vector<double> SparseVector;
//typedef MathLib::Matrix<double> DenseMatrix;
//typedef std::vector<double> DenseVector;
//
//namespace NumLib
//{
//class ISolution
//{};
//
//
//class TimeEulerSpFEMLinearSolution : public ISolution
//{
//	//Linear EQS
//	DiscretizedEQS<SparseMatrix, SparseVector> *_linearEQS;
//	LinearSolver<SparseMatrix, SparseVector> *_linearSolver;
//	// Mesh
//	DoFManager _dofManager;
//	//IC,BC,ST
//	FemLib::FemNodalFunctionScalar* _u0;
//    FemLib::FemDirichletBC<double> *_bc;
//    FemLib::FemNeumannBC<double> *_st;
//	//TIM
//	TimeStepping *_tim;
//
//	//
//	bool isAccepted;
//public:
//
//	void initialize(int *input, FemLib::FemNodalFunctionScalar* u0) 
//    {
//		//# Assume: A mesh does not change during simulation #
//		const MeshLib::IMesh *mesh = u0->getMesh();
//		//_dofManager = u0->getDOFManager();
//        this->_u0 = u0;
//        GeoLib::GeoObject *obj1;
//        MathLib::IFunction<double, GeoLib::Point> *f;
//        _bc = new FemLib::FemDirichletBC<double>(u0, obj1, f, 0);
//        _st = new FemLib::FemNeumannBC<double>(u0, obj1, f);
//		// TIM
//		_tim->initialize();
//		// EQS, solver
//		_linearEQS->initialize();
//		_linearSolver->initialize();
//	}
//	
//	TimeStep suggestNextTimeStep(TimeStep t_n) 
//    {
//		return _tim->next();
//	}
//	
//	bool isTimeStepAccepted() 
//    {
//		return isAccepted;
//	}
//	
//	FemLib::FemNodalFunctionScalar* solve(TimeStep t_n1, TimeStep t_n, FemLib::FemNodalFunctionScalar &u_n) 
//    {
//		FemLib::FemNodalFunctionScalar *u_n1 = new FemLib::FemNodalFunctionScalar(u_n); //copy mesh, etc..
//		
//		// collect data
//		TimeStep delta_t = t_n1 - t_n;
//		
//		// initialization
//		isAccepted = false;
//		_linearEQS->reset();
//        double theta = 1.0;
//		
//		// pre
//		doPreAssembly();
//		
//		// assembly
//        const MeshLib::IMesh *msh = u_n1->getMesh();
//        std::vector<double> local_u_n;
//        std::vector<size_t> dofmap;
//        for (size_t i=0; i<msh->getNumberOfElements(); i++) {
//            MeshLib::IElement *e = msh->getElemenet(i);
//            FemLib::IFiniteElement *fe = u_n1->getFiniteElement(e);
//			_dofManager.getMap(e, dofmap);
//			//local_u_n = u_n.getLocal(dofmap);
//			//
//			DiscretizedEQS<MathLib::Matrix<double>, std::vector<double>> localEQS;
//			MathLib::Matrix<double> M, K;
//			std::vector<double> F;
//			doLocalAssemblyEuler(fe, M, K, F);
//			//localEQS.getA() = 1.0/delta_t*M + theta*K;
//			//localEQS.getRHS() = (1.0/delta_t*M - (1.-theta)*K)*local_u_n + F;
//			
//			_linearEQS->add(localEQS, dofmap);
//		}
//		
//		// BC/ST
//		//if (!_st->isConstant())
//		//	_st.NodesValueList = doSetST();
//		//_linearEQS->addRHS(_st.NodesValueList, dofmap);
//		//if (!_bc->isConstant())
//		//	_bc.NodesValueList = doSetDirectBC();
//		//_linearEQS->setKnownX(_bc.NodesValueList, dofmap);
//		
//		// solve
//		_linearSolver->solve(_linearEQS);
//		u_n1->setNodalValues(_linearEQS->getX());
//
//		// 
//		doRightAfterSolvingPrimaryVariable();
//		
//		// check 
//		isAccepted = checkTimestep();
//		
//		if (isAccepted) {
//			// post
//			doPostTimeStep();
//		} else {
//			_tim->setNextTimeStep(t_n + delta_t/10.0);
//		}
//		
//		return u_n1;
//	}
//
//    void doLocalAssemblyEuler(FemLib::IFiniteElement *fe, MathLib::Matrix<double> &M, MathLib::Matrix<double> &K, std::vector<double> &F)
//    {
//
//    }
//
//    bool checkTimestep();
//
//    void doPreAssembly();
//
//    void doPostTimeStep();
//
//    void doRightAfterSolvingPrimaryVariable();
//
//};
//
//}
