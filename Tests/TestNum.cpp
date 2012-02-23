
#include <gtest/gtest.h>

#include <vector>

#include "Base/CodingTools.h"

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/SparseLinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"

#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"


#include "NumLib/Discrete/DiscreteSystem.h"
#include "NumLib/Discrete/DiscreteLinearEquation.h"
#include "NumLib/Discrete/DiscreteLinearEquationAssembler.h"
#include "NumLib/Discrete/ElementLocalAssembler.h"
#include "NumLib/Discrete/DoF.h"
#include "NumLib/Discrete/SparsityBuilder.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"
#include "NumLib/DomainDecomposition/DomainDecomposition.h"
#include "NumLib/DomainDecomposition/ogs5/par_ddc_group.h"

#include "NumLib/Output/Output.h"

#include "TestUtil.h"

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace NumLib;

TEST(Num, SingleDOF)
{
    DofMapManager dofManagerA;
    size_t dofId1 = dofManagerA.addDoF(10);
    dofManagerA.construct();
    const DofMap *dofMap1 = dofManagerA.getDofMap(dofId1); 

    ASSERT_EQ(dofManagerA.getNumberOfDof(), 1);
    ASSERT_EQ(dofManagerA.getTotalNumberOfDiscretePoints(), 10);
    ASSERT_TRUE(dofMap1!=0);
    ASSERT_EQ(dofMap1->getNumberOfDiscretePoints(), 10);
    ASSERT_EQ(dofMap1->getEqsID(0), 0);
    ASSERT_EQ(dofMap1->getEqsID(9), 9);
};

TEST(Num, NumberingDofByDof)
{
    DofMapManager dofManagerB;
    size_t dofIdB1 = dofManagerB.addDoF(10);
    size_t dofIdB2 = dofManagerB.addDoF(10);
    dofManagerB.construct();
    const DofMap *dofMapB1 = dofManagerB.getDofMap(dofIdB1); 
    const DofMap *dofMapB2 = dofManagerB.getDofMap(dofIdB2); 
    ASSERT_EQ(dofManagerB.getNumberOfDof(), 2);
    ASSERT_EQ(dofManagerB.getTotalNumberOfDiscretePoints(), 20);
    ASSERT_EQ(dofMapB1->getEqsID(0), 0);
    ASSERT_EQ(dofMapB1->getEqsID(9), 9);
    ASSERT_EQ(dofMapB2->getEqsID(0), 10);
    ASSERT_EQ(dofMapB2->getEqsID(9), 19);
};

TEST(Num, NumberingDofByPoint)
{
    DofMapManager dofManagerB;
    size_t dofIdB1 = dofManagerB.addDoF(10);
    size_t dofIdB2 = dofManagerB.addDoF(10);
    dofManagerB.construct(DofMapManager::BY_POINT);
    const DofMap *dofMapB1 = dofManagerB.getDofMap(dofIdB1); 
    const DofMap *dofMapB2 = dofManagerB.getDofMap(dofIdB2); 
    ASSERT_EQ(dofManagerB.getNumberOfDof(), 2);
    ASSERT_EQ(dofManagerB.getTotalNumberOfDiscretePoints(), 20);
    ASSERT_EQ(dofMapB1->getEqsID(0), 0);
    ASSERT_EQ(dofMapB1->getEqsID(9), 18);
    ASSERT_EQ(dofMapB2->getEqsID(0), 1);
    ASSERT_EQ(dofMapB2->getEqsID(9), 19);
}

TEST(Num, TimeStepping1)
{
    // define problems and solution strategy
    class TestProblem : public ITransientSystem
    {
    private:
        double _t0, _tn, _dt, _v;
    public:
        TestProblem(double t0, double tn, double dt) : _t0(t0), _tn(tn), _dt(dt), _v(.0) {};
        int solveTimeStep(const TimeStep &time) { _v+=1.0; return 0; }
        double suggestNext(const TimeStep &time_current)
        {
            if (time_current.getTime() < _t0) return _t0;
            else return (time_current.getTime()+_dt <= _tn) ? time_current.getTime()+_dt : time_current.getTime();
        }
        bool isAwake(const TimeStep &time)  { return (time.getTime()>=_t0 && time.getTime()<=_tn); }
        void accept(const TimeStep &time) {};
        double getValue() {return _v;};
    };
    TestProblem problem(.0, 100., 10.);

    // pass it to discretization systems
    TimeSteppingController timeStepping;
    timeStepping.addTransientSystem(problem);
    // start time stepping
    timeStepping.setBeginning(.0);
    size_t n_steps = timeStepping.solve(100.);
    //check
    ASSERT_EQ(10, n_steps);
    ASSERT_EQ(10., problem.getValue());
}

TEST(Num, OGS5DDC)
{
#if 0
    MeshLib::IMixedOrderMesh* msh;
    bool msh_order = false;

    std::vector<ITransientSystem*> problems;
    std::set<std::pair<bool,size_t>> eqs_properties;

    OGS5::CPARDomain *par;
    OGS5::CPARDomainGroup doms(*msh, eqs_properties);
    doms.addDomain(par);
    doms.setup();
    doms.solveTimeStep(100.);
#endif
}

struct NumExample1
{
    std::vector<size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const size_t dim_eqs = 9;
    std::vector<double> exH;

    NumExample1()
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

    class TestElementAssembler : public IElemenetLocalAssembler
    {
        Matrix<double> _m;
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

TEST(Discrete, Lis1)
{
    NumExample1 ex1;
    NumExample1::TestElementAssembler ele_assembler;
    CRSLisSolver lis;
    lis.getOption().ls_method = LIS_option::CG;
    lis.getOption().ls_precond = LIS_option::NONE;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);

    // define discrete system
    DiscreteSystem dis(*msh);
    {
        // create a linear problem
        IDiscreteLinearEquation *linear_eq = dis.createLinearEquation<CRSLisSolver, NumLib::SparsityBuilderFromNodeConnectivity>(lis);
        // DoF?
        DofMapManager *dofManager = linear_eq->getDofMapManger();
        dofManager->addDoF(msh->getNumberOfNodes());
        dofManager->construct(DofMapManager::BY_DOF);
        // solve the equation
        linear_eq->construct(NumLib::ElementBasedAssembler(ele_assembler));
        linear_eq->getLinearEquation()->setKnownX(ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
        linear_eq->solve();

        ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getX(), 9, 1.e-5);
    }

}

class SubDomain
{
public:
    std::vector<size_t> list_e;
    std::vector<size_t> list_n;
    MeshLib::IMesh *msh;
    size_t getDimension() {return list_n.size();};
};

class INodalDecomposedLinearEquation : public IDiscreteLinearEquation
{

};


template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
class TemplateNodalDecomposedLinearEquation : public INodalDecomposedLinearEquation
{
public:
    /// construct 
    void construct(IDiscreteLinearEquationAssembler& assemler)
    {
        DofMapManager* dofManager = getDofMapManger();
        assert(dofManager->getNumberOfDof()>0);

        if (_do_create_eqs) {
            _do_create_eqs = false;
            MathLib::RowMajorSparsity sparse;
            T_SPARSITY_BUILDER sp_builder(*getMesh(), *dofManager, sparse);
            getLinearEquation()->create(dofManager->getTotalNumberOfDiscretePoints(), &sparse);
        } else {
            getLinearEquation()->reset();
        }
        assemler.assembly(*getMesh(), *dofManager, *getLinearEquation());
    }

    /// solve
    void solve()
    {
        _eqs->solve();
    }


    /// get the solution vector
    double* getX()
    {
        return _eqs->getX();
    }

    /// get the RHS vector
    double* getRHS()
    {
        return _eqs->getRHS();
    }
    /// get a linear equation object
    MathLib::ILinearEquations* getLinearEquation() const
    {
        return _eqs;
    }
    /// get a Dof map manager
    DofMapManager* getDofMapManger() const
    {
        return _dofManager;
    }
private:
    T_LINEAR_SOLVER* _eqs;
    DofMapManager* _dofManager;

};

class DecomposedDiscreteSystem
{
public:
    DecomposedDiscreteSystem(SubDomain &subdomain) {};
    /// create a new linear equation
    template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    IDiscreteLinearEquation* createLinearEquation(T_LINEAR_SOLVER &linear_solver)
    {
        _vec_linear_sys.push_back(new TemplateNodalDecomposedLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_subdomain, linear_solver));
        return _vec_linear_sys.back();
    }

    IDiscreteLinearEquation* getLinearEquation(size_t i)
    {
        return _vec_linear_sys[i];
    }
private:
    DISALLOW_COPY_AND_ASSIGN(DecomposedDiscreteSystem);

    SubDomain* _subdomain;
    std::vector<INodalDecomposedLinearEquation*> _vec_linear_sys;
};

TEST(Discrete, Lis2)
{
    NumExample1 ex1;
    NumExample1::TestElementAssembler ele_assembler;
    CRSLisSolver lis;
    lis.getOption().ls_method = LIS_option::CG;
    lis.getOption().ls_precond = LIS_option::NONE;
    MeshLib::IMesh *msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);

    SubDomain dom1;
    int dom1_nodes[] = {0, 1, 2, 3, 4};
    int dom1_eles[] = {0, 1, 2, 3};
    dom1.list_n.assign(dom1_nodes, dom1_nodes+5);
    dom1.list_e.assign(dom1_eles, dom1_eles+4);
    //dom1.msh = msh.getSubMesh(0, 1);

    SubDomain dom2;
    int dom2_nodes[] = {5, 6, 7, 8};
    int dom2_eles[] = {0, 1, 2, 3};
    dom2.list_n.assign(dom2_nodes, dom2_nodes+4);
    dom2.list_e.assign(dom2_eles, dom2_eles+4);
    //dom2 = msh.getSubMesh(2, 3);


    {
        // define discrete system
        DecomposedDiscreteSystem dis(dom1);
        // create a linear problem
        IDiscreteLinearEquation *linear_eq = dis.createLinearEquation<CRSLisSolver, NumLib::SparsityBuilderFromNodeConnectivity>(lis);
        // DoF?
        DofMapManager *dofManager = linear_eq->getDofMapManger();
        dofManager->addDoF(msh->getNumberOfNodes());
        dofManager->construct(DofMapManager::BY_DOF);
        // solve the equation
        linear_eq->construct(NumLib::ElementBasedAssembler(ele_assembler));
        linear_eq->getLinearEquation()->setKnownX(ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
        linear_eq->solve();

        ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getX(), 9, 1.e-5);
    }

    {
        // define discrete system
        DiscreteSystem dis(*dom2.msh);
        // create a linear problem
        IDiscreteLinearEquation *linear_eq = dis.createLinearEquation<CRSLisSolver, NumLib::SparsityBuilderFromNodeConnectivity>(lis);
        // DoF?
        DofMapManager *dofManager = linear_eq->getDofMapManger();
        dofManager->addDoF(msh->getNumberOfNodes());
        dofManager->construct(DofMapManager::BY_DOF);
        // solve the equation
        linear_eq->construct(NumLib::ElementBasedAssembler(ele_assembler));
        linear_eq->getLinearEquation()->setKnownX(ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
        linear_eq->solve();

        ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getX(), 9, 1.e-5);
    }

}

