
#include <gtest/gtest.h>

#include <vector>

#include "Base/CodingTools.h"
#include "Base/BidirectionalMap.h"

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/SparseLinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"

#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"


#include "DiscreteLib/DiscreteSystem.h"
#include "DiscreteLib/DiscreteVector.h"
#include "DiscreteLib/DiscreteLinearEquation.h"
#include "DiscreteLib/ElementLocalAssembler.h"
#include "DiscreteLib/DoF.h"
#include "DiscreteLib/SparsityBuilder.h"

#include "TestUtil.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace DiscreteLib;

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

TEST(Discrete, VecSingle1)
{
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    DiscreteVector<double> *v = dis.createVector<double>(msh->getNumberOfNodes());

    double expected[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    std::copy(expected, expected+9, v->begin());
    
    ASSERT_EQ(9, v->size());
    ASSERT_DOUBLE_ARRAY_EQ(expected, *v, 9);
}

TEST(Discrete, VecDecomposed1)
{
    MeshLib::IMesh* org_msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    const size_t n_dom = 2;
    std::vector<NodeDecomposedDiscreteSystem*> vec_dis(n_dom);
    DecomposedMasterVector<double> global_v(10);
    global_v.decompose(n_dom);

    omp_set_num_threads(n_dom);

    #pragma omp parallel shared(global_v, org_msh, vec_dis, std::cout) default(none)
    {
        int iam = omp_get_thread_num();
        //std::cout << "threads " << iam << std::endl;

        DecomposedLocalVector<double>* local_v = global_v.createLocal(iam, iam*5, (iam+1)*5);
        for (size_t i=local_v->getRangeBegin(); i<local_v->getRangeEnd(); i++)
            local_v->global(i) = iam*10 + i;
    }

    double expected[] = {0, 1, 2, 3, 4, 15, 16, 17, 18, 19};
    ASSERT_DOUBLE_ARRAY_EQ(expected, global_v, 10);
}

TEST(Discrete, VecDecomposed2)
{
#if 0
    MeshLib::IMesh* org_msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    const size_t n_dom = 2;
    std::vector<NodeDecomposedDiscreteSystem*> vec_dis(n_dom);
    DecomposedDiscreteVector<double> global_v(10);
    global_v.decompose(n_dom);

    omp_set_num_threads(n_dom);

#pragma omp parallel shared(global_v, org_msh, vec_dis, std::cout) default(none)
    {
        int iam = omp_get_thread_num();
        std::cout << "threads " << iam << std::endl;

        DecomposedDiscreteVectorLeaf<double>* local_v = global_v.createLocal(iam, iam*5, (iam+1)*5);
        for (size_t i=local_v->getRangeBegin(); i<local_v->getRangeEnd(); i++)
            (*local_v)[i] = iam*10 + i;
        Base::BidirectionalMap<size_t, size_t> msh_node_id_mapping;
        std::set<size_t> ghost_nodes;
        MeshLib::IMesh* local_msh;
        if (iam==0) {
            int dom1_eles[] = {0, 1, 2, 3};
            int dom1_ghost_nodes[] = {5, 6, 7, 8};
            std::vector<size_t> dom1_e(dom1_eles, dom1_eles+4);
            MeshGenerator::generateSubMesh(*org_msh, dom1_e, local_msh, msh_node_id_mapping);
            for (size_t i=0; i<4; i++) {
                ghost_nodes.insert(msh_node_id_mapping.mapAtoB(dom1_ghost_nodes[i]));
            }
        } else {
            int dom2_eles[] = {1, 2, 3};
            int dom2_ghost_nodes[] = {0, 1, 2, 3, 4};
            std::vector<size_t> dom2_e(dom2_eles, dom2_eles+3);
            MeshGenerator::generateSubMesh(*org_msh, dom2_e, local_msh, msh_node_id_mapping);
            for (size_t i=0; i<5; i++) {
                ghost_nodes.insert(msh_node_id_mapping.mapAtoB(dom2_ghost_nodes[i]));
            }
        }

        NodeDecomposedDiscreteSystem* dis = new NodeDecomposedDiscreteSystem(*local_msh, msh_node_id_mapping, ghost_nodes);
        vec_dis[iam] = dis;
        DecomposedDiscreteVector<double> *v = dis->createVector<double>(local_msh->getNumberOfNodes());
        for (size_t i=v->getRangeBegin(); i<v->getRangeEnd(); i++)
            (*v)[i] = iam*10 + i;
    }

    for (size_t i=0; i<global_v.size(); i++) {
        std::cout << global_v[i] << " ";
    }
    std::cout << std::endl;
#endif
}

TEST(Discrete, Lis1)
{
    DiscreteExample1 ex1;
    DiscreteExample1::TestElementAssembler ele_assembler;
    CRSLisSolver lis;
    lis.getOption().ls_method = LIS_option::CG;
    lis.getOption().ls_precond = LIS_option::NONE;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);

    // define discrete system
    DiscreteSystem dis(*msh);
    {
        // DoF?
        DofMapManager dofManager;
        dofManager.addDoF(msh->getNumberOfNodes());
        dofManager.construct(DofMapManager::BY_DOF);
        // create a linear problem
        IDiscreteLinearEquation *linear_eq = dis.createLinearEquation<CRSLisSolver, SparsityBuilderFromNodeConnectivity>(lis, dofManager);
        // solve the equation
        linear_eq->setPrescribedDoF(0, ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
        linear_eq->construct(ElementBasedAssembler(ele_assembler));
        //linear_eq->getLinearEquation()->printout();
        linear_eq->solve();

        ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getLocalX(), 9, 1.e-5);
    }

}


TEST(Discrete, Lis2)
{
    struct NodeDDC
    {
        Base::BidirectionalMap<size_t, size_t> map_global2localNodeId;
        std::set<size_t> ghost_nodes;
    };
    DiscreteExample1 ex1;
    DiscreteExample1::TestElementAssembler ele_assembler;
    MeshLib::IMesh *org_msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);

    {
        // define a local problem with ddc
        MeshLib::IMesh* local_msh = 0;
        NodeDDC dom1;
        {
            int dom1_eles[] = {0, 1, 2, 3};
            int dom1_ghost_nodes[] = {5, 6, 7, 8};
            std::vector<size_t> dom1_e(dom1_eles, dom1_eles+4);

            MeshGenerator::generateSubMesh(*org_msh, dom1_e, local_msh, dom1.map_global2localNodeId);
            for (size_t i=0; i<4; i++) {
                dom1.ghost_nodes.insert(dom1.map_global2localNodeId.mapAtoB(dom1_ghost_nodes[i]));
            }
        }
        std::vector<size_t> list_dirichlet_bc_id;
        std::vector<double> list_dirichlet_bc_value;
        ex1.setLocalDirichletBC(dom1.map_global2localNodeId, list_dirichlet_bc_id, list_dirichlet_bc_value);
        CRSLisSolver lis;
        lis.getOption().ls_method = LIS_option::CG;
        lis.getOption().ls_precond = LIS_option::NONE;

        // define discrete system
        NodeDecomposedDiscreteSystem dis(*local_msh, dom1.map_global2localNodeId, dom1.ghost_nodes);
        // create dof map
        DofMapManager dofManager;
        size_t dofId = dofManager.addDoF(local_msh->getNumberOfNodes(), dom1.ghost_nodes);
        dofManager.construct(DofMapManager::BY_POINT);
        // create a linear problem
        IDiscreteLinearEquation *linear_eq = dis.createLinearEquation<CRSLisSolver, SparsityBuilderFromNodeConnectivityWithInactiveDoFs>(lis, dofManager);
        linear_eq->setPrescribedDoF(dofId, list_dirichlet_bc_id, list_dirichlet_bc_value);
        linear_eq->construct(ElementBasedAssembler(ele_assembler));
        //linear_eq->getLinearEquation()->printout();
        
        //solve
        //linear_eq->solve();

        //ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getX(), 9, 1.e-5);
    }

    {
        // define a local problem with ddc
        MeshLib::IMesh* local_msh = 0;
        NodeDDC dom2;
        {
            int dom2_eles[] = {1, 2, 3};
            int dom2_ghost_nodes[] = {0, 1, 2, 3, 4};
            std::vector<size_t> dom2_e(dom2_eles, dom2_eles+3);
            MeshGenerator::generateSubMesh(*org_msh, dom2_e, local_msh, dom2.map_global2localNodeId);
            for (size_t i=0; i<5; i++) {
                dom2.ghost_nodes.insert(dom2.map_global2localNodeId.mapAtoB(dom2_ghost_nodes[i]));
            }
        }
        std::vector<size_t> list_dirichlet_bc_id;
        std::vector<double> list_dirichlet_bc_value;
        ex1.setLocalDirichletBC(dom2.map_global2localNodeId, list_dirichlet_bc_id, list_dirichlet_bc_value);
        CRSLisSolver lis;
        lis.getOption().ls_method = LIS_option::CG;
        lis.getOption().ls_precond = LIS_option::NONE;

        // define discrete system
        NodeDecomposedDiscreteSystem dis(*local_msh, dom2.map_global2localNodeId, dom2.ghost_nodes);
        // create a dof map
        DofMapManager dofManager;
        size_t dofId = dofManager.addDoF(local_msh->getNumberOfNodes(), dom2.ghost_nodes);
        dofManager.construct(DofMapManager::BY_POINT);
        // create a linear problem
        IDiscreteLinearEquation *linear_eq = dis.createLinearEquation<CRSLisSolver, SparsityBuilderFromNodeConnectivityWithInactiveDoFs>(lis, dofManager);
        linear_eq->setPrescribedDoF(dofId, list_dirichlet_bc_id, list_dirichlet_bc_value);
        linear_eq->construct(ElementBasedAssembler(ele_assembler));
        //linear_eq->printout();
        // solve the equation
        //linear_eq->solve();

        //ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getX(), 9, 1.e-5);    
    }
}

TEST(Discrete, DoF_single)
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

TEST(Discrete, DoF_numberingByDof)
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

TEST(Discrete, DoF_numberingByPoint)
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

TEST(Discrete, DoF_ghost_nodes)
{
    {
        DofMapManager *dofManager = new DofMapManager();
        int ghost_nodes[] = {5, 6, 7, 8};
        std::set<size_t> vec_ghost_nodes(ghost_nodes, ghost_nodes+4);
        dofManager->addDoF(9, vec_ghost_nodes);
        dofManager->construct(DofMapManager::BY_POINT);
        const DofMap *dofMap = dofManager->getDofMap(0);

        ASSERT_EQ(1, dofManager->getNumberOfDof());
        ASSERT_EQ(9, dofManager->getTotalNumberOfDiscretePoints());
        ASSERT_EQ(5, dofManager->getTotalNumberOfActiveDoFs());
        ASSERT_TRUE(dofMap!=0);
        ASSERT_EQ(9, dofMap->getNumberOfDiscretePoints());
        ASSERT_EQ(5, dofMap->getNumberOfActiveDoFs());
        ASSERT_EQ(0, dofMap->getEqsID(0));
        ASSERT_EQ(-1, dofMap->getEqsID(8));

        delete dofManager;
    }
    {
        DofMapManager *dofManager = new DofMapManager();
        int ghost_nodes[] = {0, 1, 2, 3, 4};
        std::set<size_t> vec_ghost_nodes(ghost_nodes, ghost_nodes+4);
        dofManager->addDoF(8, vec_ghost_nodes);
        dofManager->construct(DofMapManager::BY_POINT);

        ASSERT_EQ(1, dofManager->getNumberOfDof());
        ASSERT_EQ(8, dofManager->getTotalNumberOfDiscretePoints());
        ASSERT_EQ(4, dofManager->getTotalNumberOfActiveDoFs());
        const DofMap *dofMap = dofManager->getDofMap(0);
        ASSERT_TRUE(dofMap!=0);
        ASSERT_EQ(8, dofMap->getNumberOfDiscretePoints());
        ASSERT_EQ(4, dofMap->getNumberOfActiveDoFs());
        ASSERT_EQ(-1, dofMap->getEqsID(0));
        ASSERT_EQ(3, dofMap->getEqsID(7));
        delete dofManager;
    }
}
