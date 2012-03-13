
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
#include "DiscreteLib/OMPDiscreteSystem.h"
#include "DiscreteLib/DiscreteVector.h"
#include "DiscreteLib/DiscreteLinearEquation.h"
#include "DiscreteLib/ElementLocalAssembler.h"
#include "DiscreteLib/DoF.h"
#include "DiscreteLib/SparsityBuilder.h"
#include "DiscreteLib/ogs5/par_ddc_group.h"
#include "DiscreteLib/SerialNodeDecompositionDiscreteSystem.h"

#include "TestUtil.h"
#include "TestExamples.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace DiscreteLib;

TEST(Discrete, DDC1)
{
    // subdomain1
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DDCGlobaLocalMappingOffset mapping(0, msh->getNumberOfNodes(), 0);
    DDCSubDomain* dom1 = new DDCSubDomain(*msh, mapping);
    // global
    DDCGlobal ddc(DecompositionType::Node);
    ddc.addSubDomain(dom1);

    // discrete system
    SerialNodeDecompositionDiscreteSystem dis(ddc);
    
    // vector
    IDiscreteVector<double>* v = dis.createVector<double>(msh->getNumberOfNodes());
    for (size_t i=0; i<v->size(); ++i)
        (*v)[i] = i;
    ASSERT_EQ(9, v->size());
    double expected[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_DOUBLE_ARRAY_EQ(expected, *v, 9);

    // linear equation
    DiscreteExample1 ex1;
    DiscreteExample1::TestElementAssembler ele_assembler;
    CRSLisSolver lis;
    lis.getOption().ls_method = LIS_option::CG;
    lis.getOption().ls_precond = LIS_option::NONE;
    DofMapManager dofManager;
    dofManager.addDoF(msh->getNumberOfNodes());
    dofManager.construct(DofMapManager::BY_DOF);
    IDiscreteLinearEquation* linear_eq = dis.createLinearEquation<CRSLisSolver, SparsityBuilderFromNodeConnectivity>(lis, dofManager);
    linear_eq->setPrescribedDoF(0, ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
    linear_eq->construct(ElementBasedAssembler(ele_assembler));
    //linear_eq->getLinearEquation()->printout();
    linear_eq->solve();

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getLocalX(), 9, 1.e-5);
}

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

TEST(Discrete, OGS51)
{
    MeshLib::IMixedOrderMesh* msh;
    std::set<std::pair<bool, size_t>> set_property;
    OGS5::CPARDomainGroup dg(*msh, set_property);
}

//# OpenMP ###################################################################################################
TEST(Discrete, OMP_vec1)
{
    MeshLib::IMesh* org_msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    const size_t n_dom = 2;
    OMPGlobalDiscreteVector<double> global_v(10, n_dom);

    omp_set_num_threads(n_dom);

    #pragma omp parallel shared(global_v, org_msh, std::cout) default(none)
    {
        int iam = omp_get_thread_num();
        //std::cout << "threads " << iam << std::endl;

        OMPLocalDiscreteVector<double>* local_v = global_v.createLocal(iam, iam*5, (iam+1)*5);
        for (size_t i=local_v->getRangeBegin(); i<local_v->getRangeEnd(); i++)
            local_v->global(i) = iam*10 + i;
    }

    double expected[] = {0, 1, 2, 3, 4, 15, 16, 17, 18, 19};
    ASSERT_DOUBLE_ARRAY_EQ(expected, global_v, 10);
}


TEST(Discrete, OMP_eqs1)
{
    DiscreteExample1 ex1;
    struct NodeDDC
    {
        Base::BidirectionalMap<size_t, size_t> map_global2localNodeId;
        std::set<size_t> ghost_nodes;
    };
    MeshLib::IMesh *org_msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);
    const size_t n_dom = 2;
    OMPMasterNodeDecomposedDiscreteSystem global_dis(*org_msh, n_dom);

    omp_set_num_threads(n_dom);

    // define decomposed discrete system
    #pragma omp parallel shared(ex1, global_dis, org_msh, std::cout) default(none)
    {
        int iam = omp_get_thread_num();
        //std::cout << "threads " << iam << std::endl;

        MeshLib::IMesh* local_msh = 0;
        NodeDDC dom;
        if (iam==0) {
            // define a local problem with ddc
            int dom1_eles[] = {0, 1, 2, 3};
            int dom1_ghost_nodes[] = {5, 6, 7, 8};
            std::vector<size_t> dom1_e(dom1_eles, dom1_eles+4);

            MeshGenerator::generateSubMesh(*org_msh, dom1_e, local_msh, dom.map_global2localNodeId);
            for (size_t i=0; i<4; i++) {
                dom.ghost_nodes.insert(dom.map_global2localNodeId.mapAtoB(dom1_ghost_nodes[i]));
            }
        } else {
            int dom2_eles[] = {1, 2, 3};
            int dom2_ghost_nodes[] = {0, 1, 2, 3, 4};
            std::vector<size_t> dom2_e(dom2_eles, dom2_eles+3);
            MeshGenerator::generateSubMesh(*org_msh, dom2_e, local_msh, dom.map_global2localNodeId);
            for (size_t i=0; i<5; i++) {
                dom.ghost_nodes.insert(dom.map_global2localNodeId.mapAtoB(dom2_ghost_nodes[i]));
            }
        }
        std::vector<size_t> list_dirichlet_bc_id;
        std::vector<double> list_dirichlet_bc_value;
        ex1.setLocalDirichletBC(dom.map_global2localNodeId, list_dirichlet_bc_id, list_dirichlet_bc_value);

        // discrete system
        OMPLocalNodeDecomposedDiscreteSystem *local_dis = global_dis.createLocal(*local_msh, dom.map_global2localNodeId, dom.ghost_nodes);

        //// create dof map
        //DofMapManager dofManager;
        //size_t dofId = dofManager.addDoF(local_msh->getNumberOfNodes(), dom.ghost_nodes);
        //dofManager.construct(DofMapManager::BY_POINT);

        // create a linear problem
        CRSLisSolver lis;
        lis.getOption().ls_method = LIS_option::CG;
        lis.getOption().ls_precond = LIS_option::NONE;
        DofMapManager dofManager;
        size_t dofId = dofManager.addDoF(org_msh->getNumberOfNodes());
        dofManager.construct(DofMapManager::BY_POINT);
        IDiscreteLinearEquation *linear_eq = local_dis->createLinearEquation<CRSLisSolver, SparsityBuilderFromNodeConnectivityWithInactiveDoFs>(lis, dofManager);
        linear_eq->setPrescribedDoF(dofId, list_dirichlet_bc_id, list_dirichlet_bc_value);
        // construct and solve
        DiscreteExample1::TestElementAssembler ele_assembler;
        linear_eq->construct(ElementBasedAssembler(ele_assembler));    
    }


     IDiscreteLinearEquation *global_eq;
     global_eq->solve();
}


