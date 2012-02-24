
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
        // DoF?
        DofMapManager dofManager;
        dofManager.addDoF(msh->getNumberOfNodes());
        dofManager.construct(DofMapManager::BY_DOF);
        // create a linear problem
        IDiscreteLinearEquation *linear_eq = dis.createLinearEquation<CRSLisSolver, NumLib::SparsityBuilderFromNodeConnectivity>(lis, dofManager);
        // solve the equation
        linear_eq->setPrescribedDoF(0, ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
        linear_eq->construct(NumLib::ElementBasedAssembler(ele_assembler));
        //linear_eq->getLinearEquation()->printout();
        linear_eq->solve();

        ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getLocalX(), 9, 1.e-5);
    }

}

struct NodeDDC
{
    Base::BidirectionalMap<size_t, size_t> map_global2localNodeId;
    std::set<size_t> ghost_nodes;
};

TEST(Discrete, Lis2)
{
    NumExample1 ex1;
    NumExample1::TestElementAssembler ele_assembler;
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
        IDiscreteLinearEquation *linear_eq = dis.createLinearEquation<CRSLisSolver, NumLib::SparsityBuilderFromNodeConnectivityWithInactiveDoFs>(lis, dofManager);
        linear_eq->setPrescribedDoF(dofId, list_dirichlet_bc_id, list_dirichlet_bc_value);
        linear_eq->construct(NumLib::ElementBasedAssembler(ele_assembler));
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
            int dom2_ghost_nodes[] = {1, 2, 3, 4};
            std::vector<size_t> dom2_e(dom2_eles, dom2_eles+3);
            MeshGenerator::generateSubMesh(*org_msh, dom2_e, local_msh, dom2.map_global2localNodeId);
            for (size_t i=0; i<4; i++) {
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
        IDiscreteLinearEquation *linear_eq = dis.createLinearEquation<CRSLisSolver, NumLib::SparsityBuilderFromNodeConnectivityWithInactiveDoFs>(lis, dofManager);
        linear_eq->setPrescribedDoF(dofId, list_dirichlet_bc_id, list_dirichlet_bc_value);
        linear_eq->construct(NumLib::ElementBasedAssembler(ele_assembler));
        //linear_eq->printout();
        // solve the equation
        //linear_eq->solve();

        //ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getX(), 9, 1.e-5);    
    }
}

