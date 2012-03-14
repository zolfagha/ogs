
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

#include "TestUtil.h"
#include "TestExamples.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace DiscreteLib;




//# DoF ###################################################################################################
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
        dofManager->addDoF(9, &vec_ghost_nodes);
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
        dofManager->addDoF(8, &vec_ghost_nodes);
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
