/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestDoFs.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

#include <vector>

#include "BaseLib/CodingTools.h"
#include "BaseLib/BidirectionalMap.h"

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquation/SparseLinearEquation.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/LinearEquation/LisLinearEquation.h"
#endif

#include "GeoLib/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"


#include "DiscreteLib/Utils/DofEquationIdTable.h"

#include "TestUtil.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace DiscreteLib;




//# DoF ###################################################################################################
TEST(Discrete, DoF_var1)
{
    DofEquationIdTable dofManagerA;
    size_t varId = dofManagerA.addVariableDoFs(0, 0, 10);
    dofManagerA.construct();

    ASSERT_EQ(1, dofManagerA.getNumberOfVariables());
    ASSERT_EQ(1, dofManagerA.getNumberOfMeshes());
    //ASSERT_EQ(dofMap1->getNumberOfDiscretePoints(), 10);
    ASSERT_EQ(0, dofManagerA.mapEqsID(varId, 0, 0));
    ASSERT_EQ(9, dofManagerA.mapEqsID(varId, 0, 9));
};

TEST(Discrete, DoF_var2_byDof)
{
    DofEquationIdTable dofManagerB;
    size_t dofIdB1 = dofManagerB.addVariableDoFs(0, 0, 10);
    size_t dofIdB2 = dofManagerB.addVariableDoFs(0, 0, 10);
    dofManagerB.construct();
    ASSERT_EQ(2, dofManagerB.getNumberOfVariables());
    //ASSERT_EQ(dofManagerB.getTotalNumberOfDiscretePoints(), 20);
    ASSERT_EQ(0, dofManagerB.mapEqsID(dofIdB1, 0, 0));
    ASSERT_EQ(9, dofManagerB.mapEqsID(dofIdB1, 0, 9));
    ASSERT_EQ(10, dofManagerB.mapEqsID(dofIdB2, 0, 0));
    ASSERT_EQ(19, dofManagerB.mapEqsID(dofIdB2, 0, 9));
};

TEST(Discrete, DoF_var2_byPoint)
{
    DofEquationIdTable dofManagerB;
    size_t dofIdB1 = dofManagerB.addVariableDoFs(0, 0, 10);
    size_t dofIdB2 = dofManagerB.addVariableDoFs(0, 0, 10);
    dofManagerB.construct(DofNumberingType::BY_POINT);
    ASSERT_EQ(2, dofManagerB.getNumberOfVariables());
    ASSERT_EQ(1, dofManagerB.getNumberOfMeshes());
    ASSERT_EQ(0, dofManagerB.mapEqsID(dofIdB1, 0, 0));
    ASSERT_EQ(18, dofManagerB.mapEqsID(dofIdB1, 0, 9));
    ASSERT_EQ(1, dofManagerB.mapEqsID(dofIdB2, 0, 0));
    ASSERT_EQ(19, dofManagerB.mapEqsID(dofIdB2, 0, 9));
}

TEST(Discrete, DoF_inactive)
{
    DofEquationIdTable *dofManager = new DofEquationIdTable();
    size_t varId = dofManager->addVariableDoFs(0, 0, 9);
    dofManager->deactivateDoFs(varId, 0, 8);
    dofManager->construct(DofNumberingType::BY_POINT);

    ASSERT_EQ(1, dofManager->getNumberOfVariables());
    ASSERT_EQ(8, dofManager->getTotalNumberOfActiveDoFs());
    ASSERT_EQ(8, dofManager->getTotalNumberOfActiveDoFsWithoutGhost());
    ASSERT_EQ(0, dofManager->mapEqsID(varId, 0, 0));
    ASSERT_EQ(-1, dofManager->mapEqsID(varId, 0, 8));

    delete dofManager;
}

TEST(Discrete, DoF_ghost_nodes)
{
    {
        DofEquationIdTable *dofManager = new DofEquationIdTable();
        int ghost_nodes[] = {5, 6, 7, 8};
        std::vector<size_t> vec_ghost_nodes(ghost_nodes, ghost_nodes+4);
        size_t varId = dofManager->addVariableDoFs(0, 0, 9);
        dofManager->setGhostPoints(0, vec_ghost_nodes);
        dofManager->construct(DofNumberingType::BY_POINT);

        ASSERT_EQ(1, dofManager->getNumberOfVariables());
        ASSERT_EQ(9, dofManager->getTotalNumberOfActiveDoFs());
        ASSERT_EQ(5, dofManager->getTotalNumberOfActiveDoFsWithoutGhost());
        ASSERT_EQ(0, dofManager->mapEqsID(varId, 0, 0));
        ASSERT_TRUE(dofManager->isGhostPoint(0, 8));
        ASSERT_EQ(8, dofManager->mapEqsID(varId, 0, 8));

        delete dofManager;
    }
    {
        DofEquationIdTable *dofManager = new DofEquationIdTable();
        int ghost_nodes[] = {0, 1, 2, 3, 4};
        std::vector<size_t> vec_ghost_nodes(ghost_nodes, ghost_nodes+5);
        size_t varId = dofManager->addVariableDoFs(0, 0, 9);
        dofManager->setGhostPoints(0, vec_ghost_nodes);
        dofManager->construct(DofNumberingType::BY_POINT);

        ASSERT_EQ(1, dofManager->getNumberOfVariables());
        ASSERT_EQ(9, dofManager->getTotalNumberOfActiveDoFs());
        ASSERT_EQ(4, dofManager->getTotalNumberOfActiveDoFsWithoutGhost());
        ASSERT_TRUE(dofManager->isGhostPoint(0, 0));
        ASSERT_EQ(8, dofManager->mapEqsID(varId, 0, 4));
        ASSERT_EQ(0, dofManager->mapEqsID(varId, 0, 5));
        delete dofManager;
    }
}
