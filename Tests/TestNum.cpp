
#include <gtest/gtest.h>

#include <vector>

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"

#include "NumLib/Discrete/DiscretizedEQS.h"
#include "NumLib/Discrete/DoF.h"

#include "NumLib/Coupling/Clock.h"
#include "NumLib/Coupling/TransientSystems.h"

using namespace GeoLib;
using namespace MeshLib;
using namespace NumLib;

TEST(Num, testDis1)
{    
    //mesh
    IMesh* msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    //define dof
    DofMapManager dofManager;
    size_t dofId = dofManager.addDoF(msh->getNumberOfNodes());
    // construct discrete eqs
    DiscretizedEQS eqs;

    //
    const DofMap *dofMap = dofManager.getDofMap(dofId);
    for (size_t i=0; i<msh->getNumberOfElements(); i++) {
        const IElement* e = msh->getElemenet(i);
        std::vector<size_t> ele_node_ids, local_dofmap;
        e->getNodeIDList(ele_node_ids);
        dofMap->getListOfEqsID(ele_node_ids, local_dofmap);
        MathLib::Matrix<double> localK;
        eqs.add(local_dofmap, localK);
    }

}

TEST(Num, testTransient1)
{
    std::vector<ITransientSystem*> vec_systems;

    Clock clock;
    for (size_t i=0; i<vec_systems.size(); i++)
        clock.addTransientSystem(vec_systems[i]);

    TimeStep t0, t_end;
    clock.setBeginning(t0);

    // start clock
    clock.moveForwardUntill(t_end); 
}

