
#pragma once

#include <vector>
#include <string>

#include "BaseLib/OrderedMap.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "MaterialLib/Fluid.h"
#include "SolutionLib/FemProblem/FemDirichletBC.h"
#include "SolutionLib/FemProblem/FemNeumannBC.h"
#include "GeoLib/GEOObjects.h"
#include "ProcessLib/Process.h"

class Ogs6FemData
{
public:
    static Ogs6FemData* getInstance();
private:
    static Ogs6FemData* _obj;

private:
    Ogs6FemData(): geo(NULL) {};

public:
    //material data
    std::vector<MaterialLib::PorousMedia*> list_pm;
    std::vector<MaterialLib::Solid*> list_solid;
    std::vector<MaterialLib::Fluid*> list_fluid;
    //geometric data
    std::string geo_unique_name;
    GeoLib::GEOObjects* geo;
    //mesh data
    //std::vector<MeshLib::IMesh*> list_mesh;
    //time group data
    std::vector<NumLib::ITimeStepFunction*> list_tim;
    //process
    BaseLib::OrderedMap<std::string, ProcessLib::Process*> list_pcs;
    //discrete system
    std::vector<DiscreteLib::DiscreteSystem*> list_dis_sys;
    //
    std::string project_name;
    std::string project_dir;
    std::string output_dir;

    ~Ogs6FemData()
    {
        BaseLib::releaseObject(geo);
        BaseLib::releaseObjectsInStdVector(list_pm);
        BaseLib::releaseObjectsInStdVector(list_solid);
        BaseLib::releaseObjectsInStdVector(list_fluid);
        //BaseLib::releaseObjectsInStdVector(list_mesh);
        BaseLib::releaseObjectsInStdVector(list_tim);
        BaseLib::releaseObjectsInStdVector(list_dis_sys);
        BaseLib::releaseObjectsInStdMap(list_pcs);
    }
};
