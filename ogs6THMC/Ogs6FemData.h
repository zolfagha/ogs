/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs6FemData.h
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <string>

#include "BaseLib/OrderedMap.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "MaterialLib/Fluid.h"
#include "MaterialLib/Compound.h"
#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"
#include "GeoLib/GEOObjects.h"
#include "ProcessLib/Process.h"
#include "OutputController.h"

/**
 * \brief Fem data storage
 *
 * This class follows singleton pattern.
 */
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
    std::vector<MaterialLib::Compound*> list_compound;
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
    std::vector<DiscreteLib::IDiscreteSystem*> list_dis_sys;
    //
    ogs6::OutputController outController;
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
