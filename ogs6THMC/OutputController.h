/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OutputController.h
 *
 * Created on 2012-07-31 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/OrderedMap.h"
#include "GeoLib/GEOObjects.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"

namespace ogs6
{

/**
 * \brief Result output controller
 */
class OutputController
{
public:
    void initialize(const BaseLib::Options &option, const std::string &output_dir, const std::string &project_name, std::vector<DiscreteLib::DiscreteSystem*> &list_dis_sys, GeoLib::GEOObjects &geo, const std::string &geo_unique_name)
    {
        //Ogs6FemData* femData = Ogs6FemData::getInstance();
        OutputBuilder outBuilder;
        OutputTimingBuilder outTimBuilder;

        const BaseLib::Options* opOutList = option.getSubGroup("OutputList");
        for (const BaseLib::Options* op = opOutList->getFirstSubGroup("Output"); op!=0; op = opOutList->getNextSubGroup())
        {
            std::string data_type = op->getOption("DataType");
            IOutput* out = outBuilder.create(data_type);
            out->setOutputPath(output_dir, project_name);

            size_t msh_id = 0; //TODO
            out->setMesh(list_dis_sys[msh_id]->getMesh());

            std::string time_type = op->getOption("TimeType");
            size_t out_steps = op->getOption<size_t>("TimeSteps");
            out->setOutputTiming(outTimBuilder.create(time_type, out_steps));

            const std::vector<std::string>* nod_var_name = op->getOptionAsArray<std::string>("NodalVariables");
            out->addNodalVariable(*nod_var_name);
            const std::vector<std::string>* ele_var_name = op->getOptionAsArray<std::string>("ElementalVariables");
            out->addElementalVariable(*ele_var_name);

            std::string geo_type = op->getOption("GeometryType");
            std::string geo_name = op->getOption("GeometryName");
            const GeoLib::GeoObject* geo_obj = geo.searchGeoByName(geo_unique_name, geo_type, geo_name);
            out->setGeometry(geo_obj);

            _list_output.push_back(out);
        }
    }

    bool isActive(const NumLib::TimeStep &time) const
    {
        bool doOutput = false;
        for (size_t i=0; i<_list_output.size(); i++) {
            if (_list_output[i]->isActive(time)) {
                doOutput = true;
                break;
            }
        }
        return doOutput;
    }

    void outputData(const NumLib::TimeStep &time)
    {
        if (!isActive(time)) return;

        for (size_t i=0; i<_list_output.size(); i++) {
            if (_list_output[i]->isActive(time)) {
                _list_output[i]->write(time, _map_name_var);
            }
        }
    }

    void setOutput(const std::string &var_name, OutputVariableInfo &var)
    {
        _map_name_var.insert(var_name, var);
    }

private:
    std::vector<IOutput*> _list_output;
    BaseLib::OrderedMap<std::string, OutputVariableInfo> _map_name_var;
};

}
