/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PVDOutput.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "PVDOutput.h"

#include <vector>

#include "CommonIO/PVDWriter.h"
#include "CommonIO/VtuWriter.h"


void PVDOutput::write(  const NumLib::TimeStep &current_time, 
                        BaseLib::OrderedMap<std::string, OutputVariableInfo> &data)
{
    // prepare vtu data
    std::vector<VtuWriter::PointData> node_values;
    std::vector<VtuWriter::CellData> ele_values;
    
    for (BaseLib::OrderedMap<std::string, OutputVariableInfo>::iterator itr = data.begin(); itr!=data.end(); ++itr) {
        if ((itr->second.object_type == OutputVariableInfo::Node && hasNodalVariable(itr->first))
            || (itr->second.object_type == OutputVariableInfo::Element && hasElementalVariable(itr->first))) {
            OutputVariableInfo &var = itr->second;
            VtuWriter::AttributeInfo attr(var.name, var.nr_of_components, var.value);
            switch (var.data_type) {
                case OutputVariableInfo::Char: attr.data_type = VtuWriter::Char; break;
                case OutputVariableInfo::Int: attr.data_type = VtuWriter::Int; break;
                case OutputVariableInfo::Real: attr.data_type = VtuWriter::Real; break;
            }
            if (var.object_type==OutputVariableInfo::Node) {
                node_values.push_back(VtuWriter::PointData(var.name, attr));
            } else if (var.object_type==OutputVariableInfo::Element) {
                ele_values.push_back(VtuWriter::CellData(var.name, attr));
            }
        }
    }
    
    // write VTU file
    std::string vtk_file_name_relative = getOutputBaseName() + "_";
    vtk_file_name_relative += BaseLib::number2str<size_t>(current_time.getTimeStepCount()) + ".vtu";
    std::string vtk_file_name_absolute = getOutputDir();
    if (vtk_file_name_absolute.length()>0) vtk_file_name_absolute += "/";
    vtk_file_name_absolute += vtk_file_name_relative;
    VtuWriter vtuWriter(false);
    vtuWriter.write(vtk_file_name_absolute, *getMesh(), node_values, ele_values);
    
    // update PVD file
    if (_pvd==NULL) {
        _pvd = new PVDWriter();
        std::string pvd_name = getOutputDir();
        if (pvd_name.length()>0) pvd_name += "/";
        pvd_name += getOutputBaseName() + ".pvd";
        //std::string pvd_name = getOutputDir()+"\\"+getOutputBaseName() + ".pvd";
        _pvd->initialize(pvd_name);
    }
    
    _pvd->update(current_time.getTime(), vtk_file_name_relative);
}
