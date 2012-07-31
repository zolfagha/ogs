
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
        if (hasVariable(itr->first)) {
            OutputVariableInfo &var = itr->second;
            if (var.object_type==OutputObjectType::Node) {
                node_values.push_back(VtuWriter::PointData(var.name, var.value));
            } else if (var.object_type==OutputObjectType::Element) {
                ele_values.push_back(VtuWriter::CellData(var.name, var.value));
            }
        }
    }
    
    // write VTU file
    std::string vtk_file_name = getOutputDir()+"\\"+getOutputBaseName() + "_";
    vtk_file_name += BaseLib::number2str<size_t>(current_time.getTimeStepCount()) + ".vtu";
    VtuWriter vtuWriter(false);
    vtuWriter.write(vtk_file_name, *getMesh(), node_values, ele_values);
    
    // update PVD file
    if (_pvd==NULL) {
        _pvd = new PVDWriter();
        std::string pvd_name = getOutputDir()+"\\"+getOutputBaseName() + ".pvd";
        _pvd->initialize(pvd_name);
    }
    
    _pvd->update(current_time.getTime(), vtk_file_name);
}
