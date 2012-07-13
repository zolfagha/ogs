
#pragma once

#include <string>

#include "logog/include/logog.hpp"
#include "tinyxml2.h"

#include "BaseLib/Options.h"

namespace ogs6
{

/*

xml
	ProcessData
		Process, type=gw, name=gw1, time_group_id=0
			file
		Process, type=element_velocity, name=ev
			file
	TimeData
		TimeGroup, id=0, start=0, end=30, unit=DAY, control=STEPS
	coupling
		P
			M, pcs=gw1
				out, name=p
			M, pcs=ev
				in, name=p
				out, name=v

 */

inline bool readSimulationProperties(const std::string &xml_file, BaseLib::Options &properties)
{
	tinyxml2::XMLDocument doc;
	int ret = doc.LoadFile(xml_file.c_str());
	if (ret != tinyxml2::XML_NO_ERROR) {
		LOGOG_CERR << "Error in reading a XML file " << xml_file << " with error code " << ret << std::endl;
		return false;
	}
	
	// quick check
	if (doc.FirstChildElement("ProcessData")==0) {
		LOGOG_CERR << "Error in reading a XML file: ProcessData tag is not found." << std::endl;
		return false;
	}
	if (doc.FirstChildElement("TimeData")==0) {
		LOGOG_CERR << "Error in reading a XML file: TimeData tag is not found." << std::endl;
		return false;
	}
	
	// read data
	tinyxml2::XMLElement* e;
	BaseLib::Options* optProcessData = properties.addSubGroup("ProcessData");
	for (e = doc.FirstChildElement("ProcessData")->FirstChildElement("Process"); e != 0; e->NextSibling()) {
		const char* pcs_type = e->Attribute("type");
		const char* pcs_name = e->Attribute("name");
		int tim_id = -1;
		e->QueryIntAttribute("time_group_id", &tim_id);
		
		
		BaseLib::Options* optProcess = optProcessData->addSubGroup(pcs_name);
		optProcess->addOption("type", pcs_type);
		optProcess->addOptionAsNum("time_group_id", tim_id);
	}

	BaseLib::Options* optTimeData = properties.addSubGroup("TimeData");
	for (e = doc.FirstChildElement("TimeData")->FirstChildElement("TimeGroup"); e != 0; e->NextSibling()) {
		int tim_id = e->IntAttribute("id");
		double tim_start = e->DoubleAttribute("start");
		double tim_end = e->DoubleAttribute("end");
		const char* tim_unit = e->Attribute("unit");
		const char* tim_control = e->Attribute("control");
		
		
		BaseLib::Options* optTimeGroup = optTimeData->addSubGroup("TimeData");
		optTimeGroup->addOptionAsNum("time_group_id", tim_id);
	}

	
	return true;
}
 
} 