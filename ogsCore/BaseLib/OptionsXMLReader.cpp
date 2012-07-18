
#include "OptionsXMLReader.h"

#include "logog/include/logog.hpp"
#include "tinyxml2.h"

namespace BaseLib
{

void addXMLtoOptions(tinyxml2::XMLElement* e_root, BaseLib::Options &properties)
{
    BaseLib::Options* optRoot = properties.addSubGroup(e_root->GetText());

    for (tinyxml2::XMLElement* e=e_root->FirstChildElement(); e!=0; e=e_root->NextSiblingElement())
    {
        if (e->FirstChildElement()!=0) {
            addXMLtoOptions(e, properties);
        } else {
            properties.addOption(e->GetText(), e->Value());
        }
    }
}

bool addXMLtoOptions(const std::string &xml_file, BaseLib::Options &properties)
{
    tinyxml2::XMLDocument doc;
    int ret = doc.LoadFile(xml_file.c_str());
    if (ret != tinyxml2::XML_NO_ERROR) {
        LOGOG_CERR << "Error in reading a XML file " << xml_file << " with error code " << ret << std::endl;
        return false;
    }

    // read data
    addXMLtoOptions(doc.RootElement(), properties);


    return true;
}

} //end
