
#include "SimulationInfo.h"

#include <cstdio>
#include <iostream>

#include "logog/include/logog.hpp"

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "Configure.h"
#include "BuildInfo.h"

namespace ogs6
{

void SimulationInfo::outputHeader ( void )
{
//    char buf[128];
//
//    const int len = 47;

    LOGOG_COUT << std::endl;
    LOGOG_COUT << _LG("          ###################################################\n");
    LOGOG_COUT << _LG("          ##                                               ##\n");
    LOGOG_COUT << _LG("          ##              OpenGeoSys-Project 6             ##\n");
#ifdef USE_LIS
    LOGOG_COUT << _LG("          ## ") << _LG(BaseLib::bothPadding("powered by LIS",45)) << _LG(" ##\n");
#endif
    LOGOG_COUT << _LG("          ##                                               ##\n");
    LOGOG_COUT << _LG("          ##   Contributors                                ##\n");
    LOGOG_COUT << _LG("          ##   * Helmholtz Centre for Environmental        ##\n");
    LOGOG_COUT << _LG("          ##     Research - UFZ                            ##\n");
    LOGOG_COUT << _LG("          ##   * TU Dresden                                ##\n");
    LOGOG_COUT << _LG("          ##   * University of Kiel                        ##\n");
    LOGOG_COUT << _LG("          ##   * University of Edinburgh                   ##\n");
    LOGOG_COUT << _LG("          ##   * University of Tuebingen (ZAG)             ##\n");
    LOGOG_COUT << _LG("          ##   * Federal Institute for Geosciences         ##\n");
    LOGOG_COUT << _LG("          ##     and Natural Resources (BGR)               ##\n");
    LOGOG_COUT << _LG("          ##   * Helmholtz Centre Potsdam GFZ              ##\n");
    LOGOG_COUT << _LG("          ##     German Research Centre for Geosciences    ##\n");
    LOGOG_COUT << _LG("          ##                                               ##\n");
    LOGOG_COUT << _LG("          ##   Program version                             ##\n");
    LOGOG_COUT << _LG("          ##   * Version: ") << _LG(BaseLib::rightPadding(OGS_VERSION, 32)) << " ##\n";
    LOGOG_COUT << _LG("          ##   * Date   : ") << _LG(BaseLib::rightPadding(OGS_DATE, 32)) << " ##\n";
    LOGOG_COUT << _LG("          ##   * Rev.   : ") << _LG(BaseLib::rightPadding(" ", 32)) << " ##\n";
    LOGOG_COUT << _LG("          ##     ") << _LG(BaseLib::rightPadding(GIT_COMMIT_INFO, 41)) << " ##\n";
    LOGOG_COUT << _LG("          ##                                               ##\n");
    LOGOG_COUT << _LG("          ###################################################\n");
    LOGOG_COUT << _LG("\n\n");
}

SimulationInfo::SimulationInfo(const std::string &project_path)
{
    this->setProjectPath(project_path);
}


void SimulationInfo::setProjectPath(const std::string& path)
{
    _project_path = path;
    _project_dir = BaseLib::getFileDirecotryPath(path);
    _project_name = BaseLib::getFileBaseName(path);
};

} //end
