
#include "SimulatorInfo.h"

#include <cstdio>
#include <iostream>
#include "logog/include/logog.hpp"
#include "Configure.h"
#include "BuildInfo.h"

namespace ogs6
{

//align string to center of the line
void DisplayCenterAlign(const char* buf, int len)
{
	//len -= 2;
	LOGOG_COUT << _LG("          ##");
	int pad_len = len - (int)strlen(buf);
	for (int i = 0; i < pad_len / 2; i++)
		LOGOG_COUT << _LG(" ");
	LOGOG_COUT << buf;
	for (int i = 0; i < pad_len - pad_len / 2; i++)
		LOGOG_COUT << _LG(" ");
	LOGOG_COUT << _LG("##\n");
}

void DisplayFill(const char* buf, int len)
{
	LOGOG_COUT << _LG("          ##");
	for (int i = 0; i < len; i++)
		LOGOG_COUT << _LG(buf);
	LOGOG_COUT << _LG("##\n");
}

void SimulatorInfo::output ( void )
{
	char buf[128];

	const int len = 47;

	LOGOG_COUT << std::endl;
	LOGOG_COUT << _LG("          ###################################################\n");
	LOGOG_COUT << _LG("          ##                                               ##\n");
	LOGOG_COUT << _LG("          ##              OpenGeoSys-Project 6             ##\n");
	LOGOG_COUT << _LG("          ##                                               ##\n");
	LOGOG_COUT << _LG("          ##  Helmholtz Center for Environmental Research  ##\n");
	LOGOG_COUT << _LG("          ##    UFZ Leipzig - Environmental Informatics    ##\n");
	LOGOG_COUT << _LG("          ##                  TU Dresden                   ##\n");
	LOGOG_COUT << _LG("          ##              University of Kiel               ##\n");
	LOGOG_COUT << _LG("          ##            University of Edinburgh            ##\n");
	LOGOG_COUT << _LG("          ##         University of Tuebingen (ZAG)         ##\n");
	LOGOG_COUT << _LG("          ##       Federal Institute for Geosciences       ##\n");
	LOGOG_COUT << _LG("          ##          and Natural Resources (BGR)          ##\n");
	LOGOG_COUT << _LG("          ##                                               ##\n");

	sprintf(buf, "Version %s  Date %s", OGS_VERSION, OGS_DATE);
	DisplayCenterAlign(buf, len);
	DisplayCenterAlign("Git revision:", len);
	DisplayCenterAlign(GIT_COMMIT_INFO, len);
	LOGOG_COUT << _LG("          ##                                               ##\n");
	LOGOG_COUT << _LG("          ###################################################\n");
	LOGOG_COUT << _LG("\n\n");
}

} //end
