
#include "Ogs6FemData.h"

Ogs6FemData* Ogs6FemData::_obj = 0;

Ogs6FemData* Ogs6FemData::getInstance()
{
	if (_obj==0) _obj = new Ogs6FemData();
	return _obj;
}
