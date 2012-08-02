
#include "Ogs5FemData.h"


#include "BaseLib/CodingTools.h"

namespace ogs5
{

Ogs5FemData::~Ogs5FemData()
{
    BaseLib::releaseObjectsInStdVector(pcs_vector);
    BaseLib::releaseObjectsInStdVector(mfp_vector);
    BaseLib::releaseObjectsInStdVector(msp_vector);
    BaseLib::releaseObjectsInStdVector(mmp_vector);
    BaseLib::releaseObjectsInStdVector(cp_vector);
    BaseLib::releaseObjectsInStdVector(bc_vector);
    BaseLib::releaseObjectsInStdVector(st_vector);
    BaseLib::releaseObjectsInStdVector(ic_vector);
    BaseLib::releaseObjectsInStdVector(out_vector);
    BaseLib::releaseObjectsInStdVector(time_vector);
    BaseLib::releaseObjectsInStdVector(num_vector);
    //geo, msh objects are passed
}

}
