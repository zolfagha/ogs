
#include "FeObjectCachePerFeType.h"

namespace FemLib
{

IFiniteElement* FeObjectCachePerFeType::getFeObject(FiniteElementType::type fe_type)
{
    IFiniteElement *fe = 0;
    if (_mapFeObj.count(fe_type)==0) {
        fe = FemElementFactory::create(fe_type, *_msh);
        _mapFeObj[fe_type] = fe;
    } else {
        fe = _mapFeObj[fe_type];
    }
    return fe;
}

}
