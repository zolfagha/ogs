
#pragma once

#include <vector>

#include "Base/CodingTools.h"

//#include "DiscreteLib/Core/DiscreteVector.h"
//#include "DiscreteLib/Core/DiscreteLinearEquation.h"

namespace DiscreteLib
{

class IDiscreteVectorBase;
class IDiscreteLinearEquation;

/**
 * \brief Data container for any discrete objects
 */
class DiscreteDataContainer
{
public:
    DiscreteDataContainer() {};
    virtual ~DiscreteDataContainer()
    {
        Base::releaseObjectsInStdVector(_vec_vectors);
        Base::releaseObjectsInStdVector(_vec_linear_sys);
    }

    size_t addVector(IDiscreteVectorBase* v)
    {
        _vec_vectors.push_back(v);
        return _vec_vectors.size()-1;
    }

    size_t getNumberOfVectors() const {return _vec_vectors.size();};

    IDiscreteVectorBase* getVector(size_t i)
    {
        return _vec_vectors[i];
    }

    size_t addLinearEquation(IDiscreteLinearEquation* eq)
    {
        _vec_linear_sys.push_back(eq);
        return _vec_linear_sys.size()-1;
    }

    size_t getNumberOfLinearEquations() const {return _vec_linear_sys.size();};

    IDiscreteLinearEquation* getLinearEquation(size_t i)
    {
        return _vec_linear_sys[i];
    }

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteDataContainer);

    std::vector<IDiscreteLinearEquation*> _vec_linear_sys;
    std::vector<IDiscreteVectorBase*> _vec_vectors;
};

} //end
