
#pragma once

#include <vector>

#include "Base/CodingTools.h"

#include "DiscreteLib/Core/IDiscreteResource.h"
#include "DiscreteLib/Core/IDiscreteLinearEquation.h"

namespace DiscreteLib
{

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

    size_t addVector(IDiscreteResource* v)
    {
        _vec_vectors.push_back(v);
        v->setObjectID(_vec_vectors.size()-1);
        return _vec_vectors.size()-1;
    }

    void eraseVector(IDiscreteResource* v)
    {
    	const size_t i = v->getObjectID();
    	if (_vec_vectors.size() > i) {
    		_vec_vectors[i] = 0;
    	}
    }

    size_t getNumberOfVectors() const {return _vec_vectors.size();};

    IDiscreteResource* getVector(size_t i)
    {
        return _vec_vectors[i];
    }

    size_t addLinearEquation(IDiscreteLinearEquation* eq)
    {
        _vec_linear_sys.push_back(eq);
        eq->setObjectID(_vec_linear_sys.size()-1);
        return _vec_linear_sys.size()-1;
    }

    void eraseLinearEquation(IDiscreteLinearEquation* eq)
    {
    	const size_t i = eq->getObjectID();
    	if (_vec_linear_sys.size() > i) {
    		_vec_linear_sys[i] = 0;
    	}
    }

    size_t getNumberOfLinearEquations() const {return _vec_linear_sys.size();};

    IDiscreteLinearEquation* getLinearEquation(size_t i)
    {
        return _vec_linear_sys[i];
    }

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteDataContainer);

    std::vector<IDiscreteLinearEquation*> _vec_linear_sys;
    std::vector<IDiscreteResource*> _vec_vectors;
};

} //end
