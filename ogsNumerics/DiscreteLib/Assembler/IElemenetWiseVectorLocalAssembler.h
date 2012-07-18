
#pragma once

#include <valarray>

namespace MeshLib
{
class IElement;
}

namespace DiscreteLib
{

/**
 * \brief Interface of all element local assembler classes
 */
template <class T>
class IElemenetWiseVectorLocalAssembler
{
public:
    typedef std::valarray<T> LocalVectorType;

    virtual ~IElemenetWiseVectorLocalAssembler() {};

    /// assemble a local linear equation for the given element
    virtual void assembly(const MeshLib::IElement &e, LocalVectorType &eqs) = 0;
};

} //end
