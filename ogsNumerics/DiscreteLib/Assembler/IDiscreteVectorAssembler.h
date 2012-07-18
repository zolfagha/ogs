
#pragma once

#include "DiscreteLib/Core/IDiscreteVector.h"

namespace MeshLib 
{
class IMesh;
}

namespace DiscreteLib
{

class DofEquationIdTable;

/**
 * \brief Interface of discrete vector assembler classes
 */
template <class T>
class IDiscreteVectorAssembler
{
public:
    typedef IDiscreteVector<T> VectorType;

    virtual ~IDiscreteVectorAssembler() {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param vec Discrete vector
    virtual void assembly( const MeshLib::IMesh &msh, const DofEquationIdTable &dofManager, VectorType &vec) = 0;
};

}
