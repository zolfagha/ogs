
#pragma once

#include <vector>
#include <valarray>

#include "DiscreteLib/EquationId/DofEquationIdTable.h"
#include "DiscreteLib/Core/DiscreteVector.h"
#include "DiscreteLib/Assembler/IDiscreteVectorAssembler.h"

#include "IStencilWiseTransientResidualLocalAssembler.h"

namespace MeshLib
{
class IMesh;
}

namespace NumLib
{

class TimeStep;


/**
 * \brief Element-based discrete system assembler classes
 */
class StencilWiseTransientResidualAssembler : public DiscreteLib::IDiscreteVectorAssembler<double>
{
public:
	typedef DiscreteLib::IDiscreteVectorAssembler<double>::VectorType GlobalVectorType;
	typedef IStencilWiseTransientResidualLocalAssembler::LocalVectorType LocalVectorType;

    /// @param time
    /// @param u0
    /// @param u1
    /// @param a
	StencilWiseTransientResidualAssembler(const TimeStep* time, const std::vector<DiscreteLib::DiscreteVector<double>*>* u0, const std::vector<DiscreteLib::DiscreteVector<double>*>* u1, IStencilWiseTransientResidualLocalAssembler* a)
        : _transient_e_assembler(a), _timestep(time), _u0(u0), _u1(u1)
    { };


    /// Conduct the element by element assembly procedure
    ///
    /// @param msh 				Mesh
    /// @param dofManager 		Dof map manager
    /// @param list_dofId 		List of Dof IDs used in this problem
    /// @param r 				Residual
    virtual void assembly( const MeshLib::IMesh &msh, const DiscreteLib::DofEquationIdTable &dofManager, GlobalVectorType &r)
    {

    }

private:
    IStencilWiseTransientResidualLocalAssembler* _transient_e_assembler;
    const TimeStep* _timestep;
    const std::vector<DiscreteLib::DiscreteVector<double>*>* _u0;
    const std::vector<DiscreteLib::DiscreteVector<double>*>* _u1;
};


}
