
#pragma once

#include <vector>

#include "MathLib/MatrixVector.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/EquationId/DofEquationIdTable.h"

#include "IDiscreteVectorAssembler.h"
#include "IElemenetWiseVectorLocalAssembler.h"

namespace MeshLib 
{
class IMesh;
}

namespace DiscreteLib
{

/**
 * \brief Element-based discrete vector assembler classes
 */
template <class T>
class ElementWiseVectorAssembler : public IDiscreteVectorAssembler<T>
{
public:
	typedef IDiscreteVectorAssembler<T>::VectorType GlobalVectorType;

    ///
	explicit ElementWiseVectorAssembler(IElemenetWiseVectorLocalAssembler<T> &a) : _e_assembler(&a) {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param vec Discrete vector
    void assembly(const MeshLib::IMesh &msh, const DofEquationIdTable &dofManager, GlobalVectorType &globalVec)
    {
    	MathLib::Vector localVec;
        std::vector<size_t> ele_node_ids, ele_node_size_order;
        std::vector<long> local_dofmap_row;
        const size_t n_ele = msh.getNumberOfElements();

        std::vector<double> local_u_n;
        for (size_t i=0; i<n_ele; i++) {
            MeshLib::IElement *e = msh.getElemenet(i);
            // get dof map
            e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);
            e->getListOfNumberOfNodesForAllOrders(ele_node_size_order);
            dofManager.mapEqsID(msh.getID(), ele_node_ids, local_dofmap_column, local_dofmap_row); //TODO order
            // local assembly
            localVec.resize(local_dofmap_column.size(), .0);
            _e_assembler->assembly(*e, localVec);
            // update global
            globalVec.addSubvector(local_dofmap_row, &localVec[0]);
        }
    };

private:
    IElemenetWiseVectorLocalAssembler<T>* _e_assembler;
};


}
