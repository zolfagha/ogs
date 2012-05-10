
#pragma once

#include <vector>
#include "DiscreteLib/Core/DiscreteVector.h"
#include "DiscreteLib/Core/DataType.h"
#include "DiscreteLib/EquationId/DofEquationIdTable.h"

namespace DiscreteLib
{

/// create a subset of vector u corresponding to the given vector index
void getLocalVector2(const DofEquationIdTable &dofManager, const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<DiscreteVector<double>*> &list_multiple_u, std::vector<double> &local_u);

/// create a subset of vector u corresponding to the given vector index
void getLocalVector(const DofEquationIdTable &dofManager, const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<DiscreteVector<double>*> &list_multiple_u, LocalVector &local_u);

} //end

