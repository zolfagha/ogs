
#pragma once

#include <string>
#include <vector>

#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"

namespace NumLib
{
class ITXFunction;
}

namespace SolutionLib
{

/**
 * \brief Variable
 */
class FemVariable
{
public:
	FemVariable(size_t id, const std::string &name) : _id(id), _name(name)
	{

	}

	//----------------------------------------------------------------------
	size_t getID() const {return _id;};
	const std::string& getName() const { return _name;}


	//----------------------------------------------------------------------
    void setIC(NumLib::ITXFunction* ic) { _f_ic = ic; };
    NumLib::ITXFunction* getIC() const { return _f_ic; };


	//----------------------------------------------------------------------
    void addDirichletBC(FemLib::FemDirichletBC* bc)
    {
        _map_bc1.push_back(bc);
    }
    size_t getNumberOfDirichletBC() const {return _map_bc1.size();};
    FemLib::FemDirichletBC* getDirichletBC(size_t bc_id) const
    {
        return _map_bc1[bc_id];
    };


	//----------------------------------------------------------------------
    void addNeumannBC(FemLib::IFemNeumannBC* bc2)
    {
        _map_bc2.push_back(bc2);
    }
    size_t getNumberOfNeumannBC() const {return _map_bc2.size();};
    FemLib::IFemNeumannBC* getNeumannBC(size_t bc_id) const
    {
        return _map_bc2[bc_id];
    };

private:
    size_t _id;
    std::string _name;
    NumLib::ITXFunction* _f_ic;
    std::vector<FemLib::FemDirichletBC*> _map_bc1;
    std::vector<FemLib::IFemNeumannBC*> _map_bc2;
};

} //end

