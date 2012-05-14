
#pragma once

#include <string>
#include <vector>

#include "FemDirichletBC.h"
#include "FemNeumannBC.h"

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
    void addDirichletBC(FemDirichletBC* bc)
    {
        _map_bc1.push_back(bc);
    }
    size_t getNumberOfDirichletBC() const {return _map_bc1.size();};
    FemDirichletBC* getDirichletBC(size_t bc_id) const
    {
        return _map_bc1[bc_id];
    };


	//----------------------------------------------------------------------
    void addNeumannBC(IFemNeumannBC* bc2)
    {
        _map_bc2.push_back(bc2);
    }
    size_t getNumberOfNeumannBC() const {return _map_bc2.size();};
    IFemNeumannBC* getNeumannBC(size_t bc_id) const
    {
        return _map_bc2[bc_id];
    };

private:
    size_t _id;
    std::string _name;
    NumLib::ITXFunction* _f_ic;
    std::vector<FemDirichletBC*> _map_bc1;
    std::vector<IFemNeumannBC*> _map_bc2;
};

} //end

