
#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"
//#include "IVBVProblem.h"
#include "TimeSteppingProblem.h"
#include "FemVariable.h"

namespace SolutionLib
{

class AbstractFemIVBVProblem : public TimeSteppingProblem
{
public:
    ///
    explicit AbstractFemIVBVProblem(    DiscreteLib::DiscreteSystem* dis)
        : _discrete_system(dis)
    {
    }

    ///
    virtual ~AbstractFemIVBVProblem()
    {
        BaseLib::releaseObjectsInStdVector(_variables);
    }

    /// get this discrete system
    DiscreteLib::DiscreteSystem* getDiscreteSystem() {return _discrete_system;};

    /// get the mesh
    MeshLib::IMesh* getMesh() {return _discrete_system->getMesh();};

    /// create FE approximation field
    FemVariable* addVariable(const std::string name)
    {
        _variables.push_back(new FemVariable(_variables.size(), name));
        return _variables.back();
    }

    /// get a variable
    FemVariable* getVariable(size_t i) const { return _variables[i]; }

    /// get the number of variables
    size_t getNumberOfVariables() const { return _variables.size(); }

private:
    DISALLOW_COPY_AND_ASSIGN(AbstractFemIVBVProblem);

private:
    DiscreteLib::DiscreteSystem* _discrete_system;
    std::vector<FemVariable*> _variables;
};


} //end
