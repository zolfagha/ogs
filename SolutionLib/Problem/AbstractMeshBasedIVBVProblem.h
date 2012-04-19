
#pragma once

#include <vector>
#include <map>

#include "Base/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"

#include "IVBVProblem.h"

namespace SolutionLib
{

/**
 * \brief Mesh based discrete IVBV problems
 */
class AbstractMeshBasedDiscreteIVBVProblem : public IVBVProblem
{
public:
	///
    AbstractMeshBasedDiscreteIVBVProblem(MeshLib::IMesh &msh) : _msh(&msh), _tim(0) {};

    ///
    virtual ~AbstractMeshBasedDiscreteIVBVProblem()
    {
        Base::releaseObject(_tim);
    }

    /// set a mesh
    void setMesh(MeshLib::IMesh &msh);

    /// get the mesh
    MeshLib::IMesh* getMesh() const {return _msh;};

    /// set a time stepping function
    void setTimeSteppingFunction(NumLib::ITimeStepFunction &f)
    {
        _tim = f.clone();
    }

    /// get the time stepping function
    NumLib::ITimeStepFunction* getTimeSteppingFunction() const {return _tim;};

private:
    MeshLib::IMesh* _msh;
    NumLib::ITimeStepFunction* _tim;
};

}
