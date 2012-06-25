
#pragma once

namespace MeshLib
{
class IMesh;
}


namespace SolutionLib
{

/**
 * \brief Mesh based discrete IVBV problems
 */
class MeshBasedProblem
{
public:
	///
    explicit MeshBasedProblem(MeshLib::IMesh* msh) : _msh(msh)
	{};

    ///
    virtual ~MeshBasedProblem() {}

    /// set a mesh
    void setMesh(MeshLib::IMesh* msh) {_msh = msh;};

    /// get the mesh
    MeshLib::IMesh* getMesh() const {return _msh;};

private:
    MeshLib::IMesh* _msh;
};

}
