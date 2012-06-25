
#pragma once

#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/IMesh.h"

#include "IFemElement.h"

namespace FemLib
{

/**
 * \brief Base class for implementation of finite element classes
 */
template <FiniteElementType::type T_FETYPE, size_t N_VARIABLES>
class TemplateFeBase : public IFiniteElement
{
public:
    TemplateFeBase(MeshLib::IMesh &msh) : _msh(&msh), _ele(0) {};
    virtual ~TemplateFeBase() {};

    void setMesh(MeshLib::IMesh &msh) {_msh = &msh;};

    const MeshLib::IMesh* getMesh() const {return _msh;};

    /// return mesh element
    MeshLib::IElement* getElement() const {return _ele;};

    /// return the number of variables
    size_t getNumberOfVariables() const { return N_VARIABLES; };

    /// return finite element type
    FiniteElementType::type getFeType() const { return T_FETYPE; };

protected:
    void setElement(MeshLib::IElement &e) {_ele = &e;};

private:
    MeshLib::IMesh* _msh;
    MeshLib::IElement* _ele;
};


}

