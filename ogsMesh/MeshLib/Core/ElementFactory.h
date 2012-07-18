#pragma once

#include "Element.h"

namespace MeshLib
{
class ElemenetFactory
{
public:
  static IElement* createNewElement(const ElementShape::type t)
  {
    if (t == ElementShape::LINE)    return new Line();
    else if (t == ElementShape::QUAD) return new Quadrirateral();
    else if (t == ElementShape::HEXAHEDRON) return new Hexahedron();
    else if (t == ElementShape::TRIANGLE) return new Triangle();
    else if (t == ElementShape::TETRAHEDRON) return new Tetrahedron();
    else if (t == ElementShape::PRISM) return new Prism();
    else if (t == ElementShape::PYRAMID) return new Pyramid();
    return NULL;
  };
};

} // end namespace

