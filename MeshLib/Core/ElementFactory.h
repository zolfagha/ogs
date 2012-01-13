#pragma once

#include <cmath>
#include <vector>
#include "Element.h"

//------------------------------------------------------------------------
namespace MeshLib
{
class ElemenetFactory
{
public:
  static IElement* createNewElement(const ElementType::type t)
  {
    if (t == ElementType::LINE)	return new Line();
    else if (t == ElementType::QUAD) return new Quadrirateral();
    else if (t == ElementType::HEXAHEDRON) return new Hexahedron();
    else if (t == ElementType::TRIANGLE) return new Triangle();
    else if (t == ElementType::TETRAHEDRON) return new Tetrahedron();
    else if (t == ElementType::PRISM) return new Prism();
    else if (t == ElementType::PYRAMID) return new Pyramid();
    return NULL;
  };
};
} // end namespace

