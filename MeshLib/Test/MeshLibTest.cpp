
#include "MeshLib/Core/Element.h"
#include "MeshLib/Core/UnstructuredMesh.h"
#include "MeshLib/Core/StructuredMesh.h"
#include "MeshLib//Core/CoordinateSystem.h"

using namespace MeshLib;

int main (int argc, char* argv[]) {

    Quadrirateral quad;

    MeshLib::UnstructuredMesh mesh;
    double len[3];
    double unit_len[3];
    MeshLib::StructuredMesh<MeshLib::ElementShape::QUAD> m(CoordinateSystemType::XY, GeoLib::Point(), len, unit_len);
    return 0;
}

