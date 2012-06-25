
#pragma once

namespace FemLib
{
/// Finite element type
struct FiniteElementType
{
    enum type {
        LINE2,
        LINE3,
        TRI3,
        TRI3CONST,
        TRI6,
        QUAD4,
        QUAD8,
        QUAD9,
        INVALID
    };
};

}
