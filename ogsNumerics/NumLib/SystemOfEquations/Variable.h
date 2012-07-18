
#pragma once

#include <string>

namespace NumLib
{

struct Variable
{
    size_t id;
    size_t n_dof;
    std::string name;

    Variable(size_t i, size_t n, std::string str) : id(i), n_dof(n), name(str) {};
};


} //end

