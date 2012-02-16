
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

namespace NumLib
{

#if 0
class ElementAssembly
{
    void get(int ele, int &A, int &b);
};

class SubDomain
{
    int list_ele;
    int list_nod;
};


class IDomainDecomposedSystem
{
public:
    virtual void addSubDomain(int list_ele) = 0;

    virtual void setProblem() = 0;

    virtual void solve() = 0;
};

class NodeBasedDomainDecomposedSystem : public IDomainDecomposedSystem
{
public:
    void addSubDomain(int list_ele);

    void setProblem();

    void solve();
};


class ElementBasedDomainDecomposedSystem : public IDomainDecomposedSystem
{
public:
    void addSubDomain(int list_ele);

    void setProblem();

    void solve();
};
#endif

}
