
#pragma once

#include "NumLib/Discrete/DiscretizedEQS.h"

namespace NumLib
{

class CouplingEquations
{
private:

public:
    size_t addSystem(size_t var_id, std::vector<int> list_submatrix, std::vector<int> list_subRHS);

};

void test()
{
/*
element e;
fe_u, fe_p;

for each gp
  MAT.eval(gp);
  Kuu, Cup;
  Kpp, Cpu;
  A(dof(u,u))+=Kuu;
  A(dof(u,p))+=Cup;

// eqs: GW
primary:   p Kpp
secondary: u Kpu, T Kpt
F1
// eqs: M
primary:   u Kuu
secondary: p Kup, T Kpt
F2
// eqs: T
primary:   T Ktt
secondary: p Ktp, u Ktu
F3

// coupling equation
K11 u1 + K12 u2 = F1
K21 u1 + K22 u2 = F2

// monolithic
A = |K11 K12|
    |K21 K22|
// partitioned
[K11] u1 = F1 - K12 u2
[K22] u2 = F2 - K21 u1

*/
}


}
