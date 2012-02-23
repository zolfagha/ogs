
#include "DenseLinearEquations.h"

#include <algorithm>

#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"


namespace MathLib
{

void DenseLinearEquations::setOption(const Base::Options &option)
{

}

void DenseLinearEquations::solve()
{
    MathLib::GaussAlgorithm solver(*this->getA());
    DenseLinearEquationsBase::VectorType *b = this->getRHSAsStdVec();
    DenseLinearEquationsBase::VectorType *x = this->getXAsStdVec();
    std::copy(b->begin(), b->end(), x->begin());
    solver.execute(&(*x)[0]);
}


} //end namespace
