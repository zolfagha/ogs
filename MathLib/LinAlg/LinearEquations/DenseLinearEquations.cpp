
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
    double *b = this->getRHS();
    double *x = this->getX();
    std::copy(b, b+this->getDimension(), x);
    solver.execute(x);
}


} //end namespace
