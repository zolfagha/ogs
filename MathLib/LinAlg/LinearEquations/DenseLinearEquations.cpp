
#include "DenseLinearEquations.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"

namespace MathLib
{

void DenseLinearEquations::setOption(const Base::Options &option)
{

}

void DenseLinearEquations::solve()
{
    MathLib::GaussAlgorithm solver(*this->getA());
    solver.execute(this->getX());
}


} //end namespace
