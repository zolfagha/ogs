
#pragma once

#include <vector>

#include "Base/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"
#include "NumLib/Coupling/ICoupledProblem.h"
#include "NumLib/TimeStepping/ITransientProblem.h"
#include "NumLib/TimeStepping/TimeStepping.h"

namespace NumLib
{


class TransientFemProblem: public ITransientProblem
{
private:
    MeshLib::IMesh* _msh;
    std::vector<FemLib::FemDirichletBC<double>*> vec_bc1;
    std::vector<FemLib::FemNeumannBC<double, double>*> vec_bc2;
    FixedTimeStepping *TIM;
    FemLib::FemNodalFunctionScalar* _u0;
    FemLib::FemNodalFunctionScalar* _u_t0;
    FemLib::FemNodalFunctionScalar* _u_t1;

public:	
    TransientFemProblem(FemLib::FemNodalFunctionScalar& u0) 
    {
        _msh = (MeshLib::IMesh*) u0.getMesh();
    };

    virtual ~TransientFemProblem() 
    {
        Base::releaseObjectsInStdVector(vec_bc1);
        Base::releaseObjectsInStdVector(vec_bc2);
    };

    void addDirichletBC(FemLib::FemDirichletBC<double>& bc1)
    {
        vec_bc1.push_back(&bc1);
    }

    void addNeumannBC(FemLib::FemNeumannBC<double, double>& bc2)
    {
        vec_bc2.push_back(&bc2);
    }



    double suggestNext(const TimeStep &time_current)
    {
        return TIM->next();
    }

    bool isAwake(const TimeStep &time)
    {
        return (TIM->next()==time.getTime());
    }

    int solveTimeStep(const TimeStep &time)
    {
    }

};

}
