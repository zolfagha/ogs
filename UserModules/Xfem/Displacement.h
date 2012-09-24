/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Displacement.h
 *
 * Created on 2012-09-20 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include "BaseLib/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "FemLib/Function/FemIntegrationPointFunction.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXVectorFunctionAsColumnData.h"

#include "ProcessLib/AbstractTimeIndependentProcess.h"

namespace xfem
{

/**
 *
 */
template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER>
class FunctionDisplacement
: public ProcessLib::AbstractTimeIndependentProcess
{
public:
    enum Out { Displacement=0 };

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef T_LINEAR_SOLVER MyLinearSolver;

    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
    typedef typename NumLib::TXVectorFunctionAsColumnData<MyNodalFunctionVector> NodalPointScalarWrapper;

    ///
    FunctionDisplacement()
    : ProcessLib::AbstractTimeIndependentProcess("XFEM_DEFORMATION", 0, 1),
      _feObjects(0), _displacement(0)//, _strain(0), _stress(0)
    {
        // set default parameter name
        this->setOutputParameterName(Displacement, "Displacement");
    };

    ///
    virtual ~FunctionDisplacement()
    {
        BaseLib::releaseObject(_feObjects, _displacement /*, _strain, _stress*/);
    }

    /// initialize this process
    virtual bool initialize(const BaseLib::Options &option);

    /// finalize but nothing to do here
    virtual void finalize() {};

    ///
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    void accept(const NumLib::TimeStep &/*time*/);

protected:
    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual void output(const NumLib::TimeStep &time);

private:
    //void calculateStressStrain();

    DISALLOW_COPY_AND_ASSIGN(FunctionDisplacement);

private:
    FemLib::LagrangianFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
    MyNodalFunctionVector* _displacement;
    std::vector<NodalPointScalarWrapper*> _vec_u_components;
//    FemLib::FEMIntegrationPointFunctionVector* _strain;
//    FemLib::FEMIntegrationPointFunctionVector* _stress;
    MeshLib::IMesh* _msh;
    MyDiscreteSystem* _dis;

};

}

#include "Displacement.hpp"

