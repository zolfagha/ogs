/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file XFEM_EXAMPLE_CRACK1.h
 *
 * Created on 2012-09-20 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include "BaseLib/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "FemLib/Function/FemNodalFunction.h"
//#include "FemLib/Function/FemIntegrationPointFunction.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
//#include "NumLib/Function/TXVectorFunctionAsColumnData.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "ProcessLib/AbstractTimeIndependentProcess.h"

namespace xfem
{

/**
 *
 */
//template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER>
class FunctionXFEM_EXAMPLE_CRACK1
: public ProcessLib::AbstractTimeIndependentProcess
{
public:
    //enum Out { Displacement=0 };

    typedef DiscreteLib::DiscreteSystem MyDiscreteSystem;
//    //typedef T_LINEAR_SOLVER MyLinearSolver;
//
//    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
//    typedef typename NumLib::TXVectorFunctionAsColumnData<MyNodalFunctionVector> NodalPointScalarWrapper;

    ///
    FunctionXFEM_EXAMPLE_CRACK1()
    : ProcessLib::AbstractTimeIndependentProcess("XFEM_DEFORMATION", 0, 0),
      _feObjects(0) //, _displacement(0)//, _strain(0), _stress(0)
    {
        // set default parameter name
        //this->setOutputParameterName(Displacement, "Displacement");
    };

    ///
    virtual ~FunctionXFEM_EXAMPLE_CRACK1()
    {
        BaseLib::releaseObject(_feObjects /*, _displacement , _strain, _stress*/);
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

    DISALLOW_COPY_AND_ASSIGN(FunctionXFEM_EXAMPLE_CRACK1);

private:
    FemLib::LagrangianFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
    MyNodalFunctionVector* _displacement;
    MyNodalFunctionVector* _exact_displacement;
//    std::vector<NodalPointScalarWrapper*> _vec_u_components;
//    FemLib::FEMIntegrationPointFunctionVector* _strain;
//    FemLib::FEMIntegrationPointFunctionVector* _stress;
    MeshLib::IMesh* _msh;
    MyDiscreteSystem* _dis;

};

}

//#include "Displacement.hpp"

