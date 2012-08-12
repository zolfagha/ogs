/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestNum.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

#include <vector>

#include "BaseLib/CodingTools.h"
#include "BaseLib/BidirectionalMap.h"

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/SparseLinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"

#include "GeoLib/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"


#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteLinearEquation.h"
#include "DiscreteLib/EquationId/DofEquationIdTable.h"
#include "DiscreteLib/Assembler/IDiscreteLinearEquationAssembler.h"
#include "DiscreteLib/Assembler/IElemenetWiseLinearEquationLocalAssembler.h"
#include "DiscreteLib/Utils/SparsityBuilder.h"
#ifdef USE_MPI
#include "DiscreteLib/ogs5/par_ddc_group.h"
#endif
#include "NumLib/Function/Function.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"

//#include "NumLib/Output/Output.h"

#include "TestUtil.h"

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace NumLib;
using namespace DiscreteLib;



TEST(Num, TimeStepping1)
{
    // define problems and solution strategy
    class TestProblem : public ITransientSystem
    {
    private:
        double _t0, _tn, _dt, _v;
    public:
        TestProblem(double t0, double tn, double dt) : _t0(t0), _tn(tn), _dt(dt), _v(.0) {};
        virtual ~TestProblem() {};
        int solveTimeStep(const TimeStep &/*time*/) { _v+=1.0; return 0; }
        double suggestNext(const TimeStep &time_current)
        {
            if (time_current.getTime() < _t0) return _t0;
            else return (time_current.getTime()+_dt <= _tn) ? time_current.getTime()+_dt : time_current.getTime();
        }
        bool isAwake(const TimeStep &time)  { return (time.getTime()>=_t0 && time.getTime()<=_tn); }
        void accept(const TimeStep &/*time*/) {};
        double getValue() {return _v;};
    };
    TestProblem problem(.0, 100., 10.);

    // pass it to discretization systems
    TimeSteppingController timeStepping;
    timeStepping.setTransientSystem(problem);
    // start time stepping
    timeStepping.setBeginning(.0);
    size_t n_steps = timeStepping.solve(100.);
    //check
    ASSERT_EQ(10, n_steps);
    ASSERT_EQ(10., problem.getValue());
}


