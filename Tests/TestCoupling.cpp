
#include <gtest/gtest.h>
#include <vector>

#include "NumLib/Coupling/TimeSteppingController.h"
#include "NumLib/Coupling/AsyncPartSolution.h"
#include "NumLib/Coupling/PartitionedAlgorithm.h"
#include "NumLib/Coupling/ICoupledProblem.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/Coupling/PartitionedProblem.h"
#include "NumLib/Coupling/IterativePartitionedMethod.h"
#include "NumLib/Coupling/TransientPartitionedMethod.h"

using namespace MathLib;
using namespace NumLib;

// 2a + 2b + .3c = 6.9
// 3a + 5b + .2c = 13.6
// .5a + .3b + 3c = 10.1
// A. a=1, b=2, c=3
class WeakCouplingEQS1 : public TemplateSteadyMonolithicProblem<2,1>
{
public:
    enum Parameters { a=2, b = 0, c = 1 };
    int solve()
    {
        double vb = _vec_parameters[WeakCouplingEQS1::b]->eval(.0);
        double vc = _vec_parameters[WeakCouplingEQS1::c]->eval(.0);
        double va = 1./2.*(6.9 - 2.*vb - 0.3*vc);
        if (_vec_parameters[WeakCouplingEQS1::a]!=0)
            delete _vec_parameters[WeakCouplingEQS1::a];
        _vec_parameters[WeakCouplingEQS1::a] = new MathLib::FunctionConstant<double,double>(va);
        return 0;
    }
};

class WeakCouplingEQS2 :  public TemplateSteadyMonolithicProblem<2,1>
{
public:
    enum Parameters { a = 0, b = 2, c = 1 };
    int solve()
    {
        double va = _vec_parameters[WeakCouplingEQS2::a]->eval(.0);
        double vc = _vec_parameters[WeakCouplingEQS2::c]->eval(.0);
        double vb = 1./5.*(13.6-3*va-0.2*vc);
        if (_vec_parameters[WeakCouplingEQS2::b]!=0)
            delete _vec_parameters[WeakCouplingEQS2::b];
        _vec_parameters[WeakCouplingEQS2::b] = new MathLib::FunctionConstant<double,double>(vb);
        return 0;
    }
};

class WeakCouplingEQS3 : public TemplateSteadyMonolithicProblem<2,1>
{
public:
    enum Parameters { a = 0, b = 1, c = 2 };
    int solve()
    {
        double va = _vec_parameters[a]->eval(.0);
        double vb = _vec_parameters[b]->eval(.0);
        double vc = 1./3.*(10.1-0.5*va-0.3*vb);
        if (_vec_parameters[c]!=0)
            delete _vec_parameters[c];
        _vec_parameters[c] = new MathLib::FunctionConstant<double,double>(vc);
        return 0;
    }
private:
    double va, vb, vc;
};

#if 0
// 2a + 2b + 3c = 15
// 3a + 5b + 2c = 19
// 5a + 3b + 3c = 20
// A. a=1, b=2, c=3
class StrongCouplingEQS1 : public IMonolithicProblem
{
public:
    enum InputParameters { b = 0, c = 1 };
    enum OutputParameters { a = 0 };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};
    StrongCouplingEQS1() {
        _vec_in_var.resize(getNumberOfInputVarameters(), 0);
        _vec_out_var.resize(getNumberOfOutputParameters(), 0);
    }
    int solve()
    {
        double vb = _vec_in_var[EQS1::b]->eval(.0);
        double vc = _vec_in_var[EQS1::c]->eval(.0);
        double va = 1./2.*(15. - 2.*vb - 3.*vc);
        if (_vec_out_var[EQS1::a]!=0)
            delete _vec_out_var[EQS1::a];
        _vec_out_var[EQS1::a] = new MathLib::FunctionConstant<double,double>(va);
        return 0;
    }
};

class StrongCouplingEQS2 : public IMonolithicProblem
{
public:
    enum InputParameters { a = 0, c = 1 };
    enum OutputParameters { b = 0 };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};
    StrongCouplingEQS2() {
        _vec_in_var.resize(getNumberOfInputVarameters(), 0);
        _vec_out_var.resize(getNumberOfOutputParameters(), 0);
    }
    int solve()
    {
        double va = _vec_in_var[EQS2::a]->eval(.0);
        double vc = _vec_in_var[EQS2::c]->eval(.0);
        double vb = 1./5.*(19.-3*va-2.*vc);
        if (_vec_out_var[EQS2::b]!=0)
            delete _vec_out_var[EQS2::b];
        _vec_out_var[EQS2::b] = new MathLib::FunctionConstant<double,double>(vb);
        return 0;
    }
};

class StrongCouplingEQS3 : public IMonolithicProblem
{
public:
    enum InputParameters { a = 0, b = 1 };
    enum OutputParameters { c = 0 };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};
    StrongCouplingEQS3() {
        _vec_in_var.resize(getNumberOfInputVarameters(), 0);
        _vec_out_var.resize(getNumberOfOutputParameters(), 0);
    }
    int solve()
    {
        double va = _vec_in_var[a]->eval(.0);
        double vb = _vec_in_var[b]->eval(.0);
        double vc = 1./3.*(20.-5*va-3.*vb);
        if (_vec_out_var[c]!=0)
            delete _vec_out_var[c];
        _vec_out_var[c] = new MathLib::FunctionConstant<double,double>(vc);
        return 0;
    }
private:
    double va, vb, vc;
};
#endif


TEST(Coupling, SteadyCouplingCheck1)
{
    WeakCouplingEQS1 eqs1;
    WeakCouplingEQS2 eqs2;
    WeakCouplingEQS3 eqs3;

    eqs1.setParameter(WeakCouplingEQS1::a, new MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(WeakCouplingEQS2::b, new MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(WeakCouplingEQS3::c, new MathLib::FunctionConstant<double,double>(.0));

    {
        // correct
        PartitionedProblem part1(BlockJacobiMethod(1.e-4, 100));
        part1.addParameter("a", eqs1, WeakCouplingEQS1::a);
        part1.addParameter("b", eqs2, WeakCouplingEQS2::b);
        part1.addParameter("c");
        part1.connectInput("b", eqs1, WeakCouplingEQS1::b);
        part1.connectInput("c", eqs1, WeakCouplingEQS1::c);
        part1.connectInput("a", eqs2, WeakCouplingEQS2::a);
        part1.connectInput("c", eqs2, WeakCouplingEQS2::c);

        PartitionedProblem part2(BlockJacobiMethod(1.e-4, 100));
        part2.addParameter("a", part1, part1.getParameterID("a"));
        part2.addParameter("b", part1, part1.getParameterID("b"));
        part2.addParameter("c", eqs3, WeakCouplingEQS3::c);
        part2.connectInput("a", eqs3, WeakCouplingEQS3::a);
        part2.connectInput("b", eqs3, WeakCouplingEQS3::b);
        part2.connectInput("c", part1, part1.getParameterID("c"));

        ASSERT_TRUE(part2.check());
    }

    {
        // no source defined for c
        PartitionedProblem part1(BlockJacobiMethod(1.e-4, 100));
        part1.addParameter("a", eqs1, WeakCouplingEQS1::a);
        part1.addParameter("b", eqs2, WeakCouplingEQS2::b);
        part1.addParameter("c");
        part1.connectInput("a", eqs2, WeakCouplingEQS2::a);
        part1.connectInput("a", eqs3, WeakCouplingEQS3::a);
        part1.connectInput("b", eqs1, WeakCouplingEQS1::b);
        part1.connectInput("b", eqs3, WeakCouplingEQS3::b);
        part1.connectInput("c", eqs1, WeakCouplingEQS1::c);
        part1.connectInput("c", eqs2, WeakCouplingEQS2::c);

        ASSERT_TRUE(part1.check());
    }

    {
        // missing connections eqs2:c
        PartitionedProblem part1(BlockJacobiMethod(1.e-4, 100));
        part1.addParameter("a", eqs1, WeakCouplingEQS1::a);
        part1.addParameter("b", eqs2, WeakCouplingEQS2::b);
        part1.addParameter("c", eqs3, WeakCouplingEQS3::c);
        part1.connectInput("a", eqs2, WeakCouplingEQS2::a);
        part1.connectInput("a", eqs3, WeakCouplingEQS3::a);
        part1.connectInput("b", eqs1, WeakCouplingEQS1::b);
        part1.connectInput("b", eqs3, WeakCouplingEQS3::b);
        part1.connectInput("c", eqs1, WeakCouplingEQS1::c);

        ASSERT_FALSE(part1.check());
    }

}


TEST(Coupling, SteadyCouplingJacobi)
{
    WeakCouplingEQS1 eqs1;
    WeakCouplingEQS2 eqs2;
    WeakCouplingEQS3 eqs3;

    eqs1.setParameter(WeakCouplingEQS1::a, new MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(WeakCouplingEQS2::b, new MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(WeakCouplingEQS3::c, new MathLib::FunctionConstant<double,double>(.0));

    PartitionedProblem part1(BlockJacobiMethod(1.e-4, 100));
    part1.addParameter("a", eqs1, WeakCouplingEQS1::a);
    part1.addParameter("b", eqs2, WeakCouplingEQS2::b);
    part1.addParameter("c");
    part1.connectInput("b", eqs1, WeakCouplingEQS1::b);
    part1.connectInput("c", eqs1, WeakCouplingEQS1::c);
    part1.connectInput("a", eqs2, WeakCouplingEQS2::a);
    part1.connectInput("c", eqs2, WeakCouplingEQS2::c);

    PartitionedProblem part2(BlockJacobiMethod(1.e-4, 100));
    part2.addParameter("a", part1, part1.getParameterID("a"));
    part2.addParameter("b", part1, part1.getParameterID("b"));
    part2.addParameter("c", eqs3, WeakCouplingEQS3::c);
    part2.connectInput("a", eqs3, WeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, WeakCouplingEQS3::b);
    part2.connectInput("c", part1, part1.getParameterID("c"));

    ASSERT_TRUE(part2.check());
    part2.solve();

    const double epsilon = 1.e-3;
    ASSERT_NEAR(1., part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(2., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(3., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);
}

TEST(Coupling, SteadyCouplingSeidel)
{
    WeakCouplingEQS1 eqs1;
    WeakCouplingEQS2 eqs2;
    WeakCouplingEQS3 eqs3;

    eqs1.setParameter(WeakCouplingEQS1::a, new MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(WeakCouplingEQS2::b, new MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(WeakCouplingEQS3::c, new MathLib::FunctionConstant<double,double>(.0));

    PartitionedProblem part1(BlockGaussSeidelMethod(1.e-5, 100));
    part1.addParameter("a", eqs1, WeakCouplingEQS1::a);
    part1.addParameter("b", eqs2, WeakCouplingEQS2::b);
    part1.addParameter("c");
    part1.connectInput("b", eqs1, WeakCouplingEQS1::b);
    part1.connectInput("c", eqs1, WeakCouplingEQS1::c);
    part1.connectInput("a", eqs2, WeakCouplingEQS2::a);
    part1.connectInput("c", eqs2, WeakCouplingEQS2::c);

    PartitionedProblem part2(BlockGaussSeidelMethod(1.e-5, 100));
    part2.addParameter("a", part1, part1.getParameterID("a"));
    part2.addParameter("b", part1, part1.getParameterID("b"));
    part2.addParameter("c", eqs3, WeakCouplingEQS3::c);
    part2.connectInput("a", eqs3, WeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, WeakCouplingEQS3::b);
    part2.connectInput("c", part1, part1.getParameterID("c"));

    ASSERT_TRUE(part2.check());
    part2.solve();

    const double epsilon = 1.e-3;
    ASSERT_NEAR(1., part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(2., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(3., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);
}

// 2a + 2*b + .3*c = 6.9*t 
// 3a + 5b + .2c = 13.6*t
// .5a + .3b + 3c = 10.1*t
// A. a=t, b=2*t, c=3*t
class TransientWeakCouplingEQS1 : public TemplateTransientMonolithicProblem<2,1>
{
    double _dt;
public:
    TransientWeakCouplingEQS1(double dt) : _dt(dt) {};
    enum Parameters { a=2, b = 0, c = 1 };
    int solveTimeStep(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double vb = _vec_parameters[WeakCouplingEQS1::b]->eval(.0);
        double vc = _vec_parameters[WeakCouplingEQS1::c]->eval(.0);
        double va = 1./2.*(6.9*t - 2.*vb - 0.3*vc);
        if (_vec_parameters[WeakCouplingEQS1::a]!=0)
            delete _vec_parameters[WeakCouplingEQS1::a];
        _vec_parameters[WeakCouplingEQS1::a] = new MathLib::FunctionConstant<double,double>(va);
        return 0;
    }
    double suggestNext(const TimeStep &ts) {return ts.getTime()+_dt;};
    bool isAwake(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double mod = t - static_cast<int>(t/_dt)*_dt;
        if (mod==.0) return true;
        return false;
    };
};

class TransientWeakCouplingEQS2 :  public TemplateTransientMonolithicProblem<2,1>
{
    double _dt;
public:
    TransientWeakCouplingEQS2(double dt) : _dt(dt) {};
    enum Parameters { a = 0, b = 2, c = 1 };
    int solveTimeStep(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double va = _vec_parameters[WeakCouplingEQS2::a]->eval(.0);
        double vc = _vec_parameters[WeakCouplingEQS2::c]->eval(.0);
        double vb = 1./5.*(13.6*t-3*va-0.2*vc);
        if (_vec_parameters[WeakCouplingEQS2::b]!=0)
            delete _vec_parameters[WeakCouplingEQS2::b];
        _vec_parameters[WeakCouplingEQS2::b] = new MathLib::FunctionConstant<double,double>(vb);
        return 0;
    }
    double suggestNext(const TimeStep &ts) {return (ts.getTime()+_dt);};
    bool isAwake(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double mod = t - static_cast<int>(t/_dt)*_dt;
        if (mod==.0) return true;
        return false;
    };
};

class TransientWeakCouplingEQS3 : public TemplateTransientMonolithicProblem<2,1>
{
    double _dt;
public:
    TransientWeakCouplingEQS3(double dt) : _dt(dt) {};
    enum Parameters { a = 0, b = 1, c = 2 };
    int solveTimeStep(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double va = _vec_parameters[a]->eval(.0);
        double vb = _vec_parameters[b]->eval(.0);
        double vc = 1./3.*(10.1*t-0.5*va-0.3*vb);
        if (_vec_parameters[c]!=0)
            delete _vec_parameters[c];
        _vec_parameters[c] = new MathLib::FunctionConstant<double,double>(vc);
        return 0;
    }
    double suggestNext(const TimeStep &ts) {return (ts.getTime()+_dt);};
    bool isAwake(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double mod = t - static_cast<int>(t/_dt)*_dt;
        if (mod==.0) return true;
        return false;
    };
};

TEST(Coupling, TransientCouplingParallelStaggered1)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(1.0);
    TransientWeakCouplingEQS3 eqs3(1.0);
    eqs1.setParameter(TransientWeakCouplingEQS1::a, new MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, new MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, new MathLib::FunctionConstant<double,double>(.0));

    AsyncPartitionedProblem apart1(ParallelStaggeredMethod(1e-5, 100));
    apart1.addParameter("a", eqs1, TransientWeakCouplingEQS1::a);
    apart1.addParameter("b", eqs2, TransientWeakCouplingEQS2::b);
    apart1.addParameter("c");
    apart1.connectInput("b", eqs1, TransientWeakCouplingEQS1::b);
    apart1.connectInput("c", eqs1, TransientWeakCouplingEQS1::c);
    apart1.connectInput("a", eqs2, TransientWeakCouplingEQS2::a);
    apart1.connectInput("c", eqs2, TransientWeakCouplingEQS2::c);

    AsyncPartitionedProblem part2(ParallelStaggeredMethod(1.e-5, 100));
    part2.addParameter("a", apart1, apart1.getParameterID("a"));
    part2.addParameter("b", apart1, apart1.getParameterID("b"));
    part2.addParameter("c", eqs3, TransientWeakCouplingEQS3::c);
    part2.connectInput("c", apart1, apart1.getParameterID("c"));
    part2.connectInput("a", eqs3, TransientWeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, TransientWeakCouplingEQS3::b);

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    ASSERT_NEAR(1., part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(2., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(3., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    timestepping.solve(2.0);
    ASSERT_NEAR(2., part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(4., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(6., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);
}

TEST(Coupling, TransientCouplingParallelStaggered2)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);
    eqs1.setParameter(TransientWeakCouplingEQS1::a, new MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, new MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, new MathLib::FunctionConstant<double,double>(.0));

    AsyncPartitionedProblem apart1(ParallelStaggeredMethod(1e-5, 100));
    apart1.addParameter("a", eqs1, TransientWeakCouplingEQS1::a);
    apart1.addParameter("b", eqs2, TransientWeakCouplingEQS2::b);
    apart1.addParameter("c");
    apart1.connectInput("b", eqs1, TransientWeakCouplingEQS1::b);
    apart1.connectInput("c", eqs1, TransientWeakCouplingEQS1::c);
    apart1.connectInput("a", eqs2, TransientWeakCouplingEQS2::a);
    apart1.connectInput("c", eqs2, TransientWeakCouplingEQS2::c);

    AsyncPartitionedProblem part2(ParallelStaggeredMethod(1.e-5, 100));
    part2.addParameter("a", apart1, apart1.getParameterID("a"));
    part2.addParameter("b", apart1, apart1.getParameterID("b"));
    part2.addParameter("c", eqs3, TransientWeakCouplingEQS3::c);
    part2.connectInput("c", apart1, apart1.getParameterID("c"));
    part2.connectInput("a", eqs3, TransientWeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, TransientWeakCouplingEQS3::b);

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    ASSERT_NEAR(6.9/2., part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    timestepping.solve(2.0);
    ASSERT_NEAR(3.65, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(3.25, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    timestepping.solve(3.0);
    ASSERT_NEAR(7.1, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(3.25, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    timestepping.solve(4.0);
    ASSERT_NEAR(4., part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(8., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(12., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    timestepping.solve(5.0);
    ASSERT_NEAR(7.45, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(8., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(12., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);
}

TEST(Coupling, TransientCouplingParallelStaggered3)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);
    eqs1.setParameter(TransientWeakCouplingEQS1::a, &MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, &MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, &MathLib::FunctionConstant<double,double>(.0));

    AsyncPartitionedProblem apart1(ParallelStaggeredMethod(1e-5, 1));
    apart1.addParameter("a", eqs1, TransientWeakCouplingEQS1::a);
    apart1.addParameter("b", eqs2, TransientWeakCouplingEQS2::b);
    apart1.addParameter("c");
    apart1.connectInput("b", eqs1, TransientWeakCouplingEQS1::b);
    apart1.connectInput("c", eqs1, TransientWeakCouplingEQS1::c);
    apart1.connectInput("a", eqs2, TransientWeakCouplingEQS2::a);
    apart1.connectInput("c", eqs2, TransientWeakCouplingEQS2::c);

    AsyncPartitionedProblem part2(ParallelStaggeredMethod(1.e-5, 1));
    part2.addParameter("a", apart1, apart1.getParameterID("a"));
    part2.addParameter("b", apart1, apart1.getParameterID("b"));
    part2.addParameter("c", eqs3, TransientWeakCouplingEQS3::c);
    part2.connectInput("c", apart1, apart1.getParameterID("c"));
    part2.connectInput("a", eqs3, TransientWeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, TransientWeakCouplingEQS3::b);

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    ASSERT_NEAR(3.45, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    eqs1.setParameter(TransientWeakCouplingEQS1::a, &MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, &MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, &MathLib::FunctionConstant<double,double>(.0));
    timestepping.solve(2.0);
    ASSERT_NEAR(6.9, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(3.37, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    eqs1.setParameter(TransientWeakCouplingEQS1::a, &MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, &MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, &MathLib::FunctionConstant<double,double>(.0));
    timestepping.solve(3.0);
    ASSERT_NEAR(6.98, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(3.37, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    eqs1.setParameter(TransientWeakCouplingEQS1::a, &MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, &MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, &MathLib::FunctionConstant<double,double>(.0));
    timestepping.solve(4.0);
    ASSERT_NEAR(10.43, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(6.692, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(11.96633, part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    eqs1.setParameter(TransientWeakCouplingEQS1::a, &MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, &MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, &MathLib::FunctionConstant<double,double>(.0));
    timestepping.solve(5.0);
    ASSERT_NEAR(8.76305, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(6.692, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(11.96633, part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);
}

TEST(Coupling, TransientCouplingSerialStaggered1)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);
    eqs1.setParameter(TransientWeakCouplingEQS1::a, new MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, new MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, new MathLib::FunctionConstant<double,double>(.0));

    AsyncPartitionedProblem apart1(SerialStaggeredMethod(1e-5, 100));
    apart1.addParameter("a", eqs1, TransientWeakCouplingEQS1::a);
    apart1.addParameter("b", eqs2, TransientWeakCouplingEQS2::b);
    apart1.addParameter("c");
    apart1.connectInput("b", eqs1, TransientWeakCouplingEQS1::b);
    apart1.connectInput("c", eqs1, TransientWeakCouplingEQS1::c);
    apart1.connectInput("a", eqs2, TransientWeakCouplingEQS2::a);
    apart1.connectInput("c", eqs2, TransientWeakCouplingEQS2::c);

    AsyncPartitionedProblem part2(SerialStaggeredMethod(1.e-5, 100));
    part2.addParameter("a", apart1, apart1.getParameterID("a"));
    part2.addParameter("b", apart1, apart1.getParameterID("b"));
    part2.addParameter("c", eqs3, TransientWeakCouplingEQS3::c);
    part2.connectInput("c", apart1, apart1.getParameterID("c"));
    part2.connectInput("a", eqs3, TransientWeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, TransientWeakCouplingEQS3::b);

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    ASSERT_NEAR(6.9/2., part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    timestepping.solve(2.0);
    ASSERT_NEAR(3.65, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(3.25, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    timestepping.solve(3.0);
    ASSERT_NEAR(7.1, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(3.25, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    timestepping.solve(4.0);
    ASSERT_NEAR(4., part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(8., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(12., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    timestepping.solve(5.0);
    ASSERT_NEAR(7.45, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(8., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(12., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);
}

TEST(Coupling, TransientCouplingSerialStaggered2)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);
    eqs1.setParameter(TransientWeakCouplingEQS1::a, &MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, &MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, &MathLib::FunctionConstant<double,double>(.0));

    AsyncPartitionedProblem apart1(SerialStaggeredMethod(1e-5, 1));
    apart1.addParameter("a", eqs1, TransientWeakCouplingEQS1::a);
    apart1.addParameter("b", eqs2, TransientWeakCouplingEQS2::b);
    apart1.addParameter("c");
    apart1.connectInput("b", eqs1, TransientWeakCouplingEQS1::b);
    apart1.connectInput("c", eqs1, TransientWeakCouplingEQS1::c);
    apart1.connectInput("a", eqs2, TransientWeakCouplingEQS2::a);
    apart1.connectInput("c", eqs2, TransientWeakCouplingEQS2::c);

    AsyncPartitionedProblem part2(SerialStaggeredMethod(1.e-5, 1));
    part2.addParameter("a", apart1, apart1.getParameterID("a"));
    part2.addParameter("b", apart1, apart1.getParameterID("b"));
    part2.addParameter("c", eqs3, TransientWeakCouplingEQS3::c);
    part2.connectInput("c", apart1, apart1.getParameterID("c"));
    part2.connectInput("a", eqs3, TransientWeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, TransientWeakCouplingEQS3::b);

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    ASSERT_NEAR(3.45, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    eqs1.setParameter(TransientWeakCouplingEQS1::a, &MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, &MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, &MathLib::FunctionConstant<double,double>(.0));
    timestepping.solve(2.0);
    ASSERT_NEAR(6.9, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(1.3, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    eqs1.setParameter(TransientWeakCouplingEQS1::a, &MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, &MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, &MathLib::FunctionConstant<double,double>(.0));
    timestepping.solve(3.0);
    ASSERT_NEAR(9.05, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(1.3, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(0., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    eqs1.setParameter(TransientWeakCouplingEQS1::a, &MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, &MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, &MathLib::FunctionConstant<double,double>(.0));
    timestepping.solve(4.0);
    ASSERT_NEAR(12.5, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(3.38, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(11.04533, part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);

    eqs1.setParameter(TransientWeakCouplingEQS1::a, &MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(TransientWeakCouplingEQS2::b, &MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(TransientWeakCouplingEQS3::c, &MathLib::FunctionConstant<double,double>(.0));
    timestepping.solve(5.0);
    ASSERT_NEAR(12.2132, part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(3.38, part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(11.04533, part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);
}
