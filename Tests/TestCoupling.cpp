
#include <gtest/gtest.h>
#include <vector>

#include "Base/Options.h"
#include "MathLib/Function/Function.h"
#include "MathLib/Coupling/MonolithicProblem.h"
#include "MathLib/Coupling/PartitionedProblem.h"
#include "MathLib/Coupling/Algorithm/PartitionedAlgorithm.h"
#include "MathLib/Coupling/Algorithm/TransientPartitionedMethod.h"
#include "NumLib/TransientCoupling/AsyncPartitionedSystem.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"

using namespace MathLib;
using namespace NumLib;

typedef MathLib::FunctionConstant<double,double> MyFunction;

class MyConvergenceCheck
{
public:
	bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff)
	{
	    for (size_t i=0; i<vars_prev.size(); i++) {
	        double v_prev = .0;
	        vars_prev.get<MyFunction>(i)->eval(v_prev);
	        double v_cur = .0;
	        vars_current.get<MyFunction>(i)->eval(v_cur);
	        v_diff = std::abs(v_cur - v_prev);
	        if (v_diff>eps) {
	            return false;
	        }
	    }
	    return true;
	}
};

// 2a + 2b + .3c = 6.9
// 3a + 5b + .2c = 13.6
// .5a + .3b + 3c = 10.1
// A. a=1, b=2, c=3
class WeakCouplingEQS1 : public TemplateSteadyMonolithicSystem
{
	enum In {b=0, c=1};
	enum Out {a=0};
public:
    WeakCouplingEQS1() 
    {
    	resizeInputParameter(2);
    	resizeOutputParameter(1);

        setOutput(a, new MathLib::FunctionConstant<double,double>(.0));
    }

    int solve()
    {
    	double vb, vc;
        getInput<MyFunction>(b)->eval(vb);
        getInput<MyFunction>(c)->eval(vc);
        double va = 1./2.*(6.9 - 2.*vb - 0.3*vc);
        setOutput(a, new MathLib::FunctionConstant<double,double>(va));
        return 0;
    }
};

class WeakCouplingEQS2 :  public TemplateSteadyMonolithicSystem
{
	enum In {a=0, c=1};
	enum Out {b=0};
public:

    WeakCouplingEQS2() 
    {
    	resizeInputParameter(2);
    	resizeOutputParameter(1);
        setOutput(b, new MathLib::FunctionConstant<double,double>(.0));
    }

    int solve()
    {
    	double va, vc;
        getInput<MyFunction>(a)->eval(va);
        getInput<MyFunction>(c)->eval(vc);
        double vb = 1./5.*(13.6-3*va-0.2*vc);
        setOutput(b, new MathLib::FunctionConstant<double,double>(vb));
        return 0;
    }
};

class WeakCouplingEQS3 : public TemplateSteadyMonolithicSystem
{
	enum In {a=0, b=1};
	enum Out {c=0};
public:

    WeakCouplingEQS3() 
    {
    	resizeInputParameter(2);
    	resizeOutputParameter(1);
        setOutput(c, new MathLib::FunctionConstant<double,double>(.0));
    }

    int solve()
    {
    	double va, vb;
        getInput<MyFunction>(a)->eval(va);
        getInput<MyFunction>(b)->eval(vb);
        double vc = 1./3.*(10.1-0.5*va-0.3*vb);
        setOutput(c, new MathLib::FunctionConstant<double,double>(vc));
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

Base::Options* defineOption()
{
	Base::Options* options = new Base::Options();
	Base::Options* coupling = options->addSubGroup("coupling");
	Base::Options* P2 = coupling->addSubGroup("P");
	{
	//P2->addOption("name", "P2");
	P2->addOption("algorithm", "Jacobi");
	P2->addOptionAsNum("max_itr", 100);
	P2->addOptionAsNum("epsilon", 1.e-4);
	std::vector<std::string> out_var;
	out_var.push_back("a");
	out_var.push_back("b");
	out_var.push_back("c");
	P2->addOptionAsArray("out", out_var);
	}
	Base::Options* P2_sub = P2->addSubGroup("problems");
	Base::Options* P1 = P2_sub->addSubGroup("P");
	{
	//P1->addOption("name", "P1");
	P1->addOption("algorithm", "Gauss");
	P1->addOptionAsNum("max_itr", 100);
	P1->addOptionAsNum("epsilon", 1.e-4);
	std::vector<std::string> out_var;
	out_var.push_back("a");
	out_var.push_back("b");
	std::vector<std::string> in_var;
	in_var.push_back("c");
	P1->addOptionAsArray("out", out_var);
	P1->addOptionAsArray("in", in_var);
	}
	Base::Options* P1_sub = P1->addSubGroup("problems");
	Base::Options* M1 = P1_sub->addSubGroup("M1");
	{
	M1->addOption("name", "EQS1");
	std::vector<std::string> out_var;
	out_var.push_back("a");
	std::vector<std::string> in_var;
	in_var.push_back("b");
	in_var.push_back("c");
	M1->addOptionAsArray("out", out_var);
	M1->addOptionAsArray("in", in_var);
	}
	Base::Options* M2 = P1_sub->addSubGroup("M2");
	{
	M2->addOption("name", "EQS2");
	std::vector<std::string> out_var;
	out_var.push_back("b");
	std::vector<std::string> in_var;
	in_var.push_back("a");
	in_var.push_back("c");
	M2->addOptionAsArray("out", out_var);
	M2->addOptionAsArray("in", in_var);
	}
	Base::Options* M3 = P2_sub->addSubGroup("M3");
	{
	M3->addOption("name", "EQS3");
	std::vector<std::string> out_var;
	out_var.push_back("c");
	std::vector<std::string> in_var;
	in_var.push_back("a");
	in_var.push_back("b");
	M3->addOptionAsArray("out", out_var);
	M3->addOptionAsArray("in", in_var);
	}

	return options;
};

class EQSFactory
{
public:
	static TemplateSteadyMonolithicSystem* create(const std::string &eqs_name)
	{
		if (eqs_name.compare("EQS1")==0) {
			return new WeakCouplingEQS1();
		} else if (eqs_name.compare("EQS2")==0) {
				return new WeakCouplingEQS2();
		} else if (eqs_name.compare("EQS3")==0) {
			return new WeakCouplingEQS3();
		}
		return 0;
	};
};

class CouplingAlgorithmFactory
{
public:
	template<class T>
	static IPartitionedAlgorithm* create(const std::string &name, size_t max_itr, double epsilon)
	{
		if (name.compare("Jacobi")==0) {
			return new BlockJacobiMethod<T>(epsilon, max_itr);
		} else if (name.compare("Gauss")==0) {
			return new BlockGaussSeidelMethod<T>(epsilon, max_itr);
		}
		return 0;
	};
};

TemplateSteadyMonolithicSystem* buildMonolithicSystem(const Base::Options *option)
{
	TemplateSteadyMonolithicSystem* eqs = EQSFactory::create(option->getOption("name"));
	const std::vector<std::string>* in_names = option->getOptionAsArray<std::string>("in");
	const std::vector<std::string>* out_names = option->getOptionAsArray<std::string>("out");
	if (in_names!=0) {
		for (size_t i=0; i<in_names->size(); i++) {
			eqs->setInputParameterName(i, (*in_names)[i]);
		}
	}
	if (out_names!=0) {
		for (size_t i=0; i<out_names->size(); i++) {
			eqs->setOutputParameterName(i, (*out_names)[i]);
		}
	}
	return eqs;
}

PartitionedProblem* buildPartitionedSystem(const Base::Options *option)
{
	PartitionedProblem* part = new PartitionedProblem();
	//para
	const std::vector<std::string>* in_names = option->getOptionAsArray<std::string>("in");
	const std::vector<std::string>* out_names = option->getOptionAsArray<std::string>("out");
	if (in_names!=0) {
        part->resizeInputParameter(in_names->size());
		for (size_t i=0; i<in_names->size(); i++) {
			part->setInputParameterName(i, (*in_names)[i]);
		}
	}
	if (out_names!=0) {
        part->resizeOutputParameter(out_names->size());
		for (size_t i=0; i<out_names->size(); i++) {
			part->setOutputParameterName(i, (*out_names)[i]);
		}
	}
	//alg
	size_t max_itr = option->getOption<size_t>("max_itr");
	double epsilon = option->getOption<double>("epsilon");
	IPartitionedAlgorithm* alg = CouplingAlgorithmFactory::create<MyConvergenceCheck>(option->getOption("algorithm"), max_itr, epsilon);
	if (alg!=0) part->setAlgorithm(*alg);
	//problems
	const Base::Options* op_problems = option->getSubGroup("problems");
	for (Base::Options::const_iterator itr=op_problems->begin(); itr!=op_problems->end(); ++itr) {
		std::string str = itr->first;
		Base::Options* op_sub = static_cast<Base::Options*>(itr->second);
		ICoupledSystem* sys = 0;
		if (str.find("M")==0) {
			sys = buildMonolithicSystem(op_sub);
		} else if (str.compare("P")==0) {
			sys = buildPartitionedSystem(op_sub);
		}
        if (sys!=0) part->addProblem(*sys);
	}
    part->connectParameters();
	return part;
}


ICoupledSystem* buildCoupling(const Base::Options *option)
{
	const Base::Options* op_cpl = option->getSubGroup("coupling");
    if (op_cpl->hasSubGroup("M")) {
        const Base::Options* op_sub = op_cpl->getSubGroup("M");
        TemplateSteadyMonolithicSystem *sys = buildMonolithicSystem(op_sub);
        return sys;
    } else if (op_cpl->hasSubGroup("P")) {
        const Base::Options* op_sub = op_cpl->getSubGroup("P");
        PartitionedProblem *sys = buildPartitionedSystem(op_sub);
        return sys;
    }
    return 0;
}

TEST(Coupling, Option)
{
	Base::Options *option = defineOption();
	ICoupledSystem *coupled_sys = buildCoupling(option);
	ASSERT_TRUE(coupled_sys->check());
    coupled_sys->solve();

    const double epsilon = 1.e-3;
    double v1, v2, v3;
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("a")))->eval(v1);
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("b")))->eval(v2);
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("c")))->eval(v3);
    ASSERT_NEAR(1., v1, epsilon);
    ASSERT_NEAR(2., v2, epsilon);
    ASSERT_NEAR(3., v3, epsilon);

}

TEST(Coupling, SteadyCouplingCheck1)
{
    WeakCouplingEQS1 eqs1;
    eqs1.setOutputParameterName(0,"a");
    eqs1.setInputParameterName(0,"b");
    eqs1.setInputParameterName(1,"c");
    WeakCouplingEQS2 eqs2;
    eqs2.setOutputParameterName(0,"b");
    eqs2.setInputParameterName(0,"a");
    eqs2.setInputParameterName(1,"c");
    WeakCouplingEQS3 eqs3;
    eqs3.setOutputParameterName(0,"c");
    eqs3.setInputParameterName(0,"a");
    eqs3.setInputParameterName(1,"b");

    ASSERT_TRUE(eqs1.isValid());
    ASSERT_TRUE(eqs2.isValid());
    ASSERT_TRUE(eqs3.isValid());

    {
        // correct
    	BlockJacobiMethod<MyConvergenceCheck> method(1.e-4, 100);
        PartitionedProblem part1;
        part1.setAlgorithm(method);
        part1.resizeInputParameter(1);
        part1.resizeOutputParameter(2);
        part1.setOutputParameterName(0, "a");
        part1.setOutputParameterName(1, "b");
        part1.setInputParameterName(0, "c");
        part1.addProblem(eqs1);
        part1.addProblem(eqs2);
        part1.connectParameters();

        PartitionedProblem part2;
        part2.setAlgorithm(method);
        part2.resizeInputParameter(0);
        part2.resizeOutputParameter(3);
        part2.setOutputParameterName(0, "a");
        part2.setOutputParameterName(1, "b");
        part2.setOutputParameterName(2, "c");
        part2.addProblem(part1);
        part2.addProblem(eqs3);
        part2.connectParameters();

        ASSERT_TRUE(part2.check());
    }

    {
        // no source defined for c
    	BlockJacobiMethod<MyConvergenceCheck> method(1.e-4, 100);
        PartitionedProblem part1;
        part1.setAlgorithm(method);
        part1.resizeInputParameter(0);
        part1.resizeOutputParameter(2);
        part1.setOutputParameterName(0, "a");
        part1.setOutputParameterName(1, "b");
        //part1.setInputParameterName(0, "c");
        part1.addProblem(eqs1);
        part1.addProblem(eqs2);
        part1.connectParameters();

        ASSERT_FALSE(part1.check());
    }

    //{
    //    // no source defined for c
    //    BlockJacobiMethod<MyConvergenceCheck> method(1.e-4, 100);
    //    PartitionedProblem part1;
    //    part1.setAlgorithm(method);
    //    part1.resizeInputParameter(1);
    //    part1.resizeOutputParameter(2);
    //    part1.setOutputParameterName(0, "a");
    //    part1.setOutputParameterName(1, "b");
    //    part1.setInputParameterName(0, "c");
    //    part1.addProblem(eqs1);
    //    part1.addProblem(eqs2);
    //    part1.connectParameters();

    //    PartitionedProblem part2;
    //    part2.setAlgorithm(method);
    //    part2.resizeInputParameter(0);
    //    part2.resizeOutputParameter(3);
    //    part2.setOutputParameterName(0, "a");
    //    part2.setOutputParameterName(1, "b");
    //    part2.setOutputParameterName(2, "c");
    //    part2.addProblem(part1);
    //    part2.connectParameters();

    //    ASSERT_FALSE(part2.check());
    //}
}

#if 0

TEST(Coupling, SteadyCouplingJacobi)
{
    WeakCouplingEQS1 eqs1;
    WeakCouplingEQS2 eqs2;
    WeakCouplingEQS3 eqs3;

    BlockJacobiMethod<MyConvergenceCheck> method(1.e-4, 100);
    PartitionedProblem part1(method);
    //PartitionedProblem part1(BlockJacobiMethod(1.e-4, 100));
    part1.addOutputParameter("a", eqs1, eqs1.a);
    part1.addOutputParameter("b", eqs2, eqs2.b);
    part1.addInputParameter("c");
    part1.connectInput("b", eqs1, eqs1.b);
    part1.connectInput("c", eqs1, eqs1.c);
    part1.connectInput("a", eqs2, eqs2.a);
    part1.connectInput("c", eqs2, eqs2.c);

    PartitionedProblem part2(method);
    //PartitionedProblem part2(BlockJacobiMethod(1.e-4, 100));
    part2.addOutputParameter("a", part1, part1.getOutputParameterID("a"));
    part2.addOutputParameter("b", part1, part1.getOutputParameterID("b"));
    part2.addOutputParameter("c", eqs3, eqs3.c);
    part2.connectInput("a", eqs3, eqs3.a);
    part2.connectInput("b", eqs3, eqs3.b);
    part2.connectInput("c", part1, part1.getOutputParameterID("c"));

    ASSERT_TRUE(part2.check());
    part2.solve();

    const double epsilon = 1.e-3;
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(v3);
    ASSERT_NEAR(1., v1, epsilon);
    ASSERT_NEAR(2., v2, epsilon);
    ASSERT_NEAR(3., v3, epsilon);
}

TEST(Coupling, SteadyCouplingSeidel)
{
    WeakCouplingEQS1 eqs1;
    WeakCouplingEQS2 eqs2;
    WeakCouplingEQS3 eqs3;

    BlockGaussSeidelMethod<MyConvergenceCheck> method(1.e-5, 100);
    PartitionedProblem part1(method);
    //PartitionedProblem part1(BlockGaussSeidelMethod(1.e-5, 100));
    part1.addOutputParameter("a", eqs1, eqs1.a);
    part1.addOutputParameter("b", eqs2, eqs2.b);
    part1.addInputParameter("c");
    part1.connectInput("b", eqs1, eqs1.b);
    part1.connectInput("c", eqs1, eqs1.c);
    part1.connectInput("a", eqs2, eqs2.a);
    part1.connectInput("c", eqs2, eqs2.c);

    PartitionedProblem part2(method);
    part2.addOutputParameter("a", part1, part1.getOutputParameterID("a"));
    part2.addOutputParameter("b", part1, part1.getOutputParameterID("b"));
    part2.addOutputParameter("c", eqs3, eqs3.c);
    part2.connectInput("a", eqs3, eqs3.a);
    part2.connectInput("b", eqs3, eqs3.b);
    part2.connectInput("c", part1, part1.getOutputParameterID("c"));

    ASSERT_TRUE(part2.check());
    part2.solve();

    const double epsilon = 1.e-3;
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(v3);
    ASSERT_NEAR(1., v1, epsilon);
    ASSERT_NEAR(2., v2, epsilon);
    ASSERT_NEAR(3., v3, epsilon);
}



// 2a + 2*b + .3*c = 6.9*t 
// 3a + 5b + .2c = 13.6*t
// .5a + .3b + 3c = 10.1*t
// A. a=t, b=2*t, c=3*t
class TransientWeakCouplingEQS1 : public TemplateTransientMonolithicSystem<2,1>
{
    double _dt;
public:
    enum Parameters { a=2, b = 0, c = 1 };
    TransientWeakCouplingEQS1(double dt) : _dt(dt) 
    {
        reset();
    };

    void reset()
    {
        setOutput(TransientWeakCouplingEQS1::a, new MathLib::FunctionConstant<double,double>(.0));
    }

    int solveTimeStep(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double vb, vc;
        getInput<MyFunction>(b)->eval(vb);
        getInput<MyFunction>(c)->eval(vc);
        double va = 1./2.*(6.9*t - 2.*vb - 0.3*vc);
        setOutput(a, new MathLib::FunctionConstant<double,double>(va));
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
    void accept(const TimeStep &time) {};
};

class TransientWeakCouplingEQS2 :  public TemplateTransientMonolithicSystem<2,1>
{
    double _dt;
public:
    enum Parameters { a = 0, b = 2, c = 1 };

    TransientWeakCouplingEQS2(double dt) : _dt(dt) 
    {
        reset();
    };

    void reset()
    {
        setOutput(b, new MathLib::FunctionConstant<double,double>(.0));
    }

    int solveTimeStep(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double va, vc;
        getInput<MyFunction>(a)->eval(va);
        getInput<MyFunction>(c)->eval(vc);
        double vb = 1./5.*(13.6*t-3*va-0.2*vc);
        setOutput(b, new MathLib::FunctionConstant<double,double>(vb));
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
    void accept(const TimeStep &time) {};
};

class TransientWeakCouplingEQS3 : public TemplateTransientMonolithicSystem<2,1>
{
    double _dt;
public:
    enum Parameters { a = 0, b = 1, c = 2 };

    TransientWeakCouplingEQS3(double dt) : _dt(dt) 
    {
        reset();
    };

    void reset()
    {
        setOutput(TransientWeakCouplingEQS3::c, new MathLib::FunctionConstant<double,double>(.0));
    }

    int solveTimeStep(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double va, vb;
        getInput<MyFunction>(a)->eval(va);
        getInput<MyFunction>(b)->eval(vb);
        double vc = 1./3.*(10.1*t-0.5*va-0.3*vb);
        setOutput(c, new MathLib::FunctionConstant<double,double>(vc));
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
    void accept(const TimeStep &time) {};
};

TEST(Coupling, TransientCouplingParallelStaggered1)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(1.0);
    TransientWeakCouplingEQS3 eqs3(1.0);

    ParallelStaggeredMethod<MyConvergenceCheck> method(1e-5, 100);
    AsyncPartitionedSystem apart1(method);
    //AsyncPartitionedSystem apart1(ParallelStaggeredMethod(1e-5, 100));
    apart1.addOutputParameter("a", eqs1, TransientWeakCouplingEQS1::a);
    apart1.addOutputParameter("b", eqs2, TransientWeakCouplingEQS2::b);
    apart1.addInputParameter("c");
    apart1.connectInput("b", eqs1, TransientWeakCouplingEQS1::b);
    apart1.connectInput("c", eqs1, TransientWeakCouplingEQS1::c);
    apart1.connectInput("a", eqs2, TransientWeakCouplingEQS2::a);
    apart1.connectInput("c", eqs2, TransientWeakCouplingEQS2::c);

    AsyncPartitionedSystem part2(method);
    part2.addOutputParameter("a", apart1, apart1.getOutputParameterID("a"));
    part2.addOutputParameter("b", apart1, apart1.getOutputParameterID("b"));
    part2.addOutputParameter("c", eqs3, TransientWeakCouplingEQS3::c);
    part2.connectInput("c", apart1, apart1.getOutputParameterID("c"));
    part2.connectInput("a", eqs3, TransientWeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, TransientWeakCouplingEQS3::b);

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(1., v1, epsilon);
    ASSERT_NEAR(2., v2, epsilon);
    ASSERT_NEAR(3., v3, epsilon);

    timestepping.solve(2.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(2., v1, epsilon);
    ASSERT_NEAR(4., v2, epsilon);
    ASSERT_NEAR(6., v3, epsilon);
}

TEST(Coupling, TransientCouplingParallelStaggered2)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);

    ParallelStaggeredMethod<MyConvergenceCheck> method(1e-5, 100);
    AsyncPartitionedSystem apart1(method);
    //AsyncPartitionedSystem apart1(ParallelStaggeredMethod(1e-5, 100));
    apart1.addOutputParameter("a", eqs1, TransientWeakCouplingEQS1::a);
    apart1.addOutputParameter("b", eqs2, TransientWeakCouplingEQS2::b);
    apart1.addInputParameter("c");
    apart1.connectInput("b", eqs1, TransientWeakCouplingEQS1::b);
    apart1.connectInput("c", eqs1, TransientWeakCouplingEQS1::c);
    apart1.connectInput("a", eqs2, TransientWeakCouplingEQS2::a);
    apart1.connectInput("c", eqs2, TransientWeakCouplingEQS2::c);

    AsyncPartitionedSystem part2(method);
    part2.addOutputParameter("a", apart1, apart1.getOutputParameterID("a"));
    part2.addOutputParameter("b", apart1, apart1.getOutputParameterID("b"));
    part2.addOutputParameter("c", eqs3, TransientWeakCouplingEQS3::c);
    part2.connectInput("c", apart1, apart1.getOutputParameterID("c"));
    part2.connectInput("a", eqs3, TransientWeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, TransientWeakCouplingEQS3::b);

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.9/2., v1, epsilon);
    ASSERT_NEAR(0., v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(2.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(3.65, v1, epsilon);
    ASSERT_NEAR(3.25, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(3.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(7.1, v1, epsilon);
    ASSERT_NEAR(3.25, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(4.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(4., v1, epsilon);
    ASSERT_NEAR(8., v2, epsilon);
    ASSERT_NEAR(12., v3, epsilon);

    timestepping.solve(5.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(7.45, v1, epsilon);
    ASSERT_NEAR(8., v2, epsilon);
    ASSERT_NEAR(12., v3, epsilon);
}

TEST(Coupling, TransientCouplingParallelStaggered3)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);
    MathLib::FunctionConstant<double,double> f_const_zero(.0);

    ParallelStaggeredMethod<MyConvergenceCheck> method(1e-5, 1);
    AsyncPartitionedSystem apart1(method);
    //AsyncPartitionedSystem apart1(ParallelStaggeredMethod(1e-5, 1));
    apart1.addOutputParameter("a", eqs1, TransientWeakCouplingEQS1::a);
    apart1.addOutputParameter("b", eqs2, TransientWeakCouplingEQS2::b);
    apart1.addInputParameter("c");
    apart1.connectInput("b", eqs1, TransientWeakCouplingEQS1::b);
    apart1.connectInput("c", eqs1, TransientWeakCouplingEQS1::c);
    apart1.connectInput("a", eqs2, TransientWeakCouplingEQS2::a);
    apart1.connectInput("c", eqs2, TransientWeakCouplingEQS2::c);

    AsyncPartitionedSystem part2(method);
    part2.addOutputParameter("a", apart1, apart1.getOutputParameterID("a"));
    part2.addOutputParameter("b", apart1, apart1.getOutputParameterID("b"));
    part2.addOutputParameter("c", eqs3, TransientWeakCouplingEQS3::c);
    part2.connectInput("c", apart1, apart1.getOutputParameterID("c"));
    part2.connectInput("a", eqs3, TransientWeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, TransientWeakCouplingEQS3::b);

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(3.45, v1, epsilon);
    ASSERT_NEAR(0., v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(2.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.9, v1, epsilon);
    ASSERT_NEAR(3.37, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(3.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.98, v1, epsilon);
    ASSERT_NEAR(3.37, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(4.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(10.43, v1, epsilon);
    ASSERT_NEAR(6.692, v2, epsilon);
    ASSERT_NEAR(11.96633, v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(5.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(8.76305, v1, epsilon);
    ASSERT_NEAR(6.692, v2, epsilon);
    ASSERT_NEAR(11.96633, v3, epsilon);
}

TEST(Coupling, TransientCouplingSerialStaggered1)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);

    SerialStaggeredMethod<MyConvergenceCheck> method(1e-5, 100);
    AsyncPartitionedSystem apart1(method);
    //AsyncPartitionedSystem apart1(SerialStaggeredMethod(1e-5, 100));
    apart1.addOutputParameter("a", eqs1, TransientWeakCouplingEQS1::a);
    apart1.addOutputParameter("b", eqs2, TransientWeakCouplingEQS2::b);
    apart1.addInputParameter("c");
    apart1.connectInput("b", eqs1, TransientWeakCouplingEQS1::b);
    apart1.connectInput("c", eqs1, TransientWeakCouplingEQS1::c);
    apart1.connectInput("a", eqs2, TransientWeakCouplingEQS2::a);
    apart1.connectInput("c", eqs2, TransientWeakCouplingEQS2::c);

    AsyncPartitionedSystem part2(method);
    part2.addOutputParameter("a", apart1, apart1.getOutputParameterID("a"));
    part2.addOutputParameter("b", apart1, apart1.getOutputParameterID("b"));
    part2.addOutputParameter("c", eqs3, TransientWeakCouplingEQS3::c);
    part2.connectInput("c", apart1, apart1.getOutputParameterID("c"));
    part2.connectInput("a", eqs3, TransientWeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, TransientWeakCouplingEQS3::b);

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.9/2., v1, epsilon);
    ASSERT_NEAR(0., v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(2.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(3.65, v1, epsilon);
    ASSERT_NEAR(3.25, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(3.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(7.1, v1, epsilon);
    ASSERT_NEAR(3.25, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(4.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(4, v1, epsilon);
    ASSERT_NEAR(8, v2, epsilon);
    ASSERT_NEAR(12., v3, epsilon);

    timestepping.solve(5.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(7.45, v1, epsilon);
    ASSERT_NEAR(8, v2, epsilon);
    ASSERT_NEAR(12., v3, epsilon);
}

TEST(Coupling, TransientCouplingSerialStaggered2)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);
    MathLib::FunctionConstant<double,double> f_const_zero(.0);

    SerialStaggeredMethod<MyConvergenceCheck> method(1e-5, 1);
    AsyncPartitionedSystem apart1(method);
    //AsyncPartitionedSystem apart1(SerialStaggeredMethod(1e-5, 1));
    apart1.addOutputParameter("a", eqs1, TransientWeakCouplingEQS1::a);
    apart1.addOutputParameter("b", eqs2, TransientWeakCouplingEQS2::b);
    apart1.addInputParameter("c");
    apart1.connectInput("b", eqs1, TransientWeakCouplingEQS1::b);
    apart1.connectInput("c", eqs1, TransientWeakCouplingEQS1::c);
    apart1.connectInput("a", eqs2, TransientWeakCouplingEQS2::a);
    apart1.connectInput("c", eqs2, TransientWeakCouplingEQS2::c);

    AsyncPartitionedSystem part2(method);
    part2.addOutputParameter("a", apart1, apart1.getOutputParameterID("a"));
    part2.addOutputParameter("b", apart1, apart1.getOutputParameterID("b"));
    part2.addOutputParameter("c", eqs3, TransientWeakCouplingEQS3::c);
    part2.connectInput("c", apart1, apart1.getOutputParameterID("c"));
    part2.connectInput("a", eqs3, TransientWeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, TransientWeakCouplingEQS3::b);

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.9/2., v1, epsilon);
    ASSERT_NEAR(0., v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(2.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.9, v1, epsilon);
    ASSERT_NEAR(1.3, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(3.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(9.05, v1, epsilon);
    ASSERT_NEAR(1.3, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(4.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(12.5, v1, epsilon);
    ASSERT_NEAR(3.38, v2, epsilon);
    ASSERT_NEAR(11.04533, v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(5.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0., v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0., v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0., v3);
    ASSERT_NEAR(12.2132, v1, epsilon);
    ASSERT_NEAR(3.38, v2, epsilon);
    ASSERT_NEAR(11.04533, v3, epsilon);
}

#endif
