
#include <gtest/gtest.h>

#include "MathLib/Coupling/PartitionedProblem.h"
#include "MathLib/Coupling/Algorithm/BlockJacobiMethod.h"
#include "MathLib/SystemOfEquations/SystemOfEquations.h"
#include "MathLib/SystemOfEquations/CoupledProblemConstructor.h"

#include "TestUtil.h"

using namespace MathLib;


void defineProblem1(SystemOfEquations &sysEqs)
{
	// 2a + 2b + .3c = 6.9
	// 3a + 5b + .2c = 13.6
	// .5a + .3b + 3c = 10.1
	// A. a=1, b=2, c=3
	std::vector<Matrix<double>* > vecCp(3);
	for (size_t i=0; i<3; i++) vecCp[i] = new Matrix<double>(1,1);
	std::valarray<double>* Fp = new std::valarray<double>(1);
	(*vecCp[0])(0,0) = 2;
	(*vecCp[1])(0,0) = 2;
	(*vecCp[2])(0,0) = .3;
	Fp[0] = 6.9;
	std::vector<Matrix<double>* > vecCT(3);
	for (size_t i=0; i<3; i++) vecCT[i] = new Matrix<double>(1,1);
	std::valarray<double>* FT = new std::valarray<double>(1);
	(*vecCT[0])(0,0) = 3;
	(*vecCT[1])(0,0) = 5;
	(*vecCT[2])(0,0) = .2;
	FT[0] = 13.6;
	std::vector<Matrix<double>* > vecCc(3);
	for (size_t i=0; i<3; i++) vecCc[i] = new Matrix<double>(1,1);
	std::valarray<double>* Fc = new std::valarray<double>(1);
	(*vecCc[0])(0,0) = .5;
	(*vecCc[1])(0,0) = .3;
	(*vecCc[2])(0,0) = 3;
	Fc[0] = 10.1;

	Variable *p = new Variable(0, 1);
	Variable *T = new Variable(1, 1);
	Variable *c = new Variable(2, 1);
	LinearComponent *comp_pp = new LinearComponent(vecCp[0], Fp);
	LinearComponent *comp_pT = new LinearComponent(vecCp[1]);
	LinearComponent *comp_pc = new LinearComponent(vecCp[2]);
	LinearComponent *comp_Tp = new LinearComponent(vecCT[0]);
	LinearComponent *comp_TT = new LinearComponent(vecCT[1], FT);
	LinearComponent *comp_Tc = new LinearComponent(vecCT[2]);
	LinearComponent *comp_cp = new LinearComponent(vecCc[0]);
	LinearComponent *comp_cT = new LinearComponent(vecCc[1]);
	LinearComponent *comp_cc = new LinearComponent(vecCc[2], Fc);

	LinearEquation *odeFlow = new LinearEquation;
	odeFlow->addVariable(*p, *comp_pp, true);
	odeFlow->addVariable(*T, *comp_pT);
	odeFlow->addVariable(*c, *comp_pc);

	LinearEquation* odeTransport = new LinearEquation;
	odeTransport->addVariable(*p, *comp_Tp);
	odeTransport->addVariable(*T, *comp_TT, true);
	odeTransport->addVariable(*c, *comp_Tc);

	LinearEquation* odeTransport2 = new LinearEquation;
	odeTransport2->addVariable(*p, *comp_cp);
	odeTransport2->addVariable(*T, *comp_cT);
	odeTransport2->addVariable(*c, *comp_cc, true);

	sysEqs.addEquation(*odeFlow);
	sysEqs.addEquation(*odeTransport);
	sysEqs.addEquation(*odeTransport2);

	ASSERT_EQ(3, sysEqs.getNumberOfVariables());
	ASSERT_EQ(3, sysEqs.getNumberOfEquations());

}

TEST(Math, SystemOfEqs1)
{
	SystemOfEquations sysEqs;
	defineProblem1(sysEqs);

	CoupledProblemConstructor probgen;
	MyFunction* f = new MyFunction(*sysEqs.getListOfVariables(), *sysEqs.getListOfEquations());

    std::map<size_t, std::valarray<double> > inactive_x;
	DenseLinearEquations eqs1;
	f->eval(inactive_x, eqs1);
	eqs1.solve();
	double *u = eqs1.getX();

	double expected[3];
	expected[0] = 1.0;
	expected[1] = 2.0;
	expected[2] = 3.0;

	ASSERT_DOUBLE_ARRAY_EQ(expected, u, 3);
}

#if 0
TEST(Math, SystemOfEqs2)
{
	SystemOfEquations sysEqs;
	defineProblem1(sysEqs);

	{
	std::vector<Variable*> active_vars1;
	active_vars1.push_back(sysEqs.getVariable(0));
	active_vars1.push_back(sysEqs.getVariable(1));
	CoupledProblemConstructor probgen;
	MyCouplingEQS* f1 = probgen.createPartitionedProblem(sysEqs, active_vars1);

	std::valarray<double> x0(1);
    std::map<size_t, std::valarray<double> > inactive_x;
    inactive_x[2] = x0;

	DenseLinearEquations eqs1;
	f1->eval(inactive_x, eqs1);
	eqs1.solve();
	double *u1 = eqs1.getX();

	double expected[2];
	expected[0] = 3.45 - 3.25 / 2.0;
	expected[1] = 3.25 / 2.0;

	ASSERT_DOUBLE_ARRAY_EQ(expected, u1, 2);
	}

	{
	std::vector<Variable*> active_vars2;
	active_vars2.push_back(sysEqs.getVariable(2));
	CoupledProblemConstructor probgen;
	MyCouplingEQS* f2 = probgen.createPartitionedProblem(sysEqs, active_vars2);

	std::valarray<double> x0(1);
    std::map<size_t, std::valarray<double> > inactive_x;
    inactive_x[0] = x0;
    inactive_x[1] = x0;

	DenseLinearEquations eqs2;
	f2->eval(inactive_x, eqs2);
	eqs2.solve();
	double *u2 = eqs2.getX();

	double expected[1];
	expected[0] = 10.1/3.0;

	ASSERT_DOUBLE_ARRAY_EQ(expected, u2, 1);
	}
}
#endif

class MyConvergenceCheck4Array
{
public:
	bool isConverged(NamedParameterSet& vars_prev, NamedParameterSet& vars_current, double eps, double &v_diff)
	{
	    for (size_t i=0; i<vars_prev.size(); i++) {
	        MyCouplingEQS::ArrayType *v_prev = vars_prev.get<MyCouplingEQS::ParameterType>(i)->getArray();
	        MyCouplingEQS::ArrayType *v_cur = vars_current.get<MyCouplingEQS::ParameterType>(i)->getArray();
	        v_diff = std::abs((*v_cur)[0] - (*v_prev)[0]);
	        if (v_diff>eps) {
	            return false;
	        }
	    }
	    return true;
	}
};

#if 0
TEST(Math, SystemOfEqs3)
{
	SystemOfEquations sysEqs;
	defineProblem1(sysEqs);

	std::vector<Variable*> active_vars1;
	active_vars1.push_back(sysEqs.getVariable(0));
	active_vars1.push_back(sysEqs.getVariable(1));
	std::vector<Variable*> active_vars2;
	active_vars2.push_back(sysEqs.getVariable(2));

	CoupledProblemConstructor probgen;
	MyCouplingEQS* eqs1 = probgen.createPartitionedProblem(sysEqs, active_vars1);
	MyCouplingEQS* eqs2 = probgen.createPartitionedProblem(sysEqs, active_vars2);

    MyCouplingEQS::ArrayType vec0(1);
    //eqs1->setOutput(0, new MyCouplingEQS::ParameterType(vec0));
    //eqs1->setOutput(1, new MyCouplingEQS::ParameterType(vec0));
    //eqs2->setOutput(2, new MyCouplingEQS::ParameterType(vec0));

	BlockJacobiMethod<MyConvergenceCheck4Array> method(1.e-4, 100);
    PartitionedProblem part1(method);
    part1.addParameter("a", *eqs1, 0);
    part1.addParameter("b", *eqs1, 1);
    part1.addParameter("c", *eqs2, 2);
    part1.connectInput("c", *eqs1, 2);
    part1.connectInput("a", *eqs2, 0);
    part1.connectInput("b", *eqs2, 1);

    ASSERT_TRUE(part1.check());
    part1.solve();

    const double epsilon = 1.e-3;
    const MyCouplingEQS::ParameterType *o1 = part1.getOutput<MyCouplingEQS::ParameterType>(part1.getParameterID("a"));
    const MyCouplingEQS::ParameterType *o2 = part1.getOutput<MyCouplingEQS::ParameterType>(part1.getParameterID("b"));
    const MyCouplingEQS::ParameterType *o3 = part1.getOutput<MyCouplingEQS::ParameterType>(part1.getParameterID("c"));
    MyCouplingEQS::ArrayType* v1 = o1->getArray();
    MyCouplingEQS::ArrayType* v2 = o2->getArray();
    MyCouplingEQS::ArrayType* v3 = o3->getArray();
    ASSERT_NEAR(1., (*v1)[0], epsilon);
    ASSERT_NEAR(2., (*v2)[0], epsilon);
    ASSERT_NEAR(3., (*v3)[0], epsilon);

}
#endif

