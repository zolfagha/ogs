
#include <gtest/gtest.h>

#include "BaseLib/Options.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "NumLib/Coupling/PartitionedProblem.h"
#include "NumLib/Coupling/Algorithm/BlockJacobiMethod.h"
#include "NumLib/SystemOfEquations/SystemOfEquations.h"
#include "NumLib/SystemOfEquations/CoupledProblemConstructor.h"
#include "NumLib/SystemOfEquations/CouplingStrucutreBuilder4SysEqs.h"

#include "TestUtil.h"

using namespace MathLib;
using namespace NumLib;


void defineProblem1(SystemOfEquations &sysEqs)
{
	// 2a + 2b + .3c = 6.9
	// 3a + 5b + .2c = 13.6
	// .5a + .3b + 3c = 10.1
	// A. a=1, b=2, c=3
	std::vector<LocalMatrix* > vecCp(3);
	for (size_t i=0; i<3; i++) vecCp[i] = new LocalMatrix(1,1);
	LocalVector* Fp = new LocalVector(1);
	(*vecCp[0])(0,0) = 2;
	(*vecCp[1])(0,0) = 2;
	(*vecCp[2])(0,0) = .3;
	(*Fp)[0] = 6.9;
	std::vector<LocalMatrix* > vecCT(3);
	for (size_t i=0; i<3; i++) vecCT[i] = new LocalMatrix(1,1);
	LocalVector* FT = new LocalVector(1);
	(*vecCT[0])(0,0) = 3;
	(*vecCT[1])(0,0) = 5;
	(*vecCT[2])(0,0) = .2;
	(*FT)[0] = 13.6;
	std::vector<LocalMatrix* > vecCc(3);
	for (size_t i=0; i<3; i++) vecCc[i] = new LocalMatrix(1,1);
	LocalVector* Fc = new LocalVector(1);
	(*vecCc[0])(0,0) = .5;
	(*vecCc[1])(0,0) = .3;
	(*vecCc[2])(0,0) = 3;
	(*Fc)[0] = 10.1;

	Variable *p = new Variable(0, 1, "p");
	Variable *T = new Variable(1, 1, "T");
	Variable *c = new Variable(2, 1, "c");
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

class MyConvergenceCheck4Array : public NumLib::IConvergenceCheck
{
public:
    bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff)
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

BaseLib::Options* defineCouplingM3()
{
    BaseLib::Options* options = new BaseLib::Options();
    BaseLib::Options* coupling = options->addSubGroup("coupling");
    BaseLib::Options* M1 = coupling->addSubGroup("M1");
	std::vector<std::string> var;
	var.push_back("p");
	var.push_back("T");
	var.push_back("c");
    M1->addOptionAsArray("variable", var);

    return options;
};

BaseLib::Options* defineCouplingP1_M2M1()
{
    BaseLib::Options* options = new BaseLib::Options();
    BaseLib::Options* coupling = options->addSubGroup("coupling");
    BaseLib::Options* P1 = coupling->addSubGroup("P");
    {
        //P1->addOption("name", "P1");
        P1->addOption("algorithm", "Gauss");
        P1->addOption("convergence", "MyConvergenceCheck");
        P1->addOptionAsNum("max_itr", 100);
        P1->addOptionAsNum("epsilon", 1.e-4);
    }
    BaseLib::Options* P1_sub = P1->addSubGroup("problems");
    BaseLib::Options* M1 = P1_sub->addSubGroup("M1");
    BaseLib::Options* M2 = P1_sub->addSubGroup("M2");
    std::vector<std::string> var1;
    var1.push_back("p");
    var1.push_back("T");
    M1->addOptionAsArray("variable", var1);
    std::vector<std::string> var2;
    var2.push_back("c");
    M2->addOptionAsArray("variable", var2);

    return options;
};

BaseLib::Options* defineCouplingP1_3M1()
{
    BaseLib::Options* options = new BaseLib::Options();
    BaseLib::Options* coupling = options->addSubGroup("coupling");
    BaseLib::Options* P1 = coupling->addSubGroup("P");
    {
        //P1->addOption("name", "P1");
        P1->addOption("algorithm", "Gauss");
        P1->addOption("convergence", "MyConvergenceCheck");
        P1->addOptionAsNum("max_itr", 100);
        P1->addOptionAsNum("epsilon", 1.e-4);
    }
    BaseLib::Options* P1_sub = P1->addSubGroup("problems");
    BaseLib::Options* M1 = P1_sub->addSubGroup("M1");
    BaseLib::Options* M2 = P1_sub->addSubGroup("M2");
    BaseLib::Options* M3 = P1_sub->addSubGroup("M3");

    std::vector<std::string> var1;
    var1.push_back("p");
    M1->addOptionAsArray("variable", var1);
    std::vector<std::string> var2;
    var2.push_back("T");
    M2->addOptionAsArray("variable", var2);
    std::vector<std::string> var3;
    var3.push_back("c");
    M3->addOptionAsArray("variable", var3);

    return options;
};

BaseLib::Options* defineCouplingP1_P2M1()
{
    BaseLib::Options* options = new BaseLib::Options();
    BaseLib::Options* coupling = options->addSubGroup("coupling");
    BaseLib::Options* P2 = coupling->addSubGroup("P");
    {
        //P1->addOption("name", "P1");
        P2->addOption("algorithm", "Gauss");
        P2->addOption("convergence", "MyConvergenceCheck");
        P2->addOptionAsNum("max_itr", 100);
        P2->addOptionAsNum("epsilon", 1.e-4);
    }
    BaseLib::Options* P2_sub = P2->addSubGroup("problems");
    BaseLib::Options* P1 = P2_sub->addSubGroup("P");
    {
        //P1->addOption("name", "P1");
        P1->addOption("algorithm", "Gauss");
        P1->addOption("convergence", "MyConvergenceCheck");
        P1->addOptionAsNum("max_itr", 100);
        P1->addOptionAsNum("epsilon", 1.e-4);
    }
    BaseLib::Options* P1_sub = P1->addSubGroup("problems");
    BaseLib::Options* M1 = P1_sub->addSubGroup("M1");
    BaseLib::Options* M2 = P1_sub->addSubGroup("M2");
    BaseLib::Options* M3 = P2_sub->addSubGroup("M3");

    std::vector<std::string> var1;
    var1.push_back("p");
    M1->addOptionAsArray("variable", var1);
    std::vector<std::string> var2;
    var2.push_back("T");
    M2->addOptionAsArray("variable", var2);
    std::vector<std::string> var3;
    var3.push_back("c");
    M3->addOptionAsArray("variable", var3);

    return options;
};

TEST(Math, SystemOfEqs_M1)
{
    SystemOfEquations sysEqs;
    defineProblem1(sysEqs);
    //initial value
    MyCouplingEQS::ArrayType vec0(1);
    vec0 *= .0;
    std::vector<MyCouplingEQS::ArrayType*> ini_para(sysEqs.getNumberOfVariables(), &vec0);
    // grouping variables
    std::vector<std::vector<Variable*> > list_active_vars;
    list_active_vars.resize(1);
    list_active_vars[0].push_back(sysEqs.getVariable(0));
    list_active_vars[0].push_back(sysEqs.getVariable(1));
    list_active_vars[0].push_back(sysEqs.getVariable(2));
    std::vector<IPartitionedAlgorithm* > list_part_alg;

    // create coupling structure
    CoupledProblemConstructor probgen;
    std::vector<MyCouplingEQS*> list_sub_problem;
    ICoupledSystem *part1 = probgen.build(sysEqs, list_active_vars, list_part_alg, list_sub_problem);
    // configure
    for (size_t i=0; i<list_sub_problem.size(); i++) {
        list_sub_problem[i]->setInitial(ini_para);
    }

    //solve
    ASSERT_TRUE(part1->check());
    part1->solve();

    const double epsilon = 1.e-3;
    const MyCouplingEQS::ParameterType *o1 = (const MyCouplingEQS::ParameterType*)part1->getOutput(part1->getOutputParameterID("p"));
    const MyCouplingEQS::ParameterType *o2 = (const MyCouplingEQS::ParameterType*)part1->getOutput(part1->getOutputParameterID("T"));
    const MyCouplingEQS::ParameterType *o3 = (const MyCouplingEQS::ParameterType*)part1->getOutput(part1->getOutputParameterID("c"));
    MyCouplingEQS::ArrayType* v1 = o1->getArray();
    MyCouplingEQS::ArrayType* v2 = o2->getArray();
    MyCouplingEQS::ArrayType* v3 = o3->getArray();
    ASSERT_NEAR(1., (*v1)[0], epsilon);
    ASSERT_NEAR(2., (*v2)[0], epsilon);
    ASSERT_NEAR(3., (*v3)[0], epsilon);
}

TEST(Math, SystemOfEqs_P2)
{
	SystemOfEquations sysEqs;
	defineProblem1(sysEqs);
    //initial value
    MyCouplingEQS::ArrayType vec0(1);
    vec0 *= .0;
    std::vector<MyCouplingEQS::ArrayType*> ini_para(sysEqs.getNumberOfVariables(), &vec0);
    //coupling parameter
    MyConvergenceCheck4Array checker;
    BlockJacobiMethod method(checker, 1.e-4, 100);
    // grouping variables
    std::vector<std::vector<Variable*> > list_active_vars;
    list_active_vars.resize(2);
	list_active_vars[0].push_back(sysEqs.getVariable(0));
	list_active_vars[0].push_back(sysEqs.getVariable(1));
	list_active_vars[1].push_back(sysEqs.getVariable(2));
    std::vector<IPartitionedAlgorithm* > list_part_alg;
    list_part_alg.push_back(&method);

    // create coupling structure
	CoupledProblemConstructor probgen;
    std::vector<MyCouplingEQS*> list_sub_problem;
    ICoupledSystem *part1 = probgen.build(sysEqs, list_active_vars, list_part_alg, list_sub_problem);
    // configure
    for (size_t i=0; i<list_sub_problem.size(); i++) {
        list_sub_problem[i]->setInitial(ini_para);
    }

    //solve
    ASSERT_TRUE(part1->check());
    part1->solve();

    const double epsilon = 1.e-3;
    const MyCouplingEQS::ParameterType *o1 = (const MyCouplingEQS::ParameterType*)part1->getOutput(part1->getOutputParameterID("p"));
    const MyCouplingEQS::ParameterType *o2 = (const MyCouplingEQS::ParameterType*)part1->getOutput(part1->getOutputParameterID("T"));
    const MyCouplingEQS::ParameterType *o3 = (const MyCouplingEQS::ParameterType*)part1->getOutput(part1->getOutputParameterID("c"));
    MyCouplingEQS::ArrayType* v1 = o1->getArray();
    MyCouplingEQS::ArrayType* v2 = o2->getArray();
    MyCouplingEQS::ArrayType* v3 = o3->getArray();
    ASSERT_NEAR(1., (*v1)[0], epsilon);
    ASSERT_NEAR(2., (*v2)[0], epsilon);
    ASSERT_NEAR(3., (*v3)[0], epsilon);
}

TEST(Math, SystemOfEqs_P3)
{
    SystemOfEquations sysEqs;
    defineProblem1(sysEqs);
    //initial value
    MyCouplingEQS::ArrayType vec0(1);
    vec0 *= .0;
    std::vector<MyCouplingEQS::ArrayType*> ini_para(sysEqs.getNumberOfVariables(), &vec0);
    //coupling parameter
    MyConvergenceCheck4Array checker;
    BlockJacobiMethod method(checker, 1.e-4, 100);
    // grouping variables
    std::vector<std::vector<Variable*> > list_active_vars;
    list_active_vars.resize(3);
    list_active_vars[0].push_back(sysEqs.getVariable(0));
    list_active_vars[1].push_back(sysEqs.getVariable(1));
    list_active_vars[2].push_back(sysEqs.getVariable(2));
    std::vector<IPartitionedAlgorithm* > list_part_alg;
    list_part_alg.push_back(&method);

    // create coupling structure
    CoupledProblemConstructor probgen;
    std::vector<MyCouplingEQS*> list_sub_problem;
    ICoupledSystem *part1 = probgen.build(sysEqs, list_active_vars, list_part_alg, list_sub_problem);
    // configure
    for (size_t i=0; i<list_sub_problem.size(); i++) {
        list_sub_problem[i]->setInitial(ini_para);
    }

    //solve
    ASSERT_TRUE(part1->check());
    part1->solve();

    const double epsilon = 1.e-3;
    const MyCouplingEQS::ParameterType *o1 = (const MyCouplingEQS::ParameterType*)part1->getOutput(part1->getOutputParameterID("p"));
    const MyCouplingEQS::ParameterType *o2 = (const MyCouplingEQS::ParameterType*)part1->getOutput(part1->getOutputParameterID("T"));
    const MyCouplingEQS::ParameterType *o3 = (const MyCouplingEQS::ParameterType*)part1->getOutput(part1->getOutputParameterID("c"));
    MyCouplingEQS::ArrayType* v1 = o1->getArray();
    MyCouplingEQS::ArrayType* v2 = o2->getArray();
    MyCouplingEQS::ArrayType* v3 = o3->getArray();
    ASSERT_NEAR(1., (*v1)[0], epsilon);
    ASSERT_NEAR(2., (*v2)[0], epsilon);
    ASSERT_NEAR(3., (*v3)[0], epsilon);
}

class MyConvergenceChecker4ArrayFactory
{
public:
	IConvergenceCheck* create(const std::string &)
	{
		return new MyConvergenceCheck4Array();
	};
};


TEST(Math, SystemOfEqs_AutoM3)
{
    SystemOfEquations sysEqs;
    defineProblem1(sysEqs);
    //initial value
    MyCouplingEQS::ArrayType vec0(1);
    vec0 *= .0;
    std::vector<MyCouplingEQS::ArrayType*> ini_para(sysEqs.getNumberOfVariables(), &vec0);

    BaseLib::Options* option = defineCouplingM3();
    std::vector<MyCouplingEQS*> list_sub_problem;
    CoupledProblemFactory eqs_fac(sysEqs, ini_para, list_sub_problem);
    MyConvergenceChecker4ArrayFactory checker;
    CouplingStrucutreBuilder4SysEqs cpl_builder;
    ICoupledSystem* cpl_sys = cpl_builder.build(option, eqs_fac, checker);

    //for (size_t i=0; i<list_sub_problem.size(); i++) {
    //    list_sub_problem[i]->setInitial(ini_para);
    //}

    //solve
    ASSERT_TRUE(cpl_sys->check());
    cpl_sys->solve();

    const double epsilon = 1.e-3;
    const MyCouplingEQS::ParameterType *o1 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("p"));
    const MyCouplingEQS::ParameterType *o2 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("T"));
    const MyCouplingEQS::ParameterType *o3 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("c"));
    MyCouplingEQS::ArrayType* v1 = o1->getArray();
    MyCouplingEQS::ArrayType* v2 = o2->getArray();
    MyCouplingEQS::ArrayType* v3 = o3->getArray();
    ASSERT_NEAR(1., (*v1)[0], epsilon);
    ASSERT_NEAR(2., (*v2)[0], epsilon);
    ASSERT_NEAR(3., (*v3)[0], epsilon);
}

TEST(Math, SystemOfEqs_AutoP1_M2M1)
{
    SystemOfEquations sysEqs;
    defineProblem1(sysEqs);
    //initial value
    MyCouplingEQS::ArrayType vec0(1);
    vec0 *= .0;
    std::vector<MyCouplingEQS::ArrayType*> ini_para(sysEqs.getNumberOfVariables(), &vec0);

    BaseLib::Options* option = defineCouplingP1_M2M1();
    std::vector<MyCouplingEQS*> list_sub_problem;
    CoupledProblemFactory eqs_fac(sysEqs, ini_para, list_sub_problem);
    MyConvergenceChecker4ArrayFactory checker;
    CouplingStrucutreBuilder4SysEqs cpl_builder;
    ICoupledSystem* cpl_sys = cpl_builder.build(option, eqs_fac, checker);

    //for (size_t i=0; i<list_sub_problem.size(); i++) {
    //    list_sub_problem[i]->setInitial(ini_para);
    //}

    //solve
    ASSERT_TRUE(cpl_sys->check());
    cpl_sys->solve();

    const double epsilon = 1.e-3;
    const MyCouplingEQS::ParameterType *o1 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("p"));
    const MyCouplingEQS::ParameterType *o2 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("T"));
    const MyCouplingEQS::ParameterType *o3 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("c"));
    MyCouplingEQS::ArrayType* v1 = o1->getArray();
    MyCouplingEQS::ArrayType* v2 = o2->getArray();
    MyCouplingEQS::ArrayType* v3 = o3->getArray();
    ASSERT_NEAR(1., (*v1)[0], epsilon);
    ASSERT_NEAR(2., (*v2)[0], epsilon);
    ASSERT_NEAR(3., (*v3)[0], epsilon);
}

TEST(Math, SystemOfEqs_AutoP1_3M1)
{
    SystemOfEquations sysEqs;
    defineProblem1(sysEqs);
    //initial value
    MyCouplingEQS::ArrayType vec0(1);
    vec0 *= .0;
    std::vector<MyCouplingEQS::ArrayType*> ini_para(sysEqs.getNumberOfVariables(), &vec0);

    BaseLib::Options* option = defineCouplingP1_3M1();
    std::vector<MyCouplingEQS*> list_sub_problem;
    CoupledProblemFactory eqs_fac(sysEqs, ini_para, list_sub_problem);
    MyConvergenceChecker4ArrayFactory checker;
    CouplingStrucutreBuilder4SysEqs cpl_builder;
    ICoupledSystem* cpl_sys = cpl_builder.build(option, eqs_fac, checker);

    //for (size_t i=0; i<list_sub_problem.size(); i++) {
    //    list_sub_problem[i]->setInitial(ini_para);
    //}

    //solve
    ASSERT_TRUE(cpl_sys->check());
    cpl_sys->solve();

    const double epsilon = 1.e-3;
    const MyCouplingEQS::ParameterType *o1 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("p"));
    const MyCouplingEQS::ParameterType *o2 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("T"));
    const MyCouplingEQS::ParameterType *o3 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("c"));
    MyCouplingEQS::ArrayType* v1 = o1->getArray();
    MyCouplingEQS::ArrayType* v2 = o2->getArray();
    MyCouplingEQS::ArrayType* v3 = o3->getArray();
    ASSERT_NEAR(1., (*v1)[0], epsilon);
    ASSERT_NEAR(2., (*v2)[0], epsilon);
    ASSERT_NEAR(3., (*v3)[0], epsilon);
}

TEST(Math, SystemOfEqs_AutoP1_P2M1)
{
    SystemOfEquations sysEqs;
    defineProblem1(sysEqs);
    //initial value
    MyCouplingEQS::ArrayType vec0(1);
    vec0 *= .0;
    std::vector<MyCouplingEQS::ArrayType*> ini_para(sysEqs.getNumberOfVariables(), &vec0);

    BaseLib::Options* option = defineCouplingP1_P2M1();
    std::vector<MyCouplingEQS*> list_sub_problem;
    CoupledProblemFactory eqs_fac(sysEqs, ini_para, list_sub_problem);
    MyConvergenceChecker4ArrayFactory checker;
    CouplingStrucutreBuilder4SysEqs cpl_builder;
    ICoupledSystem* cpl_sys = cpl_builder.build(option, eqs_fac, checker);

    //for (size_t i=0; i<list_sub_problem.size(); i++) {
    //    list_sub_problem[i]->setInitial(ini_para);
    //}

    //solve
    ASSERT_TRUE(cpl_sys->check());
    cpl_sys->solve();

    const double epsilon = 1.e-3;
    const MyCouplingEQS::ParameterType *o1 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("p"));
    const MyCouplingEQS::ParameterType *o2 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("T"));
    const MyCouplingEQS::ParameterType *o3 = (const MyCouplingEQS::ParameterType*)cpl_sys->getOutput(cpl_sys->getOutputParameterID("c"));
    MyCouplingEQS::ArrayType* v1 = o1->getArray();
    MyCouplingEQS::ArrayType* v2 = o2->getArray();
    MyCouplingEQS::ArrayType* v3 = o3->getArray();
    ASSERT_NEAR(1., (*v1)[0], epsilon);
    ASSERT_NEAR(2., (*v2)[0], epsilon);
    ASSERT_NEAR(3., (*v3)[0], epsilon);
}
