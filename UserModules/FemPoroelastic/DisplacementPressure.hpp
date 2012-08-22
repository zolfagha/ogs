

#include "logog.hpp"

#include "MeshLib/Tools/Tools.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "SolutionLib/Fem/FemVariable.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "Ogs6FemData.h"

#include "../FemDeformationTotalForm/FemLinearElasticTools.h"

template <class T1, class T2>
typename FunctionDisplacementPressure<T1,T2>::MyVariable* FunctionDisplacementPressure<T1,T2>::getDisplacementComponent(MyVariable *u_x, MyVariable* u_y, MyVariable* u_z, const std::string &var_name)
{
    if (var_name.find("_X")!=std::string::npos) {
        return u_x;
    } else if (var_name.find("_Y")!=std::string::npos) {
        return u_y;
    } else {
        return u_z;
    }
}

template <class T1, class T2>
bool FunctionDisplacementPressure<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOption<size_t>("MeshID");
    size_t time_id = option.getOption<size_t>("TimeGroupID");

    //--------------------------------------------------------------------------
    // set up mesh and FE objects
    //--------------------------------------------------------------------------
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    MeshLib::generateHigherOrderUnstrucuredMesh(*(MeshLib::UnstructuredMesh*)msh, 2);
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangianFeObjectContainer(*msh);

    //--------------------------------------------------------------------------
    // set up problem
    //--------------------------------------------------------------------------
    _problem = new MyProblemType(dis);
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];
    _problem->setTimeSteppingFunction(*tim);

    //--------------------------------------------------------------------------
    // set up variables
    //--------------------------------------------------------------------------
    // definitions
    MyVariable* u_x = _problem->addVariable("u_x", FemLib::PolynomialOrder::Quadratic);
    MyVariable* u_y = _problem->addVariable("u_y", FemLib::PolynomialOrder::Quadratic);
    MyVariable* p = _problem->addVariable("p", FemLib::PolynomialOrder::Linear);
    _var_p_id = p->getID();
    // IC
    NumLib::TXFunctionBuilder f_builder;
    MyNodalFunctionScalar* u0 = new MyNodalFunctionScalar();
    u0->initialize(*dis, u_x->getCurrentOrder(), 0);
    u_x->setIC(u0);
    u_y->setIC(u0);
    MyNodalFunctionScalar* p0 = new MyNodalFunctionScalar();
    p0->initialize(*dis, p->getCurrentOrder(), 0);
    p->setIC(p0);
    // BC
    const BaseLib::Options* opBCList = option.getSubGroup("BCList");
    for (const BaseLib::Options* opBC = opBCList->getFirstSubGroup("BC"); opBC!=0; opBC = opBCList->getNextSubGroup())
    {
        std::string var_name = opBC->getOption("Variable");
        std::string geo_type = opBC->getOption("GeometryType");
        std::string geo_name = opBC->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = femData->geo->searchGeoByName(femData->geo_unique_name, geo_type, geo_name);
        std::string dis_name = opBC->getOption("DistributionType");
        double dis_v = opBC->getOption<double>("DistributionValue");
        NumLib::ITXFunction* f_bc =  f_builder.create(dis_name, dis_v);
        if (var_name.find(this->getOutputParameterName(Displacement))!=std::string::npos) {
            getDisplacementComponent(u_x, u_y, 0, var_name)->addDirichletBC(new SolutionLib::FemDirichletBC(msh, geo_obj, f_bc));
        } else if (var_name.find(this->getOutputParameterName(Pressure))!=std::string::npos) {
            p->addDirichletBC(new SolutionLib::FemDirichletBC(msh, geo_obj, f_bc));
        }
    }
    // ST
    const BaseLib::Options* opSTList = option.getSubGroup("STList");
    for (const BaseLib::Options* opST = opSTList->getFirstSubGroup("ST"); opST!=0; opST = opSTList->getNextSubGroup())
    {
        std::string var_name = opST->getOption("Variable");
        std::string geo_type = opST->getOption("GeometryType");
        std::string geo_name = opST->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = femData->geo->searchGeoByName(femData->geo_unique_name, geo_type, geo_name);
        std::string st_type = opST->getOption("STType");
        std::string dis_name = opST->getOption("DistributionType");
        double dis_v = opST->getOption<double>("DistributionValue");
        if (st_type.compare("NEUMANN")==0) {
            dis_v *= -1; //TODO
        }
        NumLib::ITXFunction* f_st =  f_builder.create(dis_name, dis_v);
        if (f_st!=NULL) {
            SolutionLib::IFemNeumannBC *femSt = 0;
            if (st_type.compare("NEUMANN")==0) {
                femSt = new SolutionLib::FemNeumannBC(msh, _feObjects, geo_obj, f_st);
            } else if (st_type.compare("SOURCESINK")==0) {
                femSt = new SolutionLib::FemSourceTerm(msh, geo_obj, f_st);
            }
            if (var_name.find(this->getOutputParameterName(Displacement))!=std::string::npos) {
                getDisplacementComponent(u_x, u_y, 0, var_name)->addNeumannBC(femSt);
            } else if (var_name.find(this->getOutputParameterName(Pressure))!=std::string::npos) {
                p->addNeumannBC(femSt);
            }
        } else {
            WARN("Distribution type %s is specified but not found. Ignore this ST.", dis_name.c_str());
        }
    }

    //--------------------------------------------------------------------------
    // set up equations
    //--------------------------------------------------------------------------
    MyEquationType* eqs = _problem->createEquation();
    std::vector<size_t> vec_orders;
    for (size_t i=0; i<_problem->getNumberOfVariables(); i++)
        vec_orders.push_back(_problem->getVariable(i)->getCurrentOrder());

    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects, _problem->getNumberOfVariables(), vec_orders);
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(_feObjects);
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(_feObjects);
    eqs->initialize(linear_assembler, r_assembler, j_eqs);


    //--------------------------------------------------------------------------
    // set up this solution algorithm
    //--------------------------------------------------------------------------
    _solution = new MySolutionType(dis, _problem);
    typename MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);
    _solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);

    //--------------------------------------------------------------------------
    // set up data for output
    //--------------------------------------------------------------------------
    // create u variable which is vector
    NumLib::LocalVector tmp_u0(3);
    tmp_u0 *= .0;
    _displacement = new MyNodalFunctionVector();
    _displacement->initialize(*dis, u0->getOrder(), tmp_u0);
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        _displacement->getValue(i)(0) = _solution->getCurrentSolution(0)->getValue(i);
        _displacement->getValue(i)(1) = _solution->getCurrentSolution(1)->getValue(i);
    }
    for (size_t i=0; i<2; i++) {
        _vec_u_components.push_back(new NodalPointScalarWrapper(_displacement, i));
    }

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Displacement), OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i=0; i<_vec_u_components.size(); i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Displacement) + getDisplacementComponentPostfix(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    OutputVariableInfo outP(this->getOutputParameterName(Pressure), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(_var_p_id));
    femData->outController.setOutput(outP.name, outP);

    //--------------------------------------------------------------------------
    // set initial output
    //--------------------------------------------------------------------------
    // initial output parameter
    this->setOutput(Displacement, _displacement);
    this->setOutput(Pressure, _solution->getCurrentSolution(_var_p_id));


    return true;
}

template <class T1, class T2>
void FunctionDisplacementPressure<T1,T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        _displacement->getValue(i)(0) = _solution->getCurrentSolution(0)->getValue(i);
        _displacement->getValue(i)(1) = _solution->getCurrentSolution(1)->getValue(i);
    }
    setOutput(Displacement, _displacement);
    this->setOutput(Pressure, _solution->getCurrentSolution(_var_p_id));

    //calculateStressStrain();
}

template <class T1, class T2>
void FunctionDisplacementPressure<T1,T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Displacement), OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i=0; i<_vec_u_components.size(); i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Displacement) + getDisplacementComponentPostfix(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    OutputVariableInfo outP(this->getOutputParameterName(Pressure), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(_var_p_id));
    femData->outController.setOutput(outP.name, outP);
};

