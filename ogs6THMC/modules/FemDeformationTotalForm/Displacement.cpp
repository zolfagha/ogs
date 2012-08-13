
#include "Displacement.h"

#include "logog.hpp"

#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"

#include "Ogs6FemData.h"

#include "FemLinearElasticTools.h"

SolutionLib::FemVariable* getDisplacementComponent(SolutionLib::FemVariable *u_x, SolutionLib::FemVariable* u_y, SolutionLib::FemVariable* u_z, const std::string &var_name)
{
    if (var_name.compare("DISPLACEMENT_X1")==0) {
        return u_x;
    } else if (var_name.compare("DISPLACEMENT_Y1")==0) {
        return u_y;
    } else {
        return u_z;
    }
}

bool FunctionDisplacement::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOption<size_t>("MeshID");
    size_t time_id = option.getOption<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    DiscreteLib::DiscreteSystem* dis = femData->list_dis_sys[msh_id];
    MeshLib::IMesh* msh = dis->getMesh();
    _feObjects = new FemLib::LagrangianFeObjectContainer(*msh);

    // equations
    MyEquationType::LinearAssemblerType* linear_assembler = new MyEquationType::LinearAssemblerType(_feObjects);
    MyEquationType::ResidualAssemblerType* r_assembler = new MyEquationType::ResidualAssemblerType(_feObjects);
    MyEquationType::JacobianAssemblerType* j_eqs = new MyEquationType::JacobianAssemblerType(_feObjects);
    MyEquationType* eqs = new  MyEquationType(linear_assembler, r_assembler, j_eqs);

    // set up problem
    _problem = new MyProblemType(dis);
    _problem->setEquation(eqs);
    _problem->setTimeSteppingFunction(*tim);
    // set up variable
    SolutionLib::FemVariable* u_x = _problem->addVariable("u_x");
    SolutionLib::FemVariable* u_y = _problem->addVariable("u_y");
    // IC
    NumLib::TXFunctionBuilder f_builder;
    FemLib::FemNodalFunctionScalar* u0 = new FemLib::FemNodalFunctionScalar(*dis, FemLib::PolynomialOrder::Linear, 0);
    u_x->setIC(u0);
    u_y->setIC(u0);
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
        getDisplacementComponent(u_x, u_y, 0, var_name)->addDirichletBC(new SolutionLib::FemDirichletBC(msh, geo_obj, f_bc));
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
        NumLib::ITXFunction* f_st =  f_builder.create(dis_name, dis_v);
        if (f_st!=NULL) {
            SolutionLib::IFemNeumannBC *femSt = 0;
            if (st_type.compare("NEUMANN")==0) {
                femSt = new SolutionLib::FemNeumannBC(msh, _feObjects, geo_obj, f_st);
            } else if (st_type.compare("SOURCESINK")==0) {
                femSt = new SolutionLib::FemSourceTerm(msh, geo_obj, f_st);
            }
            getDisplacementComponent(u_x, u_y, 0, var_name)->addNeumannBC(femSt);
        } else {
            WARN("Distribution type %s is specified but not found. Ignore this ST.", dis_name.c_str());
        }
    }

    // set up solution
    _solution = new MySolutionType(dis, _problem);
    MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);

    // set initial output
    OutputVariableInfo var("DISPLACEMENT", OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);

    // initial output parameter
    this->setOutput(Displacement, _problem->getVariable(0)->getIC());

    NumLib::LocalVector tmp_u0(3);
    tmp_u0 *= .0;
    _displacement = new FemLib::FemNodalFunctionVector(*dis, FemLib::PolynomialOrder::Quadratic, tmp_u0);
//    _strain = new FemLib::FEMIntegrationPointFunctionVector(*dis);
//    _stress = new FemLib::FEMIntegrationPointFunctionVector(*dis);

    return true;
}

void FunctionDisplacement::updateOutputParameter(const NumLib::TimeStep &time)
{
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        _displacement->getValue(i)(0) = _solution->getCurrentSolution(0)->getValue(i);
        _displacement->getValue(i)(1) = _solution->getCurrentSolution(1)->getValue(i);
    }
    setOutput(Displacement, _displacement);

    //calculateStressStrain();
}

void FunctionDisplacement::output(const NumLib::TimeStep &time)
{
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var("DISPLACEMENT", OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);
};

//void FunctionDisplacement::calculateStressStrain()
//{
//    const MeshLib::IMesh *msh = _problem->getMesh();
//    FemLib::FemNodalFunctionScalar *ux = _solution->getCurrentSolution(0);
//    FemLib::FemNodalFunctionScalar *uy = _solution->getCurrentSolution(1);
//    FemLib::FEMIntegrationPointFunctionVector* strain = _strain;
//    FemLib::FEMIntegrationPointFunctionVector* stress = _stress;
//    FemLib::LagrangianFeObjectContainer* feObjects = ux->getFeObjectContainer();
//
//    const size_t dim = msh->getDimension();
//    const size_t n_strain_components = (dim==2 ? 4 : 6);
//
//    Ogs6FemData* femData = Ogs6FemData::getInstance();
//
//    //calculate strain, stress
//    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++)
//    {
//        // element setup
//        MeshLib::IElement* e = msh->getElemenet(i_e);
//        const size_t nnodes = e->getNumberOfNodes();
//        FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
//        FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
//        const size_t n_gp = integral->getNumberOfSamplingPoints();
//        strain->setNumberOfIntegationPoints(i_e, n_gp);
//        stress->setNumberOfIntegationPoints(i_e, n_gp);
//        const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e->getID());
//
//        MaterialLib::Solid *solidphase = femData->list_solid[e->getGroupID()];
//        // set D
//        NumLib::LocalMatrix matD = NumLib::LocalMatrix::Zero(n_strain_components, n_strain_components);
//        //matD *= .0;
//        NumLib::LocalMatrix nv(1,1);
//        NumLib::LocalMatrix E(1,1);
//        solidphase->poisson_ratio->eval(e_pos, nv);
//        solidphase->Youngs_modulus->eval(e_pos, E);
//        double Lambda, G, K;
//        MaterialLib::calculateLameConstant(nv(0,0), E(0,0), Lambda, G, K);
//        MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, matD);
//
//        // local u
//        NumLib::LocalVector local_u(dim*nnodes);
//        for (size_t j=0; j<nnodes; j++) {
//            const size_t node_id = e->getNodeID(j);
//            local_u[j*dim] = ux->getValue(node_id);
//            local_u[j*dim+1] = uy->getValue(node_id);
//        }
//
//        // for each integration points
//        NumLib::LocalMatrix matB(n_strain_components, nnodes*dim);
//        NumLib::LocalMatrix matN(dim, nnodes*dim);
//        double r[3] = {};
//        double x[3] = {};
//        for (size_t ip=0; ip<n_gp; ip++) {
//            integral->getSamplingPoint(ip, r);
//            fe->computeBasisFunctions(r);
//            NumLib::LocalMatrix &N = *fe->getBasisFunction();
//            const NumLib::LocalMatrix &dN = *fe->getGradBasisFunction();
//            fe->getRealCoordinates(x);
//
//            // set N,B
//            setNu_Matrix(dim, nnodes, N, matN);
//            setB_Matrix(dim, nnodes, dN, matB);
//
//            // strain
//            NumLib::LocalVector gp_strain(n_strain_components);
//            gp_strain.noalias() = matB * local_u;
//            strain->setIntegrationPointValue(i_e, ip, gp_strain);
//
//            // stress
//            NumLib::LocalVector gp_stress(n_strain_components);
//            gp_stress = matD * gp_strain;
//            stress->setIntegrationPointValue(i_e, ip, gp_stress);
//
////                std::cout << "strain=\n" << gp_strain << std::endl;
////                std::cout << "D=\n" << matD << std::endl;
////                std::cout << "stress=\n" << gp_stress << std::endl;
//        }
//    }
//
//    setOutput(Strain, strain);
//    setOutput(Stress, stress);
//}
