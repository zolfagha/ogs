
OGS_ADD_PROCESS_SYS_SOLVER(LIQUID_FLOW, FunctionLiquidPressure, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS_SOLVER(GROUNDWATER_FLOW, FunctionHead, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
////OGS_ADD_PROCESS_SYS_SOLVER(GROUNDWATER_FLOW, FunctionHead, DiscreteLib::OMPDiscreteSystem, MathLib::LisLinearEquation);
//OGS_ADD_PROCESS_SYS(HEAD_TO_ELEMENT_VELOCITY, FunctionHeadToElementVelocity, DiscreteLib::DiscreteSystem);
OGS_ADD_PROCESS_SYS(PRESSURE_TO_ELEMENT_VELOCITY, FunctionPressureToElementVelocity, DiscreteLib::DiscreteSystem);
//OGS_ADD_PROCESS_SYS(PRESSURE_TO_HEAD, FunctionPressureToHead, DiscreteLib::DiscreteSystem);
//OGS_ADD_PROCESS_SYS_SOLVER(MASS_TRANSPORT, FunctionConcentration, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
////OGS_ADD_PROCESS_SYS_SOLVER(MASS_TRANSPORT, FunctionConcentration, DiscreteLib::DiscreteSystem, MathLib::DenseLinearEquation);
//OGS_ADD_PROCESS_SYS_SOLVER(DEFORMATION, FunctionDisplacement, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
//OGS_ADD_PROCESS_SYS(ELEMENT_STRESS_STRAIN, FunctionElementStressStrain, DiscreteLib::DiscreteSystem);
//OGS_ADD_PROCESS_SYS(NODAL_STRESS_STRAIN, FunctionNodalStressStrain, DiscreteLib::DiscreteSystem);
//OGS_ADD_PROCESS_SYS_SOLVER(DEFORMATION_FLOW, FunctionDisplacementPressure, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS_SOLVER(XFEM_EXAMPLE_CRACK1, xfem::FunctionDisplacement, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
