
OGS_ADD_PROCESS_SYS_SOLVER(GROUNDWATER_FLOW, FunctionHead, DiscreteLib::DiscreteSystem, MathLib::DenseLinearEquation);
//OGS_ADD_PROCESS_SYS_SOLVER(GROUNDWATER_FLOW, FunctionHead, DiscreteLib::OMPDiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS(ELEMENT_VELOCITY, FunctionHeadToElementVelocity, DiscreteLib::DiscreteSystem);
//OGS_ADD_PROCESS_SYS_SOLVER(MASS_TRANSPORT, FunctionConcentration,  DiscreteLib::DiscreteSystem, MathLib::DenseLinearEquation);
//OGS_ADD_PROCESS_SYS_SOLVER(KIN_REACT_GIA,  FunctionConcentrations, DiscreteLib::DiscreteSystem, MathLib::DenseLinearEquation  );
// OGS_ADD_PROCESS_SYS_SOLVER(DEFORMATION, FunctionDisplacement, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
//OGS_ADD_PROCESS_SYS_SOLVER(DEFORMATION, FunctionDisplacement, DiscreteLib::DiscreteSystem, MathLib::DenseLinearEquation);
//OGS_ADD_PROCESS_SYS(ELEMENT_STRESS_STRAIN, FunctionElementStressStrain, DiscreteLib::DiscreteSystem);
//OGS_ADD_PROCESS_SYS(NODAL_STRESS_STRAIN, FunctionNodalStressStrain, DiscreteLib::DiscreteSystem);
