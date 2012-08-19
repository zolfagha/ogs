
OGS_ADD_PROCESS_SYS_SOLVER(GROUNDWATER_FLOW, FunctionHead, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS_SOLVER(GROUNDWATER_FLOW2, FunctionHead, DiscreteLib::SerialNodeDdcSharedDiscreteSystem, MathLib::DenseLinearEquation);
OGS_ADD_PROCESS_SYS(ELEMENT_VELOCITY, FunctionElementVelocity, DiscreteLib::DiscreteSystem);
OGS_ADD_PROCESS_SYS_SOLVER(MASS_TRANSPORT, FunctionConcentration, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS_SOLVER(DEFORMATION, FunctionDisplacement, DiscreteLib::DiscreteSystem, MathLib::LisLinearEquation);
OGS_ADD_PROCESS_SYS(ELEMENT_STRESS_STRAIN, FunctionElementStressStrain, DiscreteLib::DiscreteSystem);
OGS_ADD_PROCESS_SYS(NODAL_STRESS_STRAIN, FunctionNodalStressStrain, DiscreteLib::DiscreteSystem);
