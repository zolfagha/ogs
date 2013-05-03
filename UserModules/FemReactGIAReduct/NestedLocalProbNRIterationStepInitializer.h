/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NestedLocalProbNRIterationStepInitializer.h
 *
 * This class is used in Newton-Raphson iterations to include the 
 * calcualation of local ODE problems
 *
 * after the file DiscreteNonlinearSolverFactory.h
 * 
 * Created on 2012-10-05 by Haibing Shao
 */

#ifndef NESTED_Local_Prob_NR_ITERATION_STEP_INITIALIZER_H
#define NESTED_Local_Prob_NR_ITERATION_STEP_INITIALIZER_H

#include "UserModules/FemKinReactGIA/NonLinearReactiveTransportTimeODELocalAssembler.h"  // TO BE CHANGED
#include "UserModules/FemKinReactGIA/NonLinearReactiveTransportJacabianLocalAssembler.h" // TO BE CHANGED
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"

template <class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>
class NestedLocalProbNRIterationStepInitializer
{
public:
    NestedLocalProbNRIterationStepInitializer(
		    NonLinearReactiveTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, T_NODAL_FUNCTION_SCALAR>* assemblerR,
            NonLinearReactiveTransportJacobianLocalAssembler<T_NODAL_FUNCTION_SCALAR, T_FUNCTION_DATA> *assemblerJ)
            : _assemblerR(assemblerR), _assemblerJ(assemblerJ) {};

	template<class T_X, class F_RESIDUALS, class F_DX>
    void pre_process(const T_X &/*dx*/, const T_X &/*x_new*/, F_RESIDUALS & f_residuals, F_DX &/*f_dx*/)
    {
		// local ode problem
		// _assemblerJ->get_function_data()->calc_nodal_xi_immob_ode( f_residuals.getTimeStepObj()->getTimeStepSize() ); 

        // update the rate calculation
        // _assemblerJ->get_function_data()->update_node_kin_reaction_rates(); 

        // calculate change of rates over change of xi_mob
		// _assemblerJ->get_function_data()->update_node_kin_reaction_drates_dxi(); 
    };

    template<class T_X, class F_RESIDUALS, class F_DX>
    void post_process(const T_X &/*dx*/, const T_X & x_new, F_RESIDUALS &/*f_residuals*/, F_DX &/*f_dx*/)
    {
        // very important is, now we have x_new, 
        // the data in the FunctionConcentrations needs to be updated, 
        // so that the in the next Newton iteration, concentrations are calculated
        // based on the new x_new values
        // _assemblerJ->get_function_data()->update_xi_mob_nodal_values( x_new ); 
    };

private:
    /**
      * pointer to assemblerR
      */ 
    NonLinearReactiveTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, T_NODAL_FUNCTION_SCALAR>* _assemblerR;

    /**
      * pointer to assemblerJ
      */ 
    NonLinearReactiveTransportJacobianLocalAssembler<T_NODAL_FUNCTION_SCALAR, T_FUNCTION_DATA>* _assemblerJ;
};

#endif  // end of ifndef