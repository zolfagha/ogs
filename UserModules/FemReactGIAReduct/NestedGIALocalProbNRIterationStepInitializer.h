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

#ifndef NESTED_GIA_Local_Prob_NR_ITERATION_STEP_INITIALIZER_H
#define NESTED_GIA_Local_Prob_NR_ITERATION_STEP_INITIALIZER_H

//#include "NonLinearGIATimeODELocalAssembler.h"  // TO BE CHANGED
//#include "NonLinearGIAJacabianLocalAssembler.h" // TO BE CHANGED
#include "ReductConc.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"

template <class FunctionReductConc>
class NestedGIALocalProbNRIterationStepInitializer
{
public:
    NestedGIALocalProbNRIterationStepInitializer(FunctionReductConc & ReductConc)
      : _ReductConc(ReductConc) {};

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
    	_ReductConc.update_xi_global_nodal_values( x_new );

    	// call the local problem
    	_ReductConc.calc_nodal_local_problem(_ReductConc.getTimeStep()->getTimeStepSize(), 1.0E-8, 1.0E-12, 50);
    };

private:
    /**
      * pointer to assemblerJ
      */ 
    //NonLinearGIAJacobianLocalAssembler<T_NODAL_FUNCTION_SCALAR, T_FUNCTION_DATA>* _assemblerJ;
    FunctionReductConc & _ReductConc;
};

#endif  // end of ifndef
