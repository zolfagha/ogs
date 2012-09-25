/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXKinReductionTransformFunction.h
 *
 * Created on 2012-09-25 by Haibing Shao
 */

#ifndef TX_KIN_REDUCTION_TRANSFORM_FUNCTION_H 
#define TX_KIN_REDUCTION_TRANSFORM_FUNCTION_H

#include "NumLib/Function/TXPosition.h"
#include "NumLib/Function/ITXFunction.h"
#include "ChemLib/chemReductionKin.h" 

namespace NumLib
{

/**
 *
 */
class TXKinReductionTransformFunction : public ITXFunction
{
public:
	typedef ogsChem::LocalVector LocalVector;

    TXKinReductionTransformFunction(std::vector<ITXFunction*> & vec_ITXFunc, 
		                            ogsChem::chemReductionKin* myKinReduction)
		: _vec_ITXFunc( vec_ITXFunc ), _myKinReduction( myKinReduction ) 
	{
		// initialize the memory
		local_conc       = LocalVector::Zero( _myKinReduction->get_n_Comp() ); 
		local_eta_mob    = LocalVector::Zero( _myKinReduction->get_n_eta_mob() ); 
		local_eta_immob  = LocalVector::Zero( _myKinReduction->get_n_eta() - _myKinReduction->get_n_eta_mob() ) ); 
		local_xi         = LocalVector::Zero( _myKinReduction->get_n_xi() );
	};

    virtual ~TXKinReductionTransformFunction() {};

    virtual void eval(const TXPosition &x, DataType &val) const
    {
		// for each concentrations, the value must first 
		// evaluated and stored in a vector
		size_t i; 
		for ( i=0; i < _myKinReduction->get_n_Comp() ; i++ )
		{			
			this->_vec_ITXFunc[i]->eval( x, local_conc[i] ); 
		}
		
		// convert concentrations to eta and xi
		this->_myKinReduction->Conc2EtaXi( local_conc, local_eta_mob, local_eta_immob, local_xi ); 
		// combine all the eta and xi vector to one output matrix
		val = LocalVector::Zero( local_eta_mob.size() + local_eta_immob.size() + local_xi.size() ); 
		val.head(local_eta_mob.size()) = local_eta_mob; 
		val.segment( local_eta_mob.size(), local_eta_immob.size() ) = local_eta_immob; 
		val.tail(local_xi.size()) = local_xi; 
	}

    virtual TXKinReductionTransformFunction* clone() const
    {
        return new TXKinReductionTransformFunction(_vec_ITXFunc, _myKinReduction);
    }

private:
    std::vector<ITXFunction*> _vec_ITXFunc;
    ogsChem::chemReductionKin* _myKinReduction;

	LocalVector local_conc; 
	LocalVector local_eta_mob; 
	LocalVector local_eta_immob; 
	LocalVector local_xi;  
};

} //end


#endif  // end of ifndef

