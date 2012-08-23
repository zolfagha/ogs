/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReductionKin.h
 *
 * Created on 2012-08-23 by Haibing Shao
 */
#ifndef CHEMREDUCTIONKIN_H
#define CHEMREDUCTIONKIN_H

#include "chemconst.h"

namespace ogsChem
{

class chemReductionKin
{
public:
	
	/**
      * constructor of the class
      */
	chemReductionKin(void);
	
	/**
      * destructor of the class
      */
	~chemReductionKin(void); 

    /**
      * convert concentration vector to eta and xi vector
      */
	void Conc2EtaXi();

	/**
      * convert eta and xi vector to concentration vector
      */
	void EtaXi2Conc(); 

private:
    /**
      * stoichiometric matrix S
      */	
	LocalMatrix _matStoi; 

	/**
      * sub-matrixes of S
      */	
	LocalMatrix _matS_1, _matS_2, _matS_1_bar, _matS_2_bar; 

	/**
      * complementary orthorgonal matrixes
      */
	LocalMatrix _matS_1_ast, _matS_2_ast, _matS_1_bar_ast, _matS_2_bar_ast; 

};

}

#endif