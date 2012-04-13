
#pragma once

#include <valarray>
#include <vector>
#include <algorithm>

#include "MathLib/Function/IFunction.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "SystemOfEquations.h"


namespace MathLib
{

class MyFunction : public TemplateFunction<std::map<size_t, std::valarray<double> >, DenseLinearEquations >
{
public:
	MyFunction(const std::vector<Variable*> &variables, const std::vector<LinearEquation*> &equations) : _variables(variables), _equations(equations) {};
	virtual ~MyFunction() {};

    virtual void eval(const std::map<size_t, std::valarray<double> > &inactive_x, DenseLinearEquations &eqs)
    {
    	size_t n_dof = 0;
    	for (size_t i=0; i<_variables.size(); i++) n_dof += _variables[i]->n_dof;

    	eqs.create(n_dof);

//    	Matrix<double> A(n_dof, n_dof);
//    	std::valarray<double> b(n_dof);
    	Matrix<double>* A = eqs.getA();
    	double* b = eqs.getRHS();

		size_t row_offset = 0;
    	for (size_t i=0; i<_equations.size(); i++) {
    		LinearEquation* eq = _equations[i];
			const size_t n_local_row_dof = eq->n_dof;
			size_t col_offset = 0;
    		for (size_t j=0; j<eq->variables.size(); j++) {
    			Variable* var = eq->variables[j];
    			LinearComponent* comp = eq->componets[j];
    			const size_t n_local_col_dof = var->n_dof;
                std::vector<size_t> pos_row(n_local_row_dof);
                for (size_t k=0; k<n_local_row_dof; k++) pos_row[k] = row_offset + k;

                if (isActiveVariable(var)) {
                    std::vector<size_t> pos_col(n_local_col_dof);
                    for (size_t k=0; k<n_local_col_dof; k++) pos_col[k] = col_offset + k;
                    if (comp->stiffness!=0) {
                        A->add(pos_row, pos_col, *comp->stiffness);
                    }
                    col_offset += n_local_col_dof;
                } else {
                    std::map<size_t, std::valarray<double> >::const_iterator itr = inactive_x.find(var->id);
                    if (itr!=inactive_x.end()) {
                        const std::valarray<double> &x = itr->second;
                        if (comp->stiffness!=0) {
                            std::valarray<double> vec(n_local_row_dof);
                            comp->stiffness->axpy(1.0, &x[0], .0, &vec[0]);
                            for (size_t k=0; k<n_local_row_dof; k++)
                                b[pos_row[k]] -= vec[k];
                        }
                    }
                }

    			if (comp->source!=0) {
					for (size_t k=0; k<n_local_row_dof; k++)
						b[pos_row[k]] += (*comp->source)[k];
    			}
    		}
    		row_offset += n_local_row_dof;
    	}

    }

    virtual MyFunction* clone() const
	{
    	return new MyFunction(_variables, _equations);
	}

private:
    bool isActiveVariable(Variable* var)
    {
        for (size_t i=0; i<_variables.size(); i++)
            if (_variables[i]->id==var->id) return true;
        return false;
    }


private:
	const std::vector<Variable*> _variables;
	const std::vector<LinearEquation*> _equations;
};

} //end
