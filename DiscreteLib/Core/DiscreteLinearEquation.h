
#pragma once

#include <vector>

#include "DiscreteVector.h"

namespace DiscreteLib
{

class DofEquationIdTable;
class IDiscreteLinearEquationAssembler;

/** 
 * \brief Interface for discrete linear equations
 */
class IDiscreteLinearEquation
{
public:
	typedef DiscreteVector<double> GlobalVectorType;

	virtual ~IDiscreteLinearEquation() {};
    /// 
    virtual void initialize() = 0;
    /// construct 
    virtual void construct(IDiscreteLinearEquationAssembler& assemler) = 0;
    /// solve
    virtual void solve() = 0;
    /// get solution
    virtual void getGlobalX(std::vector<double> &x) = 0;
    virtual double* getLocalX() = 0;
    virtual void getX(GlobalVectorType &v) = 0;
//    virtual DiscreteVector<double>* getX() = 0;
    ///// get RHS 
    //virtual double* getRHS() = 0;
    ///// get a linear equation object
    //virtual MathLib::ILinearEquations* getLinearEquation() const = 0;
    /// get a Dof map manager
    virtual DofEquationIdTable* getDofMapManger() const = 0;
    /// set prescribed dof
    virtual void setPrescribedDoF(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_prescribed_values) = 0;
    /// set additional RHS values
    virtual void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_rhs_values, double fkt) = 0;
    /// set additional RHS values
    virtual void addRHS(const GlobalVectorType &v, double fkt) = 0;
};


}
