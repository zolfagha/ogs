/**************************************************************************
   FEMLib - Object: NUM
   Task: class implementation
   Programing:
   11/2004 OK Implementation
   last modified:
**************************************************************************/
#ifndef rf_num_new_INC
#define rf_num_new_INC

#include <fstream>
#include <string>
#include <vector>

namespace NumLib
{

class CNumerics
{
private:
	// cf. Computational Geomachanics pp.62 WW
	double* DynamicDamping;
	/// For GMRES solver. 30.06.2010. WW
	long m_cols;
public:
	// method
	std::string method_name;              //OK
	// PCS
	std::string pcs_type_name;
	// RENUMBER
	int renumber_method;
	int renumber_parameter;
	// LS - Linear Solver
	int ls_method;
	int ls_max_iterations;
	int ls_error_method;
	double ls_error_tolerance;
	double ls_theta;
	int ls_precond;
	int ls_storage_method;
	std::string ls_extra_arg; //NW
	// LS - Linear Solver
	std::string nls_method_name;
	int nls_method;                       // Picard or Newton
	int nls_error_method;                 //WW
	int nls_max_iterations;
	double nls_error_tolerance;
	double nls_error_tolerance_local;
	double nls_relaxation;
	// CPL WW
	double cpl_tolerance;
	int cpl_iterations;
	std::string cpl_variable;             // MB
	// ELE
	int ele_gauss_points;                 // probably element-type-wise
	int ele_mass_lumping;
	int ele_upwind_method;                //CB
	double ele_upwinding;
	int ele_supg_method;                  //NW
	int ele_supg_method_length;           //NW
	int ele_supg_method_diffusivity;      //NW
	//FEM-FCT
	int fct_method;                       //NW
	unsigned int fct_prelimiter_type;     //NW
	double fct_const_alpha;               //NW
	// Deformation
	int GravityProfile;
	// LAGRANGE method //OK
	double lag_quality;
	int lag_max_steps;
	double lag_local_eps;
	int lag_time_weighting;
	double lag_min_weight;
	int lag_use_matrix;
	int lag_vel_method;
	//
	// Dynamics
	bool CheckDynamic();
	double GetDynamicDamping_beta1 () const {return DynamicDamping[0]; }
	double GetDynamicDamping_beta2 () const {return DynamicDamping[1]; }
	double GetDynamicDamping_bbeta () const {return DynamicDamping[2]; }
	//
	/// For GMRES. WW
	long Get_m() const {return m_cols; }
	CNumerics(std::string);
	~CNumerics(void);
	//std::ios::pos_type Read(std::ifstream*);
	//void Write(std::fstream*);
};

}


#endif
