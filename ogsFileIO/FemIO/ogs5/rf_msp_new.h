/**************************************************************************
   FEMLib - Object: MAT-SP
   Task: class implementation
   Programing:
   08/2004 WW Implementation
   last modified:
**************************************************************************/
#ifndef rf_msp_new_INC
#define rf_msp_new_INC

// C++ STL
//#include <fstream>
#include <string>
#include <vector>
#include <iostream>

#define MSP_FILE_EXTENSION ".msp"

namespace Math_Group
{class Matrix;
}

    using Math_Group::Matrix;
/*---------------------------------------------------------------*/
class CSolidProperties
{
public:
	//
	CSolidProperties();
	~CSolidProperties();

	std::ios::pos_type Read(std::ifstream*);

public:
	// Material parameters
	double PoissonRatio;
	int Youngs_mode;
	int excavation;                       //12.2009. WW
	bool excavated;                       //12.2009. To be ..... WW
	Matrix* data_Youngs;
	double ThermalExpansion;
	//
	double s_tol;                         //16.06.2008 WW
	double f_tol;                         //16.06.2008 WW
	double biot_const;
	int bishop_model;                     //05.2011 WX
	double bishop_model_value;            //05.2011 WX
	double grav_const;                    //WW
	Matrix* data_Density;
	//
	Matrix* data_Capacity;
	Matrix* data_Conductivity;
	//
	Matrix* data_Plasticity;
	Matrix* data_Creep;
	//
	int Density_mode;
	//
	int Capacity_mode;
	int Conductivity_mode;
	int Plasticity_type;
	double primary_variable[10];          //CMCD
	double primary_variable_t0[10];       //CMCD
	double primary_variable_t1[10];       //CMCD
	// Creep property
	// 1. Stationary Norton model
	int Creep_mode;
	//
	bool axisymmetry;

	int mode;                             //CMCD
	// Swelling pressure
	int SwellingPressureType;
	double Max_SwellingPressure;
	//
	std::string CurveVariable_Conductivity;
	int CurveVariableType_Conductivity;
	// Secondary data
	// Elasticity
	double E;                             // Youngs moduls calculated from data_Youngs
	double Lambda;
	double G;                             // Shear stress modulus
	double K;                             // Bulk modulus

	// Rotation matrices and their transpose: UJG 25.11.2009
	Matrix* Crotm;                        // If this is needed by permaebility calculation, we keep it. Otherwise remove it. (To do, UJG/WW)
	Matrix* D_tran;

	// Plasticity
	double dl2;
	// 2. Single yield surface
	Matrix* d2G_dSdS;
	Matrix* d2G_dSdM;
	Matrix* LocalJacobi;                  // To store local Jacobi matrix
	Matrix* inv_Jac;                      // To store the inverse of the  Jacobi matrix
	Matrix* sumA_Matrix;
	double* rhs_l;                        // To store local unknowns of 15
	double* x_l;                          // To store local unknowns of 15
	int* Li;

	// Direct stress integration for Drucker-Prager
	double* devS;
	double* dFds;
	double* dGds;
	double* D_dFds;
	double* D_dGds;
	double* dFtds;                        //WX: 08.2010
	double* dGtds;                        //WX: 08.2010
	Matrix* ConstitutiveMatrix;           //WX: 08.2010
	// Thermal properties
	int thermal_conductivity_tensor_type;
	int thermal_conductivity_tensor_dim;
	double thermal_conductivity_tensor[9];
	std::string thermal_conductivity_tensor_type_name;

	std::string name;
	// IO
	std::string file_base_name;


	// Plasticity
	// 1. Drucker-Prager
	double Al;
	double Xi;
	double Y0;
	double BetaN;
	double Hard;
	double Hard_Loc;
	double tension;     //WX:08.2010 Tension strength

	//4. Mohr-Coulomb	//WX: 11.2010. Mohr-Coulomb model
	double Ntheta;
	double Nphi;
	double csn;

	//5. Hoek-Brown WX
	double HoekB_a;
	double HoekB_s;
	double HoekB_mb;
	double HoekB_sigci;
	double HoekB_tens;
	double HoekB_cohe;
	//WW
	std::vector<std::string>  capacity_pcs_name_vector;
	//WW
	std::vector<std::string>  conductivity_pcs_name_vector;
};

extern bool MSPRead(const std::string &file_base_name, std::vector<CSolidProperties*> &msp_vector);

#endif
