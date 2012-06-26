/**************************************************************************
   FEMLib-Object: MAT-MP
   Task: MediumProperties
   Programing:
   01/2004 OK Implementation
**************************************************************************/

#ifndef rfmat_mp_new_INC
#define rfmat_mp_new_INC
/* Schutz gegen mehrfaches Einfuegen */

// C++ STL
#include <list>
#include <string>
#include <vector>
//#include <fstream>

// GeoLib
//#include "GeoType.h"
#include "makros.h" // JT

// PCSLib
//#include "rf_pcs.h"

class CMediumProperties
{
public:
private:
	//WW
	// Data base
	void SetMediumPropertiesDefaultsClay(void); //CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsSilt(void); //CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsSand(void); //CMCD 9/2004 GeoSys 4
	                                            //CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsGravel(void);
	//CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsCrystalline(void);
	//CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsBordenAquifer(void);
	// Porosity
	//CMCD 9/2004 GeoSys 4
	double PorosityEffectiveStress(long,double );
	double PorosityVolumetricFreeSwellingConstantIonicstrength(long,double,double);
	//MX 1/2005
	double PorosityEffectiveConstrainedSwelling(long,double,double, double*);
	//MX 1/2005
	double PorosityVolumetricFreeSwelling(long,double,double);
	//MX 1/2005
	double PorosityEffectiveConstrainedSwellingConstantIonicStrength(long,
	                                                                 double,
	                                                                 double,
	                                                                 double*);
	// Permeability
	// Permeabilty stress corrector WW
	double* c_coefficient;
	unsigned geo_dimension;
	int permeability_stress_mode;
	//
	//CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod1(long,double );
	//CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod2(long,double );
	//CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod3(long,double );
	//CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod4(long,double, double );
	friend class CMediumPropertiesGroup;
public:
	// Methods
	CMediumProperties(void);              // constructor
	~CMediumProperties(void);             // destructor
	CMediumProperties* Get(std::string);
	CMediumProperties* GetDB(std::string);
	CMediumProperties* GetByGroupNumber(int);
	void Set(std::string,std::string,double);
	void SetDB(std::string,std::string,double);
	int GetPropertyType(std::string);
	std::ios::pos_type Read(std::ifstream*);
	void Write(std::fstream*);

	// Properties
private:
	// PCS
	std::string pcs_type_name;            //YD
public:
	std::vector<std::string>pcs_name_vector;
private:
	std::vector<std::string> porosity_pcs_name_vector;

public:
	//GEO
	std::string geo_name;
	std::vector<std::string>geo_name_vector; //OK
	double geo_area;
	std::string geo_area_file;            //OK

	double density;
	std::string name;
	int number;
	int porosity_model;                   // porosity
	int porosity_curve;
	double porosity_model_values[15];
	double porosity;
	double KC_porosity_initial;           // HS 11.2008
	double KC_permeability_initial;       // HS 11.2008
	std::string porosity_file;            //OK/MB
	int tortuosity_model;
	double tortuosity_model_values[10];
	double tortuosity;
	int flowlinearity_model;
	double flowlinearity_model_values[10];
	int storage_model;                    // storativity
	double storage_model_values[10];
	double storage;
	int conductivity_model;
	double conductivity;
	int unconfined_flow_group;
	int permeability_model;               // permeability
	double permeability;
	double permeability_tensor[9];       
	std::string permeability_tensor_type_name;
	std::string tortuosity_tensor_type_name;
	int permeability_tensor_type;
	int tortuosity_tensor_type;

	int permeability_pressure_model;
	double permeability_pressure_model_values[10];
	double permeability_pressure_rel;
	int permeability_strain_model;        //WX: permeability function strain model. 05.2010
	int permeability_strain_model_value[3]; //WX:permeability fuction strain model value. 05.2010
	//
	// Relative permeability (JT)
	int permeability_saturation_model[MAX_FLUID_PHASES];
	double minimum_relative_permeability;
	double residual_saturation[MAX_FLUID_PHASES];
	double maximum_saturation[MAX_FLUID_PHASES];
	double saturation_exponent[MAX_FLUID_PHASES];
	double perm_saturation_value[MAX_FLUID_PHASES];
	//
	std::string permeability_file;        //SB //OK/MB string permeability_dis_type_file;
	std::string tortuosity_file;          // PCH
	bool entry_pressure_conversion;		//JT
	int capillary_pressure_model;
	int permeability_porosity_model;
	double permeability_porosity_model_values[10];
	double storativity;
	double capillary_pressure_values[5];		//JT2012
	double heat_capacity;                 // thermal properties
	int mass_dispersion_model;
	double mass_dispersion_longitudinal;
	double mass_dispersion_transverse;
	double lgpn;                          //local grid peclet number
	int heat_dispersion_model;
	double heat_dispersion_longitudinal;
	double heat_dispersion_transverse;
	double heat_conductivity_tensor[9];
	int fct_number;                       // functions
	int heat_diffusion_model;
	int evaporation;                      // if it is 647 then evaporation ON, else OFF: and corresponding heat loss will compensated by heat ST
	double heatflux;
	double vaporfraction;
	//aux
	int m_color[3];
	bool selected;
	int mode;
	// surface water
	// JOD
	double friction_coefficient, friction_exp_slope, friction_exp_depth;
	// JOD
	double overland_width, rill_height, rill_epsilon;
	bool channel;
	double argument;                      //OK
	//Dual Richards transfer coefficient  YD
	double transfer_coefficient;
	//double unsaturated_hydraulic_conductivity;
	double specific_storage;
	int vol_mat_model;                    // CB
	double vol_mat;                       //SB
	int vol_bio_model;                    // CB
	double vol_bio;                       //SB
	double foc;                           // organic carbon content
};

class CMediumPropertiesGroup                      //YD
{
public:
	CMediumPropertiesGroup() {OrigSize = 0; }
	std::string pcs_name;
	std::string pcs_type_name;
	std::vector<CMediumProperties*>mmp_group_vector;
private:
	int OrigSize;                         // For excavation simulation.
};

//YD
extern std::list<CMediumPropertiesGroup*>mmp_group_list;
extern void MMPGroupDelete(/*string pcs_type_name*/);

extern std::vector<CMediumProperties*>mmp_vector;
extern void MATLoadDB(std::string);
extern std::list<std::string> keywd_list;         //keyword-referenzliste "kw"
extern void comp_keywd_list(std::string);
extern void read_keywd_list(void);
extern std::list<std::string>mat_name_list;

extern void MMPWrite(std::string);
extern bool MMPRead(std::string);
extern void MMPWriteTecplot(std::string);
extern void MMPDelete();
extern CMediumProperties* MMPGet(const std::string&);
extern bool MMPExist(std::ifstream* mmp_file);    //OK
extern bool MMPExist();                           //OK

#define MMP_FILE_EXTENSION ".mmp"
#endif
