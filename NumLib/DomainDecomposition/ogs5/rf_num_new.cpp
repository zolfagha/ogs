/**************************************************************************
   FEMLib - Object: NUM
   Task:
   Programing:
   11/2004 OK Implementation
   last modified:
**************************************************************************/
// There is a name conflict between stdio.h and the MPI C++ binding
// with respect to the names SEEK_SET, SEEK_CUR, and SEEK_END.  MPI
// wants these in the MPI namespace, but stdio.h will #define these
// to integer values.  #undef'ing these can cause obscure problems
// with other include files (such as iostream), so we instead use
// #error to indicate a fatal error.  Users can either #undef
// the names before including mpi.h or include mpi.h *before* stdio.h
// or iostream.
#ifdef USE_MPI                                    //WW
#include "mpi.h"
#include "par_ddc.h"
//#undef SEEK_SET
//#undef SEEK_END
//#undef SEEK_CUR
#endif

//#include "makros.h"
// C++ STL
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
using namespace std;
// FEM-Makros
//#include "files0.h"
//#include "makros.h"
extern ios::pos_type GetNextSubKeyword(ifstream* file,string* line, bool* keyword);
// GeoSys-GeoLib
// GeoSys-FEMLib
#include "rf_num_new.h"
//#include "StringTools.h"
//#include "mathlib.h"

namespace NumLib
{
namespace OGS5
{

/**************************************************************************
   FEMLib-Method:
   Task: constructor
   Programing:
   11/2004 OK Implementation
   10/2005 OK pcs_type_name
   07/2007 OK DEFORMATION
**************************************************************************/
CNumerics::CNumerics(string name)
{
	pcs_type_name = name;                 //OK
	// GLOBAL
	renumber_method = 0;
	// LS - Linear Solver
	ls_method = 2;                        //OK41
	ls_max_iterations = 1000;
	ls_error_method = 1;
	ls_error_tolerance = 1e-12;
	ls_theta = 1.0;
	ls_precond = 1;
	ls_storage_method = 2;                //OK41
	m_cols = 5;                           // 06.2010. WW
	ls_extra_arg = ""; //NW
	// NLS - Nonlinear Solver
	nls_method_name = "PICARD";
	nls_method = -1;                       // Default linear, 0: Picard. 1: Newton
	nls_max_iterations = 1;               //OK
	nls_error_tolerance = 1.0e-4;
	nls_error_tolerance_local = 1.0e-10;  //For element level
	nls_relaxation = 0.0;
	// cpl WW
	cpl_iterations = 1;                   //OK
	cpl_tolerance = 1.0e-3;
	cpl_variable = "FLUX";
	// ELE
	ele_gauss_points = 3;
	ele_mass_lumping = 0;
	ele_upwind_method = 0;                //CB
	ele_upwinding = 0;
	ele_supg_method = 0;                  //NW
	ele_supg_method_length = 0;           //NW
	ele_supg_method_diffusivity = 0;      //NW
	fct_method = -1;                      //NW
	fct_prelimiter_type = 0;              //NW
	fct_const_alpha = -1.0;               //NW
	//----------------------------------------------------------------------
	// Deformation
	GravityProfile = 0;
	DynamicDamping = NULL;                //WW
	if(pcs_type_name.compare("DEFORMATION") == 0)
	{
		ls_method = 2;
		ls_error_method = 2;
		ls_error_tolerance = 1e-12;
		ls_max_iterations = 2000;
		ls_precond = 100;
		ls_storage_method = 4;
	}
	//----------------------------------------------------------------------
	if(pcs_type_name.compare("RICHARDS_FLOW") == 0)
	{
		ele_mass_lumping = 1;
		ele_upwinding = 0.5;
		ls_max_iterations = 2000;
		ls_error_method = 2;
		ls_error_tolerance = 1e-10;
		ls_precond = 4;
		ls_storage_method = 4;
		nls_max_iterations = 25;
		nls_error_tolerance = 1.0e-3;
	}
	//
}

/**************************************************************************
   FEMLib-Method:
   Task: deconstructor
   Programing:
   11/2004 OK Implementation
**************************************************************************/
CNumerics::~CNumerics(void)
{
	if(DynamicDamping)
		delete [] DynamicDamping;
	DynamicDamping = NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 WW Implementation
**************************************************************************/
bool CNumerics::CheckDynamic()
{
	if(DynamicDamping)
		return true;
	else
		return false;
}

#if 0
/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   11/2004 OK Implementation
**************************************************************************/
ios::pos_type CNumerics::Read(ifstream* num_file)
{
	string line_string;
	bool new_keyword = false;
	bool new_subkeyword = false;
	ios::pos_type position;
	ios::pos_type position_subkeyword;
	std::stringstream line;
	//========================================================================
	// Schleife ueber alle Phasen bzw. Komponenten
	while(!new_keyword)
	{
		if(new_subkeyword)
			num_file->seekg(position,ios::beg);
		new_subkeyword = false;
		position = GetNextSubKeyword(num_file,&line_string,&new_keyword);
		if(new_keyword)
			return position;
		//....................................................................
		// subkeyword found
		if(line_string.find("$PCS_TYPE") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> pcs_type_name;
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if(line_string.find("$RENUMBER") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> renumber_method;
			if(renumber_method == 2)
				line >> renumber_parameter;
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if(line_string.find("$NON_LINEAR_SOLVER") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> nls_method_name;
			line >> nls_error_tolerance;
			nls_method = 0;
			if(nls_method_name.find("NEWTON") != string::npos)
			{
				nls_method = 1;
				line >> nls_error_tolerance_local;
			}
			// 07.2010. WW
			///  Jacobian free Newton-Krylov method
			if(nls_method_name.find("JFNK") != string::npos)
			{
				nls_method = 2;
				line >> nls_error_tolerance_local;
			}
			line >> nls_max_iterations;
			line >> nls_relaxation;
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if(line_string.find("$LINEAR_SOLVER") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> ls_method;
			line >> ls_error_method >> ls_error_tolerance;
			line >> ls_max_iterations;
			line >> ls_theta;
			line >> ls_precond;
			line >> ls_storage_method;
			/// For GMRES. 06.2010. WW
			if(ls_method == 13)
				line >> m_cols;
			line.clear();
			continue;
		}
		//....................................................................
		if(line_string.find("$EXTERNAL_SOLVER_OPTION") != string::npos) // subkeyword found
		{
			ls_extra_arg = GetLineFromFile1(num_file);
			trim(ls_extra_arg);
			continue;
		}
		//....................................................................
		// subkeyword found
		if(line_string.find("$ELE_GAUSS_POINTS") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> ele_gauss_points; // probably element-type-wise
			line.clear();
			continue;
		}
		// subkeyword found
		if(line_string.find("$ELE_MASS_LUMPING") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> ele_mass_lumping;
			line.clear();
			continue;
		}
		// subkeyword found
		if(line_string.find("$ELE_UPWINDING") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			//CB now read also upwinding method
			line >> ele_upwinding >> ele_upwind_method;
			line.clear();
			continue;
		}
		// subkeyword found
		if(line_string.find("$ELE_SUPG") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			//NW
			line >> ele_supg_method >> ele_supg_method_length >>
			ele_supg_method_diffusivity;
			line.clear();
			cout << "->SUPG method is selected." << endl;
			continue;
		}
		// subkeyword found
		if(line_string.find("$GRAVITY_PROFILE") != string::npos)
		{
			line.str(GetLineFromFile1(num_file)); //WW
			line >> GravityProfile;
			line.clear();
			continue;
		}
		// subkeyword found
		if(line_string.find("$DYNAMIC_DAMPING") != string::npos)
		{
			line.str(GetLineFromFile1(num_file)); //WW
			DynamicDamping = new double[3];
			// Default
			DynamicDamping[0] = 0.515;
			DynamicDamping[1] = 0.51;
			DynamicDamping[2] = 0.51;
			line >> DynamicDamping[0] >> DynamicDamping[1] >> DynamicDamping[2];
			line.clear();
			continue;
		}
		// WW subkeyword found
		if(line_string.find("$COUPLING_ITERATIONS") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> cpl_variable; //pcs_name. WW MB
			line >> cpl_iterations;
			line >> cpl_tolerance;
			line.clear();
			continue;
		}
		//Flux corrected transport by Kuzmin (2009)
		// NW
		if(line_string.find("$FEM_FCT") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> fct_method; //1: linearized FCT
			line >> fct_prelimiter_type; //0: just cancel, 1: minmod, 2: superbee
			line >> fct_const_alpha; //-1: off, [0.0,1.0] 0: Upwind, 1: Galerkin
			line.clear();
			cout << "->FEM_FCT method is selected." << endl;
			continue;
		}

		//....................................................................
		/*
		    if(line_string.find("$TIME_STEPS")!=string::npos) { // subkeyword found
		      while((!new_keyword)||(!new_subkeyword)||(!num_file->eof())){
		        position = num_file->tellg();
		        line_string = GetLineFromFile1(num_file);
		        if(line_string.find("#")!=string::npos){
		          return position;
		        }
		        if(line_string.find("$")!=string::npos){
		          new_subkeyword = true;
		          break;
		   }
		   line.str(line_string);
		   line >> no_time_steps;
		   line >> time_step_length;
		   for(i=0;i<no_time_steps;i++)
		   time_step_vector.push_back(time_step_length);
		   line.clear();
		   }
		   }
		 */
		//....................................................................
	}
	return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: master write function
   Programing:
   11/2004 OK Implementation
   last modification:
**************************************************************************/
void NUMWrite(string base_file_name)
{
	CNumerics* m_num = NULL;
	string sub_line;
	string line_string;
	//========================================================================
	// File handling
	string num_file_name = base_file_name + NUM_FILE_EXTENSION;
	fstream num_file (num_file_name.data(),ios::trunc | ios::out);
	num_file.setf(ios::scientific,ios::floatfield);
	num_file.precision(12);
	if (!num_file.good())
		return;
	num_file.seekg(0L,ios::beg);
	//========================================================================
	num_file << "GeoSys-NUM: Numerics ------------------------------------------------\n";
	//========================================================================
	// OUT vector
	int num_vector_size = (int)num_vector.size();
	int i;
	for(i = 0; i < num_vector_size; i++)
	{
		m_num = num_vector[i];
		m_num->Write(&num_file);
	}
	num_file << "#STOP";
	num_file.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: write function
   Programing:
   11/2004 OK Implementation
   last modification:
**************************************************************************/
void CNumerics::Write(fstream* num_file)
{
	//KEYWORD
	*num_file << "#NUMERICS" << endl;
	//--------------------------------------------------------------------
	/*OK
	   *num_file << " $METHOD" << endl;
	   *num_file << method_name << endl;
	   if(method_name.find("LAGRANGE")!=string::npos){
	   *num_file << lag_quality << " " << lag_max_steps << " " << lag_local_eps << " ";
	   *num_file << lag_time_weighting << " " << lag_min_weight << " ";
	   *num_file << lag_use_matrix << " " << lag_vel_method;
	   *num_file << endl;
	   }
	 */
	//--------------------------------------------------------------------
	*num_file << " $PCS_TYPE" << endl;
	*num_file << "  " << pcs_type_name << endl;
	//--------------------------------------------------------------------
	*num_file << " $NON_LINEAR_SOLVER" << endl;
	*num_file << "  " << nls_method_name;
	*num_file << " "  << nls_error_tolerance;
	*num_file << " "  << nls_max_iterations;
	*num_file << " "  << nls_relaxation;
	*num_file << endl;
	//--------------------------------------------------------------------
	*num_file << " $LINEAR_SOLVER" << endl;
	*num_file << "  " << ls_method;
	*num_file << " "  << ls_error_method;
	*num_file << " "  << ls_error_tolerance;
	*num_file << " "  << ls_max_iterations;
	*num_file << " "  << ls_theta;
	*num_file << " "  << ls_precond;
	*num_file << " "  << ls_storage_method;
	*num_file << endl;
	//--------------------------------------------------------------------
	*num_file << " $ELE_GAUSS_POINTS" << endl;
	*num_file << "  " << ele_gauss_points;
	*num_file << endl;
	//--------------------------------------------------------------------
	*num_file << " $ELE_MASS_LUMPING" << endl;
	*num_file << "  " << ele_mass_lumping;
	*num_file << endl;
	//--------------------------------------------------------------------
	*num_file << " $ELE_UPWINDING" << endl;
	*num_file << "  " << ele_upwinding;
	*num_file << endl;
	//--------------------------------------------------------------------
}

#endif

}
}