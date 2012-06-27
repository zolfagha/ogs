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

#include "rf_num_new.h"

#include "makros.h"
// C++ STL
#include <cfloat>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <list>
#include <string>
#include "readNonBlankLineFromInputStream.h"
using namespace std;
//// FEM-Makros
//#include "files0.h"
//#include "makros.h"
//extern ios::pos_type GetNextSubKeyword(ifstream* file,string* line, bool* keyword);
//// GeoSys-GeoLib
//// GeoSys-FEMLib
//#ifndef NEW_EQS                                   //WW. 06.11.2008
//#include "matrix.h"
//#endif
//#include "StringTools.h"
//#include "mathlib.h"
//#include "rf_pcs.h"
//#include "tools.h"
//// GeoSys-MSHLib

extern size_t max_dim;                            //OK411 todo

//==========================================================================
vector<CNumerics*>num_vector;

std::ios::pos_type GetNextSubKeyword(std::ifstream* file,std::string* line, bool* keyword)
{
    char buffer[MAX_ZEILE];
    std::ios::pos_type position;
    position = file->tellg();
    *keyword = false;
    std::string line_complete;
    int i,j;
    // Look for next subkeyword
    while(!(line_complete.find("$") != std::string::npos) && (!file->eof()))
    {
        file->getline(buffer,MAX_ZEILE);
        line_complete = buffer;
        if(line_complete.find("#") != std::string::npos)
        {
            *keyword = true;
            return position;
        }
        //Anf�ngliche Leerzeichen �berlesen, i=Position des ersten Nichtleerzeichens im string
        i = (int) line_complete.find_first_not_of(" ",0);
        j = (int) line_complete.find(";",i); //Nach Kommentarzeichen ; suchen. j = Position des Kommentarzeichens, j=-1 wenn es keines gibt.
        if(j < 0)
            j = (int)line_complete.length();
        //if(j!=i) break;						 //Wenn das erste nicht-leerzeichen ein Kommentarzeichen ist, zeile �berlesen. Sonst ist das eine Datenzeile
        if(i != -1)
            *line = line_complete.substr(i,j - i);  //Ab erstem nicht-Leerzeichen bis Kommentarzeichen rauskopieren in neuen substring, falls Zeile nicht leer ist
    }
    return position;
}

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
	//
	// LS - Linear Solver
	ls_method = 2;                        //OK41
	ls_max_iterations = 1000;
	ls_error_method = 1;
	ls_error_tolerance = 1e-12;
	ls_theta = 1.0;
	ls_precond = 1;
	ls_storage_method = 2;					//OK41
	m_cols = 5;								// 06.2010. WW
	ls_extra_arg = ""; //NW
	//
	// NLS - Nonlinear Solver
	nls_method_name = "PICARD";
	nls_method = -1;						//Default linear, 0: Picard. 1: Newton. 2:JFNK
	nls_error_method = 1;					//JT2012
	nls_max_iterations = 1;					//OK
	nls_relaxation = 0.0;
	for(size_t i=0; i<DOF_NUMBER_MAX; i++)	//JT2012
		nls_error_tolerance[i] = -1.0;		//JT2012: should not default this. Should always be entered by user!
	//
	// CPL - Coupled processes
	cpl_error_specified = false;
	cpl_master_process = false;
	cpl_process = "INVALID_PROCESS";		//JT2012: do not couple with any process, unless indicated
	cpl_variable = "NONE";
	cpl_variable_JOD = "FLUX";
	cpl_max_iterations = 1;					//OK
	cpl_min_iterations = 1;					//JT2012
	for(size_t i=0; i<DOF_NUMBER_MAX; i++)	//JT2012
		cpl_error_tolerance[i] = -1.0;			//JT2012: should not default this. Should always be entered by user!
	//
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
   Task: OBJ read function
   Programing:
   11/2004 OK Implementation
**************************************************************************/
bool NUMRead(string file_base_name)
{
	//----------------------------------------------------------------------
	//OK  NUMDelete();
	//----------------------------------------------------------------------
	CNumerics* m_num = NULL;
	char line[MAX_ZEILE];
	bool overall_coupling_exists=false; //JT
	string sub_line;
	string line_string;
	ios::pos_type position;
	//========================================================================
	// File handling
	string num_file_name = file_base_name + NUM_FILE_EXTENSION;
	ifstream num_file (num_file_name.data(),ios::in);
	if (!num_file.good())
		return false;
	num_file.seekg(0L,ios::beg);
	//========================================================================
	// Keyword loop
	cout << "NUMRead" << endl;
	while (!num_file.eof())
	{
		num_file.getline(line,MAX_ZEILE);
		line_string = line;
		if(line_string.find("#STOP") != string::npos)
			return true;
		//
		if(line_string.find("$OVERALL_COUPLING") != string::npos){
			overall_coupling_exists = true; // JT: for error checking
		}
		//----------------------------------------------------------------------
		// keyword found
		if(line_string.find("#NUMERICS") != string::npos)
		{
			m_num = new CNumerics("default");
			position = m_num->Read(&num_file);
			num_vector.push_back(m_num);
			num_file.seekg(position,ios::beg);
			m_num->NumConfigure(overall_coupling_exists);					  // JT2012
		}                         // keyword found
	}                                     // eof
	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   11/2004 OK Implementation
**************************************************************************/
ios::pos_type CNumerics::Read(ifstream* num_file)
{
	string line_string;
	std::string error_method_name;
	std::string coupling_target;
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
			line.str(readNonBlankLineFromInputStream(*num_file));
			line >> pcs_type_name;
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if(line_string.find("$RENUMBER") != string::npos)
		{
			line.str(readNonBlankLineFromInputStream(*num_file));
			line >> renumber_method;
			if(renumber_method == 2)
				line >> renumber_parameter;
			line.clear();
			continue;
		}
		//....................................................................
		// JT->WW: Local tolerance previously found in $NON_LINEAR_SOLVER for NEWTON. Moved here for now.
		if(line_string.find("$PLASTICITY_TOLERANCE") != string::npos)
		{
			line.str(readNonBlankLineFromInputStream(*num_file));
			line >> nls_plasticity_local_tolerance;
		}
		//....................................................................
		// subkeyword found ($NON_LINEAR_ITERATION  -or-  $NON_LINEAR_ITERATIONS)
		if(line_string.find("$NON_LINEAR_ITERATION") != string::npos)
		{
			// JT:	in >> nls_method_name
			//		in >> error_method_name
			//		in >> max iter
			//		in >> relaxation
			//		in >> tolerance[1:dof]
			//
			line.str(readNonBlankLineFromInputStream(*num_file));
			line >> nls_method_name;
			line >> error_method_name;
			line >> nls_max_iterations;
			line >> nls_relaxation;
			//
			//setNonLinearErrorMethod(FiniteElement::convertErrorMethod(error_method_name));
			//switch(getNonLinearErrorMethod())
			//{
			//	case FiniteElement::ENORM: // only 1 tolerance required
			//		line >> nls_error_tolerance[0];
			//		break;
			//	//
			//	case FiniteElement::ERNORM: // only 1 tolerance required
			//		line >> nls_error_tolerance[0];
			//		break;
			//	//
			//	case FiniteElement::EVNORM: // 1 tolerance for each primary variable (for Deformation, only 1 tolerance required. Applies to x,y,z)
			//		for(int i=0; i<DOF_NUMBER_MAX; i++)
			//			line >> nls_error_tolerance[i];
			//		break;
			//	//
			//	case FiniteElement::LMAX: // 1 tolerance for each primary variable (for Deformation, only 1 tolerance required. Applies to x,y,z)
			//		for(int i=0; i<DOF_NUMBER_MAX; i++)
			//			line >> nls_error_tolerance[i];
			//		break;
			//	//
			//	case FiniteElement::BNORM: // only 1 tolerance required
			//		line >> nls_error_tolerance[0];
			//		break;
			//	//
			//	default:
			//		ScreenMessage("ERROR in NUMRead. Invalid non-linear iteration error method selected.\n");
			//		exit(1);
			//		break;
			//}

			nls_method = 0;
			if(nls_method_name.find("NEWTON") != string::npos)
				nls_method = 1;
			else if(nls_method_name.find("JFNK") != string::npos) //  Jacobian free Newton-Krylov method
				nls_method = 2;
			//
			line.clear();
			continue;
		}
		else if(line_string.find("$NON_LINEAR_SOLVER") != string::npos)
		{
			//ScreenMessage(" --\n Using old $NON_LINEAR_SOLVER keyword.\n");
			//ScreenMessage(" Eventually this will be obsolete. Consider switching to\n");
			//ScreenMessage(" $NON_LINEAR_ITERATIONS for better results and greater flexibility.\n");
			//ScreenMessage(" --\n");
			////
			//// JT:	in >> method_name
			////		in >> tolerance
			////		if(NEWTON) in >> tolerance_local
			////		in >> max iter
			////		in >> relaxation
			////
			////
			//line.str(readNonBlankLineFromInputStream(*num_file));
			//line >> nls_method_name;
			////
			//nls_method = 0;
			//if(nls_method_name.find("NEWTON") != string::npos)
			//	nls_method = 1;
			//else if(nls_method_name.find("JFNK") != string::npos) //  Jacobian free Newton-Krylov method
			//	nls_method = 2;
			////
			//if(nls_method > 0){
			//	line >> nls_error_tolerance[0];
			//	line >> nls_plasticity_local_tolerance;
			//	error_method_name = "BNORM"; // JT: this is hardwired in old version
			//}
			//else{
			//	line >> nls_error_tolerance[0];
			//	error_method_name = "LMAX"; // JT: this is hardwired in old version
			//}
			//setNonLinearErrorMethod(FiniteElement::convertErrorMethod(error_method_name));
			//
			line >> nls_max_iterations;
			line >> nls_relaxation;
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if(line_string.find("$LINEAR_SOLVER") != string::npos)
		{
			line.str(readNonBlankLineFromInputStream(*num_file));
			line >> ls_method;
			line >> ls_error_method;
			line >> ls_error_tolerance;
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
		// JT subkeyword found
		if(line_string.find("$COUPLING_ITERATIONS") != string::npos)
		{
			//ScreenMessage("$COUPLING_ITERATIONS keyword obsolete.\n");
			//ScreenMessage("Use $COUPLING_CONTROL and $COUPLED_PROCESS for process couplings.\n");
			exit(1);
		}
		//....................................................................
		// JT subkeyword found
		if(line_string.find("$COUPLING_CONTROL") != string::npos) // JT: For this process, how is coupling error handled?
		{
			// JT:	in >> error_method_name
			//		in >> tolerance[1:dof]
			//
			line.str(readNonBlankLineFromInputStream(*num_file));
			line >> error_method_name;
			//
			cpl_error_specified = true;
			//setCouplingErrorMethod(FiniteElement::convertErrorMethod(error_method_name));
			//switch(getCouplingErrorMethod())
			//{
			//	case FiniteElement::ENORM: // only 1 tolerance required
			//		line >> cpl_error_tolerance[0];
			//		break;
			//	//
			//	case FiniteElement::ERNORM: // only 1 tolerance required
			//		line >> cpl_error_tolerance[0];
			//		break;
			//	//
			//	case FiniteElement::EVNORM: // 1 tolerance for each primary variable (for Deformation, only 1 tolerance required. Applies to x,y,z)
			//		for(int i=0; i<DOF_NUMBER_MAX; i++)
			//			line >> cpl_error_tolerance[i];
			//		break;
			//	//
			//	case FiniteElement::LMAX: // 1 tolerance for each primary variable (for Deformation, only 1 tolerance required. Applies to x,y,z)
			//		for(int i=0; i<DOF_NUMBER_MAX; i++)
			//			line >> cpl_error_tolerance[i];
			//		break;
			//	//
			//	case FiniteElement::BNORM:
			//		ScreenMessage("ERROR in NUMRead. BNORM not configured for process couplings.\n");
			//		ScreenMessage("We suggest ENORM as a valid companion for NEWTON couplings.\n");
			//		exit(1);
			//		break;
			//	//
			//	default:
			//		ScreenMessage("ERROR in NUMRead. Invalid coupling error method selected.\n");
			//		exit(1);
			//		break;
			//}
			//
			line.clear();
			continue;
		}
		//....................................................................
		// JT subkeyword found
		if(line_string.find("$COUPLED_PROCESS") != string::npos) // JT: Is this process coupled to another process in an inner loop?
		{
			// in >> process name >> min iter >> max iter
			//
			line.str(readNonBlankLineFromInputStream(*num_file));
			line >> coupling_target;			// name of coupled process -OR- process variable
			line >> cpl_min_iterations;
			line >> cpl_max_iterations;
			//
			cpl_master_process = true;
			////
			//// Is coupling through a process name or a primary variable?
			//if(FiniteElement::convertPrimaryVariable(coupling_target) != FiniteElement::INVALID_PV){ // Then a valid process VARIABLE is entered. Use this.
			//	cpl_variable = coupling_target;
			//}
			//else if(PCSGet(coupling_target)){ // Then a valid process is entered
			//	cpl_process  = coupling_target;
			//}
			//else{
			//	ScreenMessage("WARNING. $COUPLED_PROCESS keyword encountered, but a valid process OR primary variable was not found.\n");
			//	cpl_master_process = false;
			//}
			//
			line.clear();
			continue;
		}
		//....................................................................
		if(line_string.find("$EXTERNAL_SOLVER_OPTION") != string::npos) // subkeyword found
		{
			ls_extra_arg = readNonBlankLineFromInputStream(*num_file);
			continue;
		}
		//....................................................................
		// subkeyword found
		if(line_string.find("$ELE_GAUSS_POINTS") != string::npos)
		{
			line.str(readNonBlankLineFromInputStream(*num_file));
			line >> ele_gauss_points; // probably element-type-wise
			line.clear();
			continue;
		}
		// subkeyword found
		if(line_string.find("$ELE_MASS_LUMPING") != string::npos)
		{
			line.str(readNonBlankLineFromInputStream(*num_file));
			line >> ele_mass_lumping;
			line.clear();
			continue;
		}
		// subkeyword found
		if(line_string.find("$ELE_UPWINDING") != string::npos)
		{
			line.str(readNonBlankLineFromInputStream(*num_file));
			//CB now read also upwinding method
			line >> ele_upwinding >> ele_upwind_method;
			line.clear();
			continue;
		}
		// subkeyword found
		if(line_string.find("$ELE_SUPG") != string::npos)
		{
			line.str(readNonBlankLineFromInputStream(*num_file));
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
			line.str(readNonBlankLineFromInputStream(*num_file)); //WW
			line >> GravityProfile;
			line.clear();
			continue;
		}
		// subkeyword found
		if(line_string.find("$DYNAMIC_DAMPING") != string::npos)
		{
			line.str(readNonBlankLineFromInputStream(*num_file)); //WW
			DynamicDamping = new double[3];
			// Default
			DynamicDamping[0] = 0.515;
			DynamicDamping[1] = 0.51;
			DynamicDamping[2] = 0.51;
			line >> DynamicDamping[0] >> DynamicDamping[1] >> DynamicDamping[2];
			line.clear();
			continue;
		}
		//Flux corrected transport by Kuzmin (2009)
		// NW
		if(line_string.find("$FEM_FCT") != string::npos)
		{
			line.str(readNonBlankLineFromInputStream(*num_file));
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
		        line_string = readNonBlankLineFromInputStream(*num_file);
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
