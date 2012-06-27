/**************************************************************************
   ROCKFLOW - Object: Process PCS
   Programing:
   02/2003 OK Implementation
   /2003 WW CRFProcessDeformation
   11/2003 OK re-organized
   07/2004 OK PCS2
   02/2005 WW/OK Element Assemblier and output
   12/2007 WW Classes of sparse matrix (jagged diagonal storage) and linear solver
           and parellelisation of them
   02/2008 PCH OpenMP parallelization for Lis matrix solver
**************************************************************************/
#include "rf_pcs.h"

//#include "FEMEnums.h"
//#include "Output.h"


#include "makros.h"
// C
#ifndef __APPLE__
#include <malloc.h>
#endif

// C++
#include <cfloat>
#include <iomanip>                                //WW
#include <iostream>
//#include <algorithm> // header of transform. WW
#include <set>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "readNonBlankLineFromInputStream.h"

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
//#include "files0.h"
//#include "par_ddc.h"
//#include "tools.h"
//#include "rf_pcs.h"
//#include "files0.h"

using namespace std;

// Globals, to be checked
int pcs_no_components = 0;
bool pcs_monolithic_flow = false;
int pcs_deformation = 0;
int dm_number_of_primary_nvals = 2;
bool show_onces_adp = true;
bool show_onces_mod = true;
bool show_onces_mod_flow = true;
bool show_onces_density = true;
int memory_opt = 0;
int problem_2d_plane_dm;
int anz_nval = 0;
int anz_nval0 = 0;                                //WW
//
int size_eval = 0;                                  //WW


bool hasAnyProcessDeactivatedSubdomains = false;  //NW

//--------------------------------------------------------
// Coupling Flag. WW
bool T_Process = false;					// Heat
bool H_Process = false;					// Fluid
bool H2_Process = false;				// Multi-phase
bool H3_Process = false;				// 3-phase
bool M_Process = false;					// Mechanical
bool RD_Process = false;				// Richards
bool MH_Process = false;				// MH monolithic scheme
bool MASS_TRANSPORT_Process = false;	// Mass transport
bool FLUID_MOMENTUM_Process = false;	// Momentum
bool RANDOM_WALK_Process = false;		// RWPT
bool PTC_FLOW_Process = false;			// PTC
bool pcs_created = false;
//
int pcs_number_deformation = -1;				// JT2012
int pcs_number_flow = -1;						// JT2012
int pcs_number_heat = -1;						// JT2012
vector<int>pcs_number_mass;						// JT2012

#define noCHECK_EQS
#define noCHECK_ST_GROUP
#define noCHECK_BC_GROUP

extern size_t max_dim;                            //OK411 todo

//////////////////////////////////////////////////////////////////////////
// PCS vector
//////////////////////////////////////////////////////////////////////////
// It is better to have space between data type and data name. WW
vector<CRFProcess*> pcs_vector;
//vector<string> ele_val_name_vector; // PCH
template <class T> T* resize(T* array, size_t old_size, size_t new_size);
//////////////////////////////////////////////////////////////////////////
// Construction / destruction
//////////////////////////////////////////////////////////////////////////

/**************************************************************************/
/* ROCKFLOW - Funktion: GetUncommentedLine
 */
/* Aufgabe:
   R���ckgabe ist ist ein string mit dem Zeileninhalt ab dem ersten Nicht-Leerzeichen
   bis zum ersten Auftreten des Kommentartzeichens ";"
   Abgeleitet aus GetLineFromFile1()                                            */
/* Programmaenderungen:
    06/2009     SB  First Version
 **************************************************************************/
std::string GetUncommentedLine(std::string line)
{
	std::string zeile = "";
	int i = 0, j = 0;
	//----------------------------------------------------------------------
	i = (int) line.find_first_not_of(" ",0); //Anf���ngliche Leerzeichen ���berlesen, i=Position des ersten Nichtleerzeichens im string
	j = (int) line.find(";",i);           //Nach Kommentarzeichen ; suchen. j = Position des Kommentarzeichens, j=-1 wenn es keines gibt.
	if((i != -1))
		zeile = line.substr(i,j - i);  //Ab erstem nicht-Leerzeichen bis Kommentarzeichen rauskopieren in neuen substring, falls Zeile nicht leer ist
	i = (int) zeile.find_last_not_of(" "); // Suche nach dem letzten Zeichen, dass kein Leerzeichen ist
	if(i >= 0)
	{
		line = zeile.substr(0,i + 1); // Leerzeichen am Ende rausschneiden
		zeile = line;
	}

	return zeile;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2004 OK Implementation
   08/2004 WW Read the deformation process
           Check the comment key '//' in .pcs
   last modified:
   10/2010 TF changed process type handling from string to enum
**************************************************************************/
bool PCSRead(std::string file_base_name)
{
	//----------------------------------------------------------------------
	char line[MAX_ZEILE];
	int indexCh1a, indexCh2a;
	std::string CommentK("//");
	std::string line_string;
	std::string pcs_file_name;
	std::ios::pos_type position;
	//========================================================================
	// File handling
	pcs_file_name = file_base_name + PCS_FILE_EXTENSION;
	std::ifstream pcs_file(pcs_file_name.data(), ios::in);
	if (!pcs_file.good())
	{
		cout << "Warning: no PCS data *.pcs file is missing" << endl;
		return false;
	}

	// rewind the file
	pcs_file.clear();
	pcs_file.seekg(0, std::ios::beg);
	//========================================================================
	// Keyword loop
	std::cout << "PCSRead ... " << std::flush;
	while (!pcs_file.eof())
	{
		pcs_file.getline(line, MAX_ZEILE);
		line_string = line;
		line_string = GetUncommentedLine(line_string);
		if (line_string.find("#STOP") != string::npos)
			break;
		indexCh1a = (int) line_string.find_first_of(CommentK.c_str());
		indexCh2a = (int) line_string.find("#PROCESS");
		//----------------------------------------------------------------------
		// keyword found
		if (indexCh2a > indexCh1a && (indexCh1a == -1))
		{
			CRFProcess* m_pcs = new CRFProcess();
			m_pcs->file_name_base = file_base_name; //OK
			position = m_pcs->Read(&pcs_file);
			m_pcs->PCSReadConfigurations(); //JT
			m_pcs->pcs_number = pcs_vector.size();

			//RelocateDeformationProcess(m_pcs);
			//			if (m_pcs->_pcs_type_name.find("DEFORMATION") != string::npos) { // TF
			if (isDeformationProcess(m_pcs->getProcessType()))
			{
				pcs_vector.push_back(m_pcs->CopyPCStoDM_PCS());
				pcs_vector[pcs_vector.size() - 1]->pcs_number = pcs_vector.size();
				delete m_pcs;
			}
			else{
				pcs_vector.push_back(m_pcs);
			}

			pcs_file.seekg(position, std::ios::beg);
		}                         // keyword found
	}                                     // eof

	std::cout << "done, read " << pcs_vector.size() << " processes" << std::endl;

	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: PCS read function
   Programing:
   06/2004 OK Implementation
   08/2004 WW Read deformation process
   11/2004 OK file streaming
   12/2005 OK MSH_TYPE
   01/2006 OK GEO_TYPE
**************************************************************************/
std::ios::pos_type CRFProcess::Read(std::ifstream* pcs_file)
{
	char line[MAX_ZEILE];
	string line_string;
	string CommentK("//");
	string hash("#");
	bool new_keyword = false;
	bool new_subkeyword = false;
	ios::pos_type position;
	ios::pos_type position_subkeyword;
	std::stringstream line_stream;
	saturation_switch = false;            // JOD for Richards
	//----------------------------------------------------------------------
	while (!new_keyword)
	{
		position = pcs_file->tellg();
		pcs_file->getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			break;
		}
		//....................................................................
		// WW Comment line
		if (line_string.find_first_of(CommentK.c_str()) != string::npos)
			return position;
		//SB check for comment sign ;
		line_string = GetUncommentedLine(line_string);
		//....................................................................
		// subkeyword found
		if (line_string.find("$PCS_TYPE") != string::npos)
			while ((!new_keyword) || (!new_subkeyword) || (!pcs_file->eof()))
			{
				position = pcs_file->tellg();
				line_string = readNonBlankLineFromInputStream(*pcs_file);
				if (line_string.find("#") != string::npos)
					return position;
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				line_stream.str(line_string);
				std::string pcs_type_name;
				line_stream >> pcs_type_name;
				pcs_type_name_vector.push_back(pcs_type_name);
				this->setProcessType (FiniteElement::convertProcessType(pcs_type_name));
				line_stream.clear();

				if (isFlowProcess (this->getProcessType()))
				{
					H_Process = true;
					this->isPCSFlow = true; //JT2012
					pcs_number_flow = pcs_vector.size(); //JT2012
					if(this->getProcessType() == FiniteElement::PS_GLOBAL ||
					   this->getProcessType() == FiniteElement::MULTI_PHASE_FLOW ||
					   pcs_type_name.find("H2") != string::npos){
						   this->isPCSMultiFlow = true;
					}
				}
				if (isDeformationProcess(this->getProcessType()))
				{
					M_Process = true;
					this->isPCSDeformation = true; //JT2012
					//JT: "pcs_number_deformation" is set in CRFProcessDeformation::Initialization()
				}
				if (this->getProcessType() == FiniteElement::MASS_TRANSPORT)
				{
					H_Process = true;
					MASS_TRANSPORT_Process = true;
					this->isPCSMass = true; //JT2012
					pcs_number_mass[pcs_no_components] = pcs_vector.size(); //JT2012
					pcs_no_components++;
					this->setProcessPrimaryVariable(FiniteElement::CONCENTRATION);
				}
				if (this->getProcessType() == FiniteElement::PTC_FLOW){
					PTC_FLOW_Process = true;
				}
				if (this->getProcessType() == FiniteElement::HEAT_TRANSPORT){
					T_Process = true;
					this->isPCSHeat = true; //JT2012
					pcs_number_heat = pcs_vector.size(); //JT2012
				}
				if (this->getProcessType() == FiniteElement::FLUID_MOMENTUM){
					FLUID_MOMENTUM_Process = true;
				}
				if (this->getProcessType() == FiniteElement::RANDOM_WALK){
					RANDOM_WALK_Process = true;
				}
			}
		//....................................................................
		// subkeyword found
		if (line_string.find("$NUM_TYPE") != string::npos)
		{
			*pcs_file >> num_type_name;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$CPL_TYPE") != string::npos)
		{
			*pcs_file >> cpl_type_name;
			if (cpl_type_name.compare("MONOLITHIC") == 0)
			{
				pcs_monolithic_flow = true;
				pcs_deformation = 11;
			}
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$TIM_TYPE") != string::npos)
		{
			*pcs_file >> tim_type_name;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$APP_TYPE") != string::npos)
		{
			*pcs_file >> rwpt_app;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$PRIMARY_VARIABLE") != string::npos)
		{
			*pcs_file >> primary_variable_name;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$ELEMENT_MATRIX_OUTPUT") != string::npos)
		{
			*pcs_file >> Write_Matrix; //WW
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		//WW
		if (line_string.find("$BOUNDARY_CONDITION_OUTPUT") != string::npos)
		{
			write_boundary_condition = true;
			continue;
		}
		//....................................................................
		//BG 05/2012
		if (line_string.find("$OutputMassOfGasInModel") != string::npos)
		{
			OutputMassOfGasInModel = true;
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$ST_RHS") != string::npos)
		{
			*pcs_file >> WriteSourceNBC_RHS; //WW
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		if (line_string.find("$PROCESSED_BC") != string::npos) //25.08.2011. WW
		{
			*pcs_file >> WriteProcessed_BC;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}

		//....................................................................
		// subkeyword found
		if (line_string.find("$MEMORY_TYPE") != string::npos)
		{
			*pcs_file >> Memory_Type; //WW
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$RELOAD") != string::npos)
		{
			*pcs_file >> reload; //WW
			if (reload == 1 || reload == 3)
				*pcs_file >> nwrite_restart;  //kg44 read number of timesteps between writing restart files
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		// subkeyword found
		if (line_string.find("$DEACTIVATED_SUBDOMAIN") != string::npos)
		{
			//WW
			*pcs_file >> NumDeactivated_SubDomains >> ws;
			Deactivated_SubDomain = new int[NumDeactivated_SubDomains];
			for (int i = 0; i < NumDeactivated_SubDomains; i++)
				*pcs_file >> Deactivated_SubDomain[i] >> ws;
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$MSH_TYPE") != string::npos)
		{
			*pcs_file >> msh_type_name >> ws;
			continue;
		}
		//....................................................................
		//		if (line_string.find("$GEO_TYPE") != string::npos) { //OK
		//			*pcs_file >> geo_type >> geo_type_name >> ws;
		//			continue;
		//		}
		//
		//....................................................................
		// subkeyword found
		if (line_string.find("$MEDIUM_TYPE") != string::npos)
		{
			while ((!new_keyword) || (!new_subkeyword) || (!pcs_file->eof()))
			{
				position_subkeyword = pcs_file->tellg();
				*pcs_file >> line_string;
				if (line_string.size() == 0)
					break;
				if (line_string.find("#") != string::npos)
				{
					new_keyword = true;
					break;
				}
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.find("CONTINUUM") != string::npos)
				{
					*pcs_file >> line_string;
					//WW
					double w_m = strtod(line_string.data(), NULL);
					continuum_vector.push_back(w_m);
					//WW
					continuum_vector.push_back(1.0 - w_m);
					break; //WW
				}
				pcs_file->ignore(MAX_ZEILE, '\n');
			}
			continue;
		}
		//OK
		if (line_string.find("$SATURATION_SWITCH") != string::npos)
		{
			saturation_switch = true;
			continue;
		}
		//SB4900
		if (line_string.find("$USE_VELOCITIES_FOR_TRANSPORT") != string::npos)
		{
			//// Only for fluid momentum process
			//if (this->getProcessType () == FiniteElement::FLUID_MOMENTUM)
			//	use_velocities_for_transport = true;
			continue;
		}
		//Interface to Eclipse and Dumux, BG, 09/2010
		//	if(line_string.find("$SIMULATOR")!=string::npos) { //OK
		if(line_string.compare("$SIMULATOR") == 0) // BG, 09/2010, coupling to Eclipse and DuMux
		{
			*pcs_file >> this->simulator;
			continue;
		}
		if(line_string.find("$SIMULATOR_PATH") == 0) // BG, 09/2010, coupling to Eclipse and DuMux
		{
			*pcs_file >> this->simulator_path;
			continue;
		}
		// BG, 09/2010, coupling to Eclipse and DuMux
		if(line_string.find("$SIMULATOR_MODEL_PATH") == 0)
		{
			*pcs_file >> this->simulator_model_path;
			continue;
		}
		// BG, 09/2010, coupling to Eclipse and DuMux
		if(line_string.find("$USE_PRECALCULATED_FILES") == 0)
		{
			this->PrecalculatedFiles = true;
			continue;
		}
		// KB, 02/2011, coupling to Eclipse and DuMux
		if(line_string.find("$SIMULATOR_WELL_PATH") == 0)
		{
			*pcs_file >> this->simulator_well_path;
			continue;
		}
		// BG, NB 11/2010, calculating phase transition for CO2
		if(line_string.find("$PHASE_TRANSITION") == 0)
		{
			string tempstring;
			*pcs_file >> tempstring;
			//if (tempstring == "CO2_H2O_NaCl")
			//	this->Phase_Transition_Model = 1;
			continue;
		}
		//WX:07.2011
		if(line_string.find("$TIME_CONTROLLED_EXCAVATION") == 0)
		{
			*pcs_file >> ExcavMaterialGroup >> ExcavDirection >>
			ExcavBeginCoordinate >> ExcavCurve;
			continue;
		}
		//....................................................................
	}
	//----------------------------------------------------------------------
	return position;
}

