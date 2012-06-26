/**************************************************************************
   FEMLib - Class: BC BoundaryConditions
   Task:
   Programing:
   02/2004 OK Implementation
   last modified
**************************************************************************/
#include "rf_bc_new.h"

// C++ STL
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <list>
#include <sstream>

#include "makros.h"
//// FileIO
////#include "BoundaryConditionIO.h"
//#include "GeoIO.h"
#include "ProcessIO.h"
#include "readNonBlankLineFromInputStream.h"
//
//// GEOLib
////#include "geo_lib.h"
////#include "geo_sfc.h"
//
//// GEOLIB
//#include "GEOObjects.h"
//
//// MSHLib
////#include "mshlib.h"
//// FEMLib
//extern void remove_white_space(std::string*);
////#include "problem.h"
//#include "gs_project.h"
//#include "tools.h"
////#include "rf_node.h"
////#include "rf_pcs.h"
////#include "rf_fct.h"
//#include "rfmat_cp.h"
////#include "geo_ply.h"
//// MathLib
//#include "InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
//#include "mathlib.h"
//
//#include "BoundaryCondition.h"

#ifndef _WIN32
#include <cstdio>
#include <cstdlib>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

double cputime(double x)
{
	struct rusage rsrc;
	double usr, sys;

	if (getrusage(RUSAGE_SELF, &rsrc) == -1)
	{
		perror("times");
		exit(1);
	}

	usr = rsrc.ru_utime.tv_sec + 1.0e-6 * rsrc.ru_utime.tv_usec;
	sys = rsrc.ru_stime.tv_sec + 1.0e-6 * rsrc.ru_stime.tv_usec;

	return usr + sys - x;
}
#endif


//==========================================================================
std::list<CBoundaryCondition*> bc_list;
std::vector<std::string> bc_db_head;
std::list<CBoundaryConditionsGroup*> bc_group_list;
std::vector<CBoundaryCondition*> bc_db_vector;

/**************************************************************************
   FEMLib-Method:
   Task: BC read function
   Programing:
   01/2004 OK Implementation
   09/2004 OK POINTS method
   11/2004 MX stream string
**************************************************************************/
std::ios::pos_type CBoundaryCondition::Read(std::ifstream* bc_file,
//                                            const GEOLIB::GEOObjects& geo_obj,
                                            const std::string& unique_fname,
                                            bool & valid)
{
	std::string line_string;
	bool new_keyword = false;
	std::ios::pos_type position;

	std::string sub_string, strbuff;
	int ibuff;                            //pos,
	double dbuff;                         //WW
	std::stringstream in;

    std::string FilePath;

	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		position = bc_file->tellg();
		line_string = readNonBlankLineFromInputStream (*bc_file);
		if (line_string.size() < 1)
			break;
		if (line_string.find("#") != std::string::npos)
		{
			new_keyword = true;
			break;
		}

		if (line_string.find("$PCS_TYPE") != std::string::npos)
			if (!FileIO::ProcessIO::readProcessInfo (*bc_file, _pcs_type))
				valid = false;

		if (line_string.find("$PRIMARY_VARIABLE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			std::string tmp;
			in >> tmp;    // _pcs_pv_name;
			in.clear();
		}

		// HS, this is new. later on we should stick to COMP_NAME, PRIMARY_VARIABLE support will be removed.
		if (line_string.find("$COMP_NAME") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			std::string tmp;
			in >> tmp;    // _pcs_pv_name;
			in.clear();
		}

		if (line_string.find("$GEO_TYPE") != std::string::npos)
			//if (!FileIO::GeoIO::readGeoInfo (this, *bc_file, geo_name, geo_obj,
			//                                 unique_fname))
			//	valid = false;
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            std::string tmp;
            in >> tmp;    // _pcs_pv_name;
            in.clear();
        }

		//PCH
		if (line_string.find("$DIS_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> line_string; //sub_line
			_periodic = false; // JOD

			// Source terms are assign to element nodes directly. 23.02.2009. WW
			if (line_string.find("DIRECT") != std::string::npos)
			{
				this->setProcessDistributionType(FiniteElement::DIRECT);
				in >> fname;
				fname = FilePath + fname;
				in.clear();
			}

			if (line_string.find("CONSTANT") != std::string::npos)
			{
				this->setProcessDistributionType(FiniteElement::CONSTANT);
				in >> geo_node_value; //sub_line
				in.clear();
			}
			// If a linear function is given. 25.08.2011. WW
			if (line_string.find("FUNCTION") != std::string::npos)
			{
				setProcessDistributionType(FiniteElement::FUNCTION);
				in.clear();
				dis_linear_f = new LinearFunctionData(*bc_file);
			}
			if (line_string.find("LINEAR") != std::string::npos)
			{
				this->setProcessDistributionType(FiniteElement::LINEAR);
				// Distribuded. WW
				size_t nLBC;
				in >> nLBC; //sub_line
				in.clear();

				for (size_t i = 0; i < nLBC; i++)
				{
					in.str(readNonBlankLineFromInputStream(*bc_file));
					in >> ibuff >> dbuff >> strbuff;
					in.clear();

					//           *bc_file>>ibuff>>dbuff;
					_PointsHaveDistribedBC.push_back(ibuff);
					_DistribedBC.push_back(dbuff);
					if (strbuff.size() > 0)
					{
						_PointsFCTNames.push_back(strbuff);
						time_dep_interpol = true;
					}
				}
				//        bc_file->ignore(MAX_ZEILE,'\n');
			}
		}

		// Time dependent function
		//..Time dependent curve ............................................
		if (line_string.find("$TIM_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> line_string;

			if (line_string.find("CURVE") != std::string::npos)
			{
				//				tim_type_name = "CURVE";
				this->setProcessDistributionType(FiniteElement::CONSTANT);
				in >> _curve_index;
				in.clear();

				//        pos1=pos2+1;
				//        sub_string = get_sub_string(buffer,"  ",pos1,&pos2);
				//		_curve_index = atoi(sub_string.c_str());
			}
			continue;
		}

		if (line_string.find("$FCT_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> fct_name; //sub_line
			in.clear();
		}

		if (line_string.find("$MSH_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> sub_string; //sub_line
			_msh_type_name = "NODE";
			if (sub_string.find("NODE") != std::string::npos)
			{
				in >> _msh_node_number;
				in.clear();
			}
		}

		if (line_string.find("$DIS_TYPE_CONDITION") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file)); // CONSTANT -21500.0
			in >> line_string;
			if (line_string.find("CONSTANT") != std::string::npos)
			{
				this->setProcessDistributionType(FiniteElement::CONSTANT);
				in >> geo_node_value;
				in.clear();
			}
			in.str(readNonBlankLineFromInputStream(*bc_file)); // 0.0 IF HEAD > 0.04
			std::string pcs_pv_name_cond; // 07/2010 TF temp string
			in >> node_value_cond >> line_string >> pcs_pv_name_cond
			>> line_string >> condition;
			in.clear();
			in.str(readNonBlankLineFromInputStream(*bc_file)); // PCS OVERLAND_FLOW
			std::string pcs_type_name_cond;
			in >> line_string >> pcs_type_name_cond;
			in.clear();
			conditional = true;
		}

		if (line_string.find("$EPSILON") != std::string::npos) // NW
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> epsilon;
			in.clear();
		}
		//....................................................................
		//aktive state of the bc is time controlled  WX
		if (line_string.find("$TIME_CONTROLLED_ACTIVE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> time_contr_curve;
			in.clear();
		}
		//....................................................................
		//bc for excated boundaries WX
		if (line_string.find("$EXCAVATION") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> bcExcav >> MatGr;
			in.clear();
		}
		//....................................................................
	}
	return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: BC constructor
   Programing:
   01/2004 OK Implementation
**************************************************************************/
CBoundaryCondition::CBoundaryCondition() :
	geo_name (""), _curve_index (-1), dis_linear_f(NULL)
{
	this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
	// FCT
	conditional = false;
	time_dep_interpol = false;
	epsilon = 1e-9;                       //NW
	time_contr_curve = -1;                //WX
	bcExcav = -1;                         //WX
	MatGr = -1;                           //WX
}



/**************************************************************************
   FEMLib-Method:
   Task: BC deconstructor
   Programing:
   01/2004 OK Implementation
**************************************************************************/
CBoundaryCondition::~CBoundaryCondition()
{
	// DIS
	node_number_vector.clear();
	geo_node_number = -1;
	geo_node_value = 0.0;

	//WW
	if(dis_linear_f)
		delete dis_linear_f;
	dis_linear_f = NULL;
}

const std::string& CBoundaryCondition::getGeoName () const
{
	return geo_name;
}


/**************************************************************************
   FEMLib-Method:
   Task: BC read function
   Programing:
   01/2004 OK Implementation
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
   05/2010 TF changes due to new GEOLIB integration, some improvements
**************************************************************************/
bool BCRead(std::string const& file_base_name, 
        //const GEOLIB::GEOObjects& geo_obj,
            const std::string& unique_name)
{
	char line[MAX_ZEILE];
	std::string line_string, bc_file_name;

	// File handling
	bc_file_name = file_base_name + BC_FILE_EXTENSION;

	std::ifstream bc_file(bc_file_name.data(), std::ios::in);
	if (!bc_file.good())
	{
		std::cout << "! Error in BCRead: No boundary conditions !" << std::endl;
		return false;
	}

	// Keyword loop
	std::cout << "BCRead ... " << std::flush;
	while (!bc_file.eof())
	{
		bc_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != std::string::npos)
		{
			std::cout << "done, read " << bc_list.size()
			          << " boundary conditions" << std::endl;
			return true;
		}
		if (line_string.find("#BOUNDARY_CONDITION") != std::string::npos)
		{
			CBoundaryCondition* bc(new CBoundaryCondition());
			bool valid (true);
			std::ios::pos_type position = bc->Read(&bc_file, unique_name, valid);
			if (valid)
				bc_list.push_back(bc);
			else
				delete bc;
			bc_file.seekg(position, std::ios::beg);
		} // keyword found
	} // eof
	return true;
}

///**************************************************************************
//   FEMLib-Method: BCWrite
//   Task: master write function
//   Programing:
//   02/2004 OK Implementation
//   last modification:
//**************************************************************************/
//
//void BCWrite(std::string const& base_file_name)
//{
//	std::string sub_line;
//	std::string line_string;
//
//	// File handling
//	std::string bc_file_name (base_file_name + BC_FILE_EXTENSION);
//	std::fstream bc_file(bc_file_name.data(), std::ios::trunc | std::ios::out);
//	bc_file.setf(std::ios::scientific, std::ios::floatfield);
//	bc_file.precision(12);
//	//OK string tec_file_name = base_file_name + ".tec";
//	//OK fstream tec_file (tec_file_name.data(),ios::trunc|ios::out);
//	//OK tec_file.setf(ios::scientific,ios::floatfield);
//	//OK tec_file.precision(12);
//	if (!bc_file.good())
//		return;
//	bc_file.seekg(0L, std::ios::beg); // rewind?
//	bc_file <<
//	"GeoSys-BC: Boundary Conditions ------------------------------------------------\n";
//	//========================================================================
//	// BC list
//	std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
//	while (p_bc != bc_list.end())
//	{
//		FileIO::BoundaryConditionIO::write(bc_file, *(*p_bc));
//		++p_bc;
//	}
//	bc_file << "#STOP";
//	bc_file.close();
//	//OK tec_file.close();
//}

/**************************************************************************
   FEMLib-Method:
   01/2004 OK Implementation
   07/2007 OK V2, global function
**************************************************************************/
//CBoundaryCondition* BCGet(const std::string &pcs_name, const std::string &geo_type_name,
//		const std::string &geo_name)
//{
//	std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
//	while (p_bc != bc_list.end()) {
//		if (((*p_bc)->pcs_type_name.compare(pcs_name) == 0)
//				&& ((*p_bc)->geo_type_name.compare(geo_type_name) == 0)
//				&& ((*p_bc)->getGeoName().compare(geo_name) == 0))
//			return *p_bc;
//		++p_bc;
//	}
//	return NULL;
//}


CBoundaryConditionsGroup::CBoundaryConditionsGroup(void)
{
	msh_node_number_subst = -1;           //
	time_dep_bc = -1;
}

