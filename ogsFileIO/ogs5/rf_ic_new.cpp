/**************************************************************************
   FEMLib - Object: Initial Conditions IC
   Task:
   Programing:
   08/2004 OK Implementation
   12/2005 OK Restart
   last modified
**************************************************************************/

//#include "gs_project.h"
#include "makros.h"
// C++ STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
using namespace std;
// FEM-Makros
//#include "files0.h"
//#include "files0.h"
//#include "mathlib.h"
// Base
//#include "StringTools.h"
// GEOLIB
//#include "GEOObjects.h"
// MSHLib
//#include "msh_lib.h"
// FEMLib
#include "rf_ic_new.h"
//#include "rf_node.h"
//#include "rf_pcs.h"
//#include "rfmat_cp.h"

//#include "InitialCondition.h"
#include "readNonBlankLineFromInputStream.h"

//==========================================================================
vector<CInitialConditionGroup*>ic_group_vector;
vector<CInitialCondition*>ic_vector;

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
**************************************************************************/
CInitialCondition::CInitialCondition() 
{
	//  geo_type_name = "DOMAIN";
	this->setProcessDistributionType(FiniteElement::CONSTANT);
	// HS: not needed, removed.
	// m_node = new CNodeValue();
	// m_node->node_value = 0.0;
	SubNumber = 0;
	this->setProcess(NULL);               //OK

}



/**************************************************************************
   FEMLib-Method:
   Task: BC deconstructor
   Programing:
   04/2004 OK Implementation
**************************************************************************/
CInitialCondition::~CInitialCondition(void)
{
/* KR
   if ( !node_value_vector.empty() )
   {
      for (size_t i=0; i < node_value_vector.size(); i++)
         delete node_value_vector[i];
   }
   node_value_vector.clear();
 */
}

/**************************************************************************
   FEMLib-Method:
   Task: IC read function
   Programing:
   08/2004 OK Implementation
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
   05/2010 TF reformated, restructured, signature changed, use new GEOLIB data structures
**************************************************************************/
bool ICRead(const std::string& file_base_name,
            //const GEOLIB::GEOObjects& geo_obj, 
            const std::string& unique_name)
{
	// File handling
	std::string ic_file_name = file_base_name + IC_FILE_EXTENSION;
	std::ifstream ic_file(ic_file_name.data(), std::ios::in);
	if (!ic_file.good())
	{
		std::cout << "WARNING: ICRead: No initial conditions !" << std::endl;
		return false;
	}

	char line[MAX_ZEILE];
	std::string line_string;
	std::ios::pos_type position;

	// Keyword loop
	std::cout << "ICRead" << std::endl;
	while (!ic_file.eof())
	{
		ic_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != string::npos)
			return true;

		// keyword found
		if (line_string.find("#INITIAL_CONDITION") != std::string::npos)
		{
			CInitialCondition* ic = new CInitialCondition();
			std::ios::pos_type pos (ic_file.tellg());
			position = ic->Read(&ic_file, unique_name);
			if (pos != position)
				ic_vector.push_back(ic);
			else
			{
				std::cerr <<
				"WARNING: in ICRead: could not read initial condition" << std::endl;
				delete ic;
			}
			ic = NULL;
			ic_file.seekg(position, ios::beg);
		}                         // keyword found
	}                                     // eof
	return true;
}


/**************************************************************************
   FEMLib-Method:
   Task: ST read function
   Programing:
   08/2004 OK Implementation (DOMAIN)
   10/2004 OK POLYLINE
   11/2004 CMCD String streaming, Gradient
   01/2005 WW Initialize by sub-domain
   05/2005 OK PRIMARY_VARIABLE
   12/2005 OK RESTART
   07/2006 WW Read data by expression
   06/2010 TF changed signature to use the new GEOLIB data structures
**************************************************************************/
ios::pos_type CInitialCondition::Read(std::ifstream* ic_file,
                                      //const GEOLIB::GEOObjects& geo_obj,
                                      const std::string& unique_geo_name)
{
	string line_string;
	std::stringstream in;
	ios::pos_type position;

	bool new_keyword = false;
	int ibuf (0);
	double d_buf (0.0);

	// read loop
	while (!new_keyword)
	{
		position = ic_file->tellg();
		line_string = readNonBlankLineFromInputStream(*ic_file);
		if (line_string.size() < 1)
			break;
		if (line_string.find("#") != string::npos)
		{
			new_keyword = true;
			break;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$PCS_TYPE") != string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*ic_file));
			std::string tmp;
			in >> tmp;    // pcs_type_name;
			this->setProcessType(FiniteElement::convertProcessType(tmp));
			in.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$PRIMARY_VARIABLE") != string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*ic_file));
			std::string tmp;
			in >> tmp;    // pcs_pv_name;
			in.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$COMP_NAME") != string::npos)
		{
            in.str(readNonBlankLineFromInputStream(*ic_file));
			std::string tmp;
			in >> tmp;
			in.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$DIS_TYPE") != string::npos)
		{
            in.str(readNonBlankLineFromInputStream(*ic_file));
			std::string tmp;
			in >> tmp;    // dis_type_name;
			this->setProcessDistributionType(FiniteElement::convertDisType(tmp));

			if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
				//KR CNodeValue* m_node = new CNodeValue();
				//KR node_value_vector.push_back(m_node);
				in >> geo_node_value;
			else if (this->getProcessDistributionType() == FiniteElement::GRADIENT)
			{
				in >> gradient_ref_depth; //CMCD
				in >> gradient_ref_depth_value; //CMCD
				in >> gradient_ref_depth_gradient; //CMCD
			}
			else if (this->getProcessDistributionType() == FiniteElement::RESTART)
				in >> rfr_file_name;
			else if (this->getProcessDistributionType() == FiniteElement::DIRECT)
			{
				in >> fname;
				fname = FilePath + fname;
			}
			in.clear();
			continue;
		}
		//....................................................................
		//subkeyword found
		if (line_string.find("$GEO_TYPE") != string::npos)
		{
            in.str(readNonBlankLineFromInputStream(*ic_file));
			std::string geo_type_name;
			in >> geo_type_name;
			if (geo_type_name.find("POINT") != string::npos)
			{
				in >> geo_name;
				in.clear();
				geo_name = ""; // REMOVE CANDIDATE
			}
			if (geo_type_name.find("POLYLINE") != string::npos)
			{
				in >> geo_name;
			}
			if (geo_type_name.find("SURFACE") != string::npos)
			{
				in >> geo_name;
			}
			//if (geo_type_name.find("VOLUME") != string::npos)
			if (geo_type_name.find("DOMAIN") != string::npos)
			{
				//  Give initial condition by patches of domain. WW
				if (geo_type_name.find("SUB") != string::npos)
				{
					*ic_file >> SubNumber;
					for (size_t i = 0; i < SubNumber; i++)
					{
						*ic_file >> ibuf >> d_buf;
						subdom_index.push_back(ibuf);
						subdom_ic.push_back(d_buf);
					}
				}
			}
			in.clear();
			continue;
		}
	}                                     // Schleife ueber alle Phasen bzw. Komponenten
	return position;
}

