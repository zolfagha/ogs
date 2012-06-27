/**************************************************************************
 FEMLib - Object: Source Terms ST
 Task:
 Programing:
 01/2004 OK Implementation
 last modified
 **************************************************************************/
#include "rf_st_new.h"

#include "makros.h"
// C++ STL
#include <fstream>
#include <cfloat>
#include <iostream>
#include <sstream>
#include <set>
#include <cstdlib>

#include "readNonBlankLineFromInputStream.h"

//#include "files0.h"
//#include "mathlib.h"
//
//// GeoSys-GeoLib
//#include "GEOObjects.h"
//// GEO
////#include "GeoType.h"
//
//// GeoSys-MshLib
//#include "fem_ele.h"
//
//#include "tools.h"                                //GetLineFromFile
///* Tools */
//#ifndef NEW_EQS                                   //WW. 06.11.2008
//#include "matrix.h"
//#endif
//
//// GeoSys-FEMLib
////OK_IC #include "rfsousin.h"
//#include "rf_tim_new.h"
//
//// Math
//#include "matrix_class.h"
//
////#include "pcs_dm.h"
//
//// FEM
////#include "problem.h"
//// For analytical source terms
//#include "rf_mfp_new.h"
//#include "rf_node.h"
//#include "rfmat_cp.h"
//
//// Base
//#include "quicksort.h"
//
//// MathLib
//#include "InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
//
//// FileIO
//#include "FEMIO/GeoIO.h"
//#include "FEMIO/ProcessIO.h"
//#include "readNonBlankLineFromInputStream.h"
//
//#include "SourceTerm.h"
//
//using FiniteElement::CElement;
//using MeshLib::CElem;
//using MeshLib::CEdge;
//using MeshLib::CNode;
//using Math_Group::vec;

#ifndef GRAVITY_CONSTANT
#define GRAVITY_CONSTANT 9.81
#endif

std::vector<CSourceTerm*> st_vector;
std::list<CSourceTermGroup*> st_group_list;
std::vector<std::string> analytical_processes;
std::vector<std::string> analytical_processes_polylines;
//std::vector<NODE_HISTORY*> node_history_vector;   //CMCD
/**************************************************************************
 FEMLib-Method:
 Task: ST constructor
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
CSourceTerm::CSourceTerm() :
	ProcessInfo(), _coupled (false), _sub_dom_idx(-1), GIS_shape_head(NULL)
                                                  // 07.06.2010, 03.2010. WW
{
   CurveIndex = -1;
   //KR critical_depth = false;
   //	COUPLING_SWITCH = false;
   geo_node_value = 0.0;
   nodes = NULL;                                  //OK
   analytical = false;                            //CMCD
   //  display_mode = false; //OK
   this->TimeInterpolation = 0;                   //BG
}


/**************************************************************************
 FEMLib-Method:
 Task: BC deconstructor
 Programing:
 04/2004 OK Implementation
 **************************************************************************/
CSourceTerm::~CSourceTerm()
{
   DeleteHistoryNodeMemory();
   //    dis_file_name.clear();
   node_number_vector.clear();
   node_value_vector.clear();
   node_renumber_vector.clear();
   PointsHaveDistribedBC.clear();
   DistribedBC.clear();
   element_st_vector.clear();
   //WW----------22.02.2007-------------------
   // TF 06/2010
   size_t size(normal2surface.size());
   for (size_t i = 0; i < size; i++)
      delete normal2surface[i];
   size = pnt_parameter_vector.size();
   for (size_t i = 0; i < size; i++)
      delete pnt_parameter_vector[i];
   if(GIS_shape_head)                             // 07.06.2010. WW
   {
      delete [] GIS_shape_head;
      GIS_shape_head = NULL;
   }
   //WW

   //WW---------------------------------------
}


const std::string& CSourceTerm::getGeoName() const
{
   return geo_name;
}


double CSourceTerm::getCoupLeakance () const
{
   return _coup_leakance;
}


/**************************************************************************
 FEMLib-Method:
 Task: ST read function
 Programing:
 01/2004 OK Implementation
 11/2004 MB neues read Konzept
 02/2005 MB River condition
 03/2005 WW Node force released by excavation
 11/2005 CMCD Analytical source for matrix
 04/2006 OK CPL
 04/2006 OK MSH_TYPE
06/2010 TF modification of the signature, added geo_obj and unique_name
**************************************************************************/
std::ios::pos_type CSourceTerm::Read(std::ifstream *st_file,
		//const GEOLIB::GEOObjects& geo_obj, 
        const std::string& unique_name)
{
   char line[MAX_ZEILE];
   std::string line_string, sub_string;
   bool new_keyword = false;

   std::stringstream in;

                                                  // JOD 4.10.01
   channel = 0, node_averaging = 0, no_surface_water_pressure = 0;
   std::ios::pos_type position;

   // read loop
   while (!new_keyword)
   {
      position = st_file->tellg();
      line_string = readNonBlankLineFromInputStream(*st_file);
      if (line_string.find("#") != std::string::npos)
      {
         new_keyword = true;
         break;
      }

      /* search for keywords */
                                                  // subkeyword found
      if (line_string.find("$PCS_TYPE") != std::string::npos)
      {
    	  //FileIO::ProcessIO::readProcessInfo (*st_file, _pcs_type);
          in.str (readNonBlankLineFromInputStream (*st_file));
         std::string tmp;
         in >> tmp;
         //setProcessType (convertProcessType (tmp));
         in.clear();
         continue;
      }

                                                  // subkeyword found
      if (line_string.find("$PRIMARY_VARIABLE") != std::string::npos)
      {
    	  in.str (readNonBlankLineFromInputStream (*st_file));
//         in.str(readNonBlankLineFromInputStream(*st_file));
         std::string tmp;
         in >> tmp;

         if ( this->getProcessType() == FiniteElement::MASS_TRANSPORT )
         {
             // HS set the pointer to MCP based on component name.
             // first do a check whether this name is existing and unique.
             //if ( cp_name_2_idx.count( tmp ) == 1 )
             //{
             //    setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess() );
             //    setProcessPrimaryVariable( FiniteElement::CONCENTRATION );
             //}
             //else
             //{
             //    DisplayErrorMsg("Error: In reading ST file, the input component names are not found in MCP file!!!");
             //    exit(1);
             //}
         }
         else
         {
             //setProcess( PCSGet( this->getProcessType() ) );
             //setProcessPrimaryVariable (FiniteElement::convertPrimaryVariable (tmp));
         }
         in.clear();
         continue;
      }

      if (line_string.find("$COMP_NAME") != std::string::npos)
      {
    	  in.str(readNonBlankLineFromInputStream (*st_file));
//         in.str(readNonBlankLineFromInputStream(*st_file));
         std::string tmp;
         in >> tmp;
         // HS set the pointer to MCP based on component name.
         // first do a check whether this name is existing and unique.
         //if ( cp_name_2_idx.count( tmp ) == 1 )
         //{
         //    setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess() );
         //    setProcessPrimaryVariable( FiniteElement::CONCENTRATION );
         //}
         //else
         //{
         //    DisplayErrorMsg("Error: In reading ST file, the input component names are not found in MCP file!!!");
         //    exit(1);
         //}
         in.clear();
         continue;
      }

      if (line_string.find("$GEO_TYPE") != std::string::npos)
      {
         //ReadGeoType(st_file, geo_obj, unique_name);
         continue;
      }

                                                  //05.09.2008 WW
      if (line_string.find("$DIS_TYPE") != std::string::npos)
      {
         //10.04.2008. WW  if(line_string.compare("$DIS_TYPE")==0) {
         if (line_string.find("CONDITION") != std::string::npos)
         {
            _coupled = true;
            ReadDistributionType(st_file);
            in.str(readNonBlankLineFromInputStream(*st_file));
            in >> line_string >> pcs_type_name_cond;
            in.clear();
            in.str(readNonBlankLineFromInputStream(*st_file));    //
            in >> pcs_pv_name_cond;
            in.clear();
//            in.str(readNonBlankLineFromInputStream(*st_file));
            in.str(readNonBlankLineFromInputStream(*st_file));
            in >> _coup_leakance >> rill_height;
            in.clear();
         }                                        //05.09.2008 WW
         else
         {
            ReadDistributionType(st_file);
            continue;
         }
      }

      if (line_string.find("$NODE_AVERAGING") != std::string::npos)
      {
         in.clear();
         node_averaging = true;
         continue;
      }
                                                  // JOD 4.10.01
      if (line_string.find("$NEGLECT_SURFACE_WATER_PRESSURE") != std::string::npos)
      {
         in.clear();
         no_surface_water_pressure = true;
         continue;
      }

      if (line_string.find("$CHANNEL") != std::string::npos)
      {
         in.clear();
         in.str(readNonBlankLineFromInputStream(*st_file));
         in >> channel_width;
         channel = 1;
         continue;
      }

      if (line_string.find("$TIM_TYPE") != std::string::npos)
      {
         in.str(readNonBlankLineFromInputStream(*st_file));
         in >> tim_type_name;
         if (tim_type_name.find("CURVE") != std::string::npos)
         {
        	 //				dis_type = 0;
            in >> CurveIndex;
         }
         in.clear();
         continue;
      }

  	  //defines if time dependent source terms are use as piecewise constant or linear interpolated; BG 05/2011
      if (line_string.find("$TIME_INTERPOLATION") != std::string::npos)
      {
         in.str(readNonBlankLineFromInputStream(*st_file));
         in >> interpolation_method;
         if (interpolation_method.find("LINEAR") != std::string::npos)
         {
            this->TimeInterpolation = 0;
         }
         if (interpolation_method.find("PIECEWISE_CONSTANT") != std::string::npos)
         {
            this->TimeInterpolation = 1;
         }
         in.clear();
         continue;
      }

      if (line_string.find("$FCT_TYPE") != std::string::npos)
      {
         in.str(readNonBlankLineFromInputStream(*st_file));
         in >> fct_name;                          //sub_line
                                                  //WW
         if (fct_name.find("METHOD") != std::string::npos)
            in >> fct_method;
         in.clear();
      }

      if (line_string.find("$MSH_TYPE") != std::string::npos)
      {
         in.str(readNonBlankLineFromInputStream(*st_file));
         std::string sub_string;
         in >> sub_string;                        //sub_line
         msh_type_name = "NODE";
         if (sub_string.find("NODE") != std::string::npos)
         {
            in >> msh_node_number;
            in.clear();
         }
         continue;
      }
   }                                              // end !new_keyword
   return position;
}


/**************************************************************************
 FEMLib-Method:
 Task: for CSourceTerm::Read
 Programing:
 11/2007 JOD Implementation
 02/2009 WW  Add a functionality to directly assign source terms to element nodes.
 **************************************************************************/
void CSourceTerm::ReadDistributionType(std::ifstream *st_file)
{
   std::stringstream in;
   // 03.2010 WW
   std::string aline;
   std::stringstream ss;
   int abuff, nLBC = 0;
   double bbuff;

   std::string dis_type_name;
   in.str(readNonBlankLineFromInputStream(*st_file));
   in >> dis_type_name;

   this->setProcessDistributionType (FiniteElement::convertDisType(dis_type_name));

   if (dis_type_name.compare (convertDisTypeToString (this->getProcessDistributionType())) != 0)
   {
      std::cerr << "Error in CSourceTerm::ReadDistributionType (): dist_type_name #" << dis_type_name << "#, new: " << convertDisTypeToString (this->getProcessDistributionType()) << std::endl;
      exit (1);
   }

   if (   this->getProcessDistributionType() == FiniteElement::CONSTANT
       || this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN
       || this->getProcessDistributionType() == FiniteElement::CONSTANT_GEO      )
   {
      in >> geo_node_value;
      in.clear();
   }

   //	if (dis_type_name.find("ANALYTICAL") != std::string::npos) {
   if (this->getProcessDistributionType() == FiniteElement::ANALYTICAL)
   {
      in >> analytical_material_group;            //Which material group is it being applied to
      in >> analytical_diffusion;                 //D value
      in >> analytical_porosity;                  //n value of matrix
      in >> analytical_tortousity;                //t value of matrix
      in >> analytical_linear_sorption_Kd;        //Linear sorption coefficient
      in >> analytical_matrix_density;            //Density of solid
      in >> number_of_terms;                      //no timesteps to consider in solution
      in >> resolution;                           //every nth term will be considered
      in >> factor;                               //to convert temperature to energy
      analytical = true;
      analytical_processes.push_back(convertPrimaryVariableToString (getProcessPrimaryVariable()));
      //		if (geo_type_name.compare("POLYLINE") == 0)
      //if (this->getGeoType() == GEOLIB::POLYLINE)
      //   analytical_processes_polylines.push_back(geo_name);
      in.clear();
   }

	// If a linear function is given. 25.08.2011. WW
	if (getProcessDistributionType() == FiniteElement::FUNCTION)
	{
	  in.clear();
	  //dis_linear_f = new LinearFunctionData(*st_file);
	}

   if (this->getProcessDistributionType() == FiniteElement::LINEAR || this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
   {
      in >> nLBC;
      in.clear();
      for (int i = 0; i < nLBC; i++)
      {
         in.str(readNonBlankLineFromInputStream(*st_file));
         in >> abuff >> bbuff;
         in.clear();
         PointsHaveDistribedBC.push_back(abuff);
         DistribedBC.push_back(bbuff);
      }

      //      Read LINENODES AND VALUES......
      in.clear();
   }

   if (this->getProcessDistributionType() == FiniteElement::CRITICALDEPTH)
   {
      //KR critical_depth = true;
      in >> geo_node_value;
      in.clear();
      in.str(readNonBlankLineFromInputStream(*st_file));
      in >> rill_height;
      in.clear();
      //		dis_type = 6;
   }

   if (this->getProcessDistributionType() == FiniteElement::NORMALDEPTH)
   {
      dis_type_name = "NORMALDEPTH";
      in >> geo_node_value;
      in.clear();
      in.str(readNonBlankLineFromInputStream(*st_file));
      in >> normaldepth_slope >> rill_height;
      in.clear();
   }


   if (this->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
   {
      dis_type_name = "GREEN_AMPT";
      in >> geo_node_value;
      in.clear();
      in.str(readNonBlankLineFromInputStream(*st_file));
      in >> sorptivity >> constant >> rainfall >> moistureDeficit;
      in.clear();
   }
   // Soure terms are assign to element nodes directly. 23.02.2009. WW
   if(dis_type_name.find("DIRECT")!=std::string::npos)
   {
      dis_type_name = "DIRECT";
      in >> fname;
      fname = FilePath+fname;
      in.clear();
   }

   // Soure terms from precipitation are assign to element nodes directly.03.2010. WW
   if(dis_type_name.find("PRECIPITATION")!=std::string::npos)
   {
      dis_type_name = "PRECIPITATION";
      in >> fname;
      fname = FilePath+fname;
      std::ifstream ins(fname.c_str());
      if(!ins.good())
      {
         std::cout<<"Could not find file "<<fname<<std::endl;
         exit(0);
      }
      double timess;
      GIS_shape_head = new double[6];             // 07.06.2010. WW
      for(int i=0; i<6; i++)
      {
         getline(ins, aline);
         ss.str(aline);
         ss>> aline >> GIS_shape_head[i];
         ss.clear();

      }
      while(!ins.eof())
      {
         getline(ins, aline);
         if(aline.find("#STOP")!=std::string::npos)
            break;
         ss.str(aline);
         ss>> timess >> aline;
         precip_times.push_back(timess);
         precip_files.push_back(aline);

         ss.clear();
      }
      in.clear();
   }
}


/**************************************************************************
 FEMLib-Method:
 Task: for CSourceTerm::Read
 Programing:
 11/2007 JOD Implementation
 06/2010 TF modification of the signature, added geo_obj and unique_name
 **************************************************************************/
void CSourceTerm::ReadGeoType(std::ifstream *st_file,
//const GEOLIB::GEOObjects& geo_obj, 
const std::string& unique_name)
{
//   FileIO::GeoIO::readGeoInfo(this, *st_file, geo_name, geo_obj, unique_name);

   if (getProcessPrimaryVariable() == FiniteElement::EXCAVATION) //WW
   {
      std::stringstream strstr;
      strstr.str(readNonBlankLineFromInputStream(*st_file));
      //size_t tmp_geo_type;
      std::string sub_string;
      strstr >> sub_string >> _sub_dom_idx;
      strstr.clear();
   }
}


/**************************************************************************
 FEMLib-Method:
 Task: ST read function
 Programing:
 01/2004 OK Implementation
 06/2010 TF modification of the signature, added geo_obj and unique_name
 **************************************************************************/
bool STRead(const std::string &file_base_name,
//const GEOLIB::GEOObjects& geo_obj, 
const std::string& unique_name)
{
   char line[MAX_ZEILE];
   std::string line_string, st_file_name;
   std::ios::pos_type position;

   // File handling
   st_file_name = file_base_name + ST_FILE_EXTENSION;
   std::ifstream st_file(st_file_name.data(), std::ios::in);

   if (!st_file.good())
   {
      std::cout << "! Warning in STRead: No source terms !" << std::endl;
      return false;
   }

   // Keyword loop
   std::cout << "STRead ... " << std::flush;
   while (!st_file.eof())
   {
      st_file.getline(line, MAX_ZEILE);
      line_string = line;
                                                  //Code included to make dynamic memory for analytical solution
      if (line_string.find("#STOP") != std::string::npos)
      {
         size_t no_source_terms(st_vector.size());
         size_t no_an_sol = 0, number_of_terms = 0;
                                                  //Need to find the memory size limits for the anal. solution.
         for (size_t i = 0; i < no_source_terms; i++)
         {
            if (st_vector[i]->isAnalytical())
            {
               no_an_sol++;
               number_of_terms = std::max(st_vector[i]->getNumberOfTerms(), number_of_terms);
            }
         }
         if (no_an_sol > 0)
         {
            for (size_t i = 0; i < no_source_terms; i++)
            {
               st_vector[i]->setNumberOfAnalyticalSolutions (no_an_sol);
               st_vector[i]->setMaxNumberOfTerms (number_of_terms);
            }
         }

         std::cout << "done, read " << st_vector.size() << " source terms" << std::endl;
         return true;
      }
      //----------------------------------------------------------------------
                                                  // keyword found
      if (line_string.find("#SOURCE_TERM") != std::string::npos)
      {
         CSourceTerm *st(new CSourceTerm());
         std::ios::pos_type pos (st_file.tellg());
         position = st->Read(&st_file, unique_name);
         if (pos != position)
         {
            st_vector.push_back(st);
         }
         else
         {
            std::cerr << "WARNING: in STRead: could not read source term" << std::endl;
            delete st;
         }
         st_file.seekg(position, std::ios::beg);
      }                                           // keyword found
   }                                              // eof

   std::cout << "done, read " << st_vector.size() << " source terms" << std::endl;

   return true;
}


