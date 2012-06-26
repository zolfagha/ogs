/**************************************************************************
   FEMLib - Object: OUT
   Task:
   Programing:
   06/2004 OK Implementation
   last modified:
**************************************************************************/
//#include "Configure.h"
#include "rf_out_new.h"
// C++ STL
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <string>
using namespace std;

#include "makros.h"
#include "Output.h"

//// Base
//#include "StringTools.h"
//
//// FEM-Makros
//#include "files0.h"
//#include "makros.h"
//// GeoSys-GeoLib
//#include "GEOObjects.h"
//#include "files0.h"
//#include "geo_ply.h"
//#include "geo_sfc.h"
//// GeoSys-FEMLib
//#include "LegacyVtkInterface.h"
//#include "fem_ele_std.h"
//#include "mathlib.h"
//#include "rf_msp_new.h"
//#include "rf_pcs.h"
//#include "rf_pcs.h"
//#include "rf_random_walk.h"
//#include "rf_tim_new.h"
//// GeoSys-MSHLib
//#include "msh_lib.h"
//
//// FileIO/FEMIO
//#include "FEMIO/GeoIO.h"
//
//#include "problem.h"
//
//// Base
//#include "StringTools.h"
//
//extern size_t max_dim;                            //OK411 todo
//
//#ifdef CHEMAPP
//#include "eqlink.h"
//#endif
//#include "vtk.h"
//// MPI Parallel
//#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
//#include "par_ddc.h"
//#endif
//#ifdef SUPERCOMPUTER
//// kg44 this is usefull for io-buffering as endl flushes the buffer
//#define endl '\n'
//#define MY_IO_BUFSIZE 4096
//#endif
//#ifdef GEM_REACT
//#include "rf_REACT_GEM.h"
//#endif
//using MeshLib::CFEMesh;
////==========================================================================
vector<COutput*>out_vector;

/**************************************************************************
   FEMLib-Method:
   Task: OUT read function
   Programing:
   06/2004 OK Implementation
   08/2004 WW Remove the old files
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
   06/2006 WW Remove the old files by new way
   06/2010 TF reformated, restructured, signature changed, use new GEOLIB data structures
**************************************************************************/
bool OUTRead(const std::string& file_base_name,
//             const GEOLIB::GEOObjects& geo_obj,
             const std::string& unique_name)
{
	char line[MAX_ZEILE];
	std::string line_string;
	ios::pos_type position;
	bool output_version = false; // 02.2011. WW

	// File handling
	std::string out_file_name = file_base_name + OUT_FILE_EXTENSION;
	std::ifstream out_file(out_file_name.data(), ios::in);
	if (!out_file.good())
		return false;
	out_file.seekg(0L, ios::beg);

	// Keyword loop
	cout << "OUTRead" << endl;
	while (!out_file.eof())
	{
		out_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != string::npos)
			return true;

		COutput* out(new COutput(out_vector.size()));
		out->getFileBaseName() = file_base_name;
		// Give version in file name
		//15.01.2008. WW
		if (line_string.find("#VERSION") != string::npos)
			output_version = true;  // 02.2011. WW
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#OUTPUT") != string::npos)
		{
			position = out->Read(out_file, unique_name);

			//if(output_version) //// 02.2011. WW
			//{
			//	std::string VersionStr = OGS_VERSION; //02.2011 WX
			//	int curPos = 0;
			//	int pos = 0;
			//	while((pos = VersionStr.find("/",curPos)) != -1)
			//	{
			//		VersionStr.replace(pos, 1, "_");
			//		curPos = pos + 1;
			//	}
			//	out->getFileBaseName().append("(V");
			//	out->getFileBaseName().append(VersionStr);
			//	out->getFileBaseName().append(")");
			//}

			out_vector.push_back(out);

			//			char number_char[3]; //OK4709
			//			sprintf(number_char, "%i", (int) out_vector.size() - 1); //OK4709
			//			out->ID = number_char; //OK4709
			//			out->setID (out_vector.size() - 1);

			out_file.seekg(position, ios::beg);
		}                         // keyword found
	}                                     // eof
	return true;
}

/**************************************************************************
   FEMLib-Method: OUTWrite
   Task: master write function
   Programing:
   06/2004 OK Implementation
   last modification:
**************************************************************************/
void OUTWrite(string base_file_name)
{
	//========================================================================
	// File handling
	string out_file_name = base_file_name + OUT_FILE_EXTENSION;
	fstream out_file (out_file_name.data(),ios::trunc | ios::out);
	out_file.setf(ios::scientific,ios::floatfield);
	out_file.precision(12);
	if (!out_file.good())
		return;
	out_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	out_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//========================================================================
	out_file << "GeoSys-OUT: Output ------------------------------------------------\n";
	//========================================================================
	// OUT vector
	size_t out_vector_size (out_vector.size());
	for(size_t i = 0; i < out_vector_size; i++)
		out_vector[i]->Write(&out_file);
	out_file << "#STOP";
	out_file.close();
}

