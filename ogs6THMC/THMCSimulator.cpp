
#include "THMCSimulator.h"

#ifdef USE_LIS
#include "lis.h"
#endif

#include "BaseLib/CodingTools.h"
#include "BaseLib/FileTools.h"
#include "SimulationInfo.h"
#include "SimulationProperties.h"
#include "Ogs6FemData.h"

namespace ogs6
{

////////////////////////////////////////////////////////////////////////////////
// Variables
////////////////////////////////////////////////////////////////////////////////
static logog::Cout *logogCout;
static logog::LogFile *logog_file;
static FormatterCustom *custom_format;
static bool isOgsInitCalled = false;
static bool isOgsExitCalled = false;
////////////////////////////////////////////////////////////////////////////////

void ogsInit(int argc, char* argv[])
{
	if (isOgsInitCalled) return;
	isOgsInitCalled = true;

    LOGOG_INITIALIZE();
    custom_format = new FormatterCustom();
    logogCout = new logog::Cout();
    logogCout->SetFormatter(*custom_format);
    logog_file = NULL;

#ifdef USE_LIS
    lis_initialize(&argc, &argv);
#endif
}

void ogsExit()
{
	if (isOgsExitCalled) return;
	isOgsExitCalled = true;

#ifdef USE_LIS
    lis_finalize();
#endif

	LOGOG_COUT << "exit ogs6." << std::endl;
    BaseLib::releaseObject(custom_format, logogCout, logog_file);
    LOGOG_SHUTDOWN();
}

THMCSimulator::THMCSimulator(int argc, char* argv[])
: _sim_info(NULL)
{
	try {
		// Command line parser
	    TCLAP::CmdLine cmd("ogs6", ' ', "0.1");
	    TCLAP::ValueArg<std::string> input_arg("i", "input", "input file", false, "", "string");
	    cmd.add( input_arg );
	    TCLAP::ValueArg<unsigned> n_cores_arg("p", "number-cores", "number of cores to use", false, 1, "number");
	    cmd.add( n_cores_arg );
	    TCLAP::ValueArg<std::string> output_arg("o", "output", "output file", false, "", "string");
	    cmd.add( output_arg );
	    TCLAP::ValueArg<unsigned> verbosity_arg("v", "verbose", "level of verbosity [0 very low information, 1 much information]", false, 0, "number");
	    cmd.add( verbosity_arg );
	    TCLAP::ValueArg<unsigned> pcs_arg("m", "modules", "list available modules [0 off, 1 on]", false, 1, "number");
	    cmd.add( pcs_arg );
	    cmd.parse( argc, argv ); // process can exit in this function

	    // initialize
		ogsInit(argc, argv);

	    // get parsed data
	    if (! output_arg.getValue().empty()) {
	    	if (!logog_file) delete logog_file;
	        logog_file = new logog::LogFile(output_arg.getValue().c_str());
	        logog_file->SetFormatter( *custom_format );
	    }

	    SimulationInfo::outputHeader();

	    unsigned flag_list_modules (pcs_arg.getValue());
	    if (flag_list_modules!=0) {
	    	ProcessBuilder::getInstance()->output();
	    }

	    if (! input_arg.getValue().empty()) {
	    	std::string proj_path = input_arg.getValue();
	    	if (checkInputFiles(proj_path)) {
			    _sim_info = new SimulationInfo(proj_path);
	    	} else {
	    		LOGOG_CERR << "Error: Cannot find file " << proj_path << std::endl;
	    	}
	    }

    } catch (TCLAP::ArgException &e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}

}

THMCSimulator::~THMCSimulator()
{
	BaseLib::releaseObject(_sim_info);
	ogsExit();
}

bool THMCSimulator::checkInputFiles(const std::string& proj_path)
{

	std::string proj_dir_path;
	std::string proj_name;

	std::string tmpFilename = proj_path;
	tmpFilename.append(".pcs");
	if(!BaseLib::IsFileExisting(tmpFilename))
	{
		LOGOG_CERR << " Error: Cannot find a PCS file - " << tmpFilename << std::endl;
		return 1;
	}

	return true;
}

int THMCSimulator::execute()
{
	if (!_sim_info) return 0;

	BaseLib::Options op;

	//-------------------------------------------------------------------------
	// Read files
	//-------------------------------------------------------------------------
	Ogs6FemData ogs6fem;
	const std::string proj_path = _sim_info->getProjectPath();

	// fem
	ogs5::Ogs5FemData ogs5femdata;
	ogs5femdata.read(proj_path);

	// ddc

	//-------------------------------------------------------------------------
	// Setup simulation
	//-------------------------------------------------------------------------
	// - create pcs
	Ogs5ToOgs6::convert(ogs5femdata, ogs6fem, op);
	// - ddc


	//-------------------------------------------------------------------------
	// Run simulation
	//-------------------------------------------------------------------------



    return 0;
}

} //end ogs6
