
#include "Simulator.h"

#include "logog/include/logog.hpp"
#include "tclap/CmdLine.h"
#ifdef USE_LIS
#include "lis.h"
#endif

#include "BaseLib/CodingTools.h"
#include "ProcessLib/ProcessBuilder.h"
#include "FormatterCustom.h"
#include "SimulatorInfo.h"

namespace ogs6
{

////////////////////////////////////////////////////////////////////////////////
// Variables
////////////////////////////////////////////////////////////////////////////////
static FormatterCustom *custom_format;
static logog::Cout *logogCout;
static logog::LogFile *logog_file;

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


Simulator::Simulator(int argc, char* argv[])
{
	ogsInit(argc, argv);

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
	    cmd.parse( argc, argv );

	    // get parsed data
	    if (! output_arg.getValue().empty()) {
	    	if (!logog_file) delete logog_file;
	        logog_file = new logog::LogFile(output_arg.getValue().c_str());
	        logog_file->SetFormatter( *custom_format );
	    }

	    SimulatorInfo sim_info;
	    sim_info.output();

	    unsigned flag_list_modules (pcs_arg.getValue());
	    if (flag_list_modules!=0) {
	        ProcessLib::ProcessBuilder::getInstance()->output();
	    }

	    if (! input_arg.getValue().empty()) {
	    	std::string project_path = input_arg.getValue();

	    	//check files

	    }

    } catch (TCLAP::ArgException &e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}

}

Simulator::~Simulator()
{
	ogsExit();
}

int Simulator::execute()
{

    return 0;
}

} //end ogs6
