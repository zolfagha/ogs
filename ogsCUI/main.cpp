
#include <string>

#include "logog/include/logog.hpp"
#include "tclap/CmdLine.h"

#include "ProcessLib/ProcessBuilder.h"
#include "ApplicationLib/ProcessList.h"
#include "FormatterCustom.h"


int main ( int argc, char *argv[] )
{

    LOGOG_INITIALIZE();

    TCLAP::CmdLine cmd("ogs6-nw", ' ', "0.1");

    TCLAP::ValueArg<std::string> matrix_arg("i", "input", "input file", false, "", "string");
    cmd.add( matrix_arg );

    TCLAP::ValueArg<unsigned> n_cores_arg("p", "number-cores", "number of cores to use", false, 1, "number");
    cmd.add( n_cores_arg );

    TCLAP::ValueArg<std::string> output_arg("o", "output", "output file", false, "", "string");
    cmd.add( output_arg );

    TCLAP::ValueArg<unsigned> verbosity_arg("v", "verbose", "level of verbosity [0 very low information, 1 much information]", false, 0, "number");
    cmd.add( verbosity_arg );

    TCLAP::ValueArg<unsigned> pcs_arg("m", "modules", "list available modules [0 off, 1 on]", false, 1, "number");
    cmd.add( pcs_arg );

    cmd.parse( argc, argv );

    // read the number of multiplication to execute
    std::string fname_mat (matrix_arg.getValue());

    FormatterCustom *custom_format (new FormatterCustom);
    logog::Cout *logogCout(new logog::Cout);
    logogCout->SetFormatter(*custom_format);

    logog::LogFile *logog_file(NULL);
    if (! output_arg.getValue().empty()) {
        logog_file = new logog::LogFile(output_arg.getValue().c_str());
        logog_file->SetFormatter( *custom_format );
    }

    // read number of threads
    unsigned n_threads (n_cores_arg.getValue());

    unsigned flag_list_modules (pcs_arg.getValue());
    if (flag_list_modules!=0) {
    	ProcessLib::ProcessBuilder::getInstance()->output();
    }


    delete custom_format;
    delete logogCout;
    delete logog_file;
    LOGOG_SHUTDOWN();

	return 0;
}
