
#include <string>

#include "logog/include/logog.hpp"
#include "logog/include/formatter.hpp"
#include "tclap/CmdLine.h"

/**
 * new formatter for logog
 */
class FormatterCustom : public logog::FormatterGCC
{
    virtual TOPIC_FLAGS GetTopicFlags( const logog::Topic &topic )
    {
        return ( Formatter::GetTopicFlags( topic ) &
                 ~( TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG ));
    }
};

int main ( int argc, char *argv[] )
{
    LOGOG_INITIALIZE();

    TCLAP::CmdLine cmd("Simple matrix vector multiplication test", ' ', "0.1");

    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-m matrix".
    TCLAP::ValueArg<std::string> matrix_arg("m", "matrix", "input matrix file", true, "", "string");

    // Add the argument mesh_arg to the CmdLine object. The CmdLine object
    // uses this Arg to parse the command line.
    cmd.add( matrix_arg );

    TCLAP::ValueArg<unsigned> n_cores_arg("p", "number-cores", "number of cores to use", false, 1, "number");
    cmd.add( n_cores_arg );

    TCLAP::ValueArg<unsigned> n_mults_arg("n", "number-of-multiplications", "number of multiplications to perform", true, 10, "number");
    cmd.add( n_mults_arg );

    TCLAP::ValueArg<std::string> output_arg("o", "output", "output file", false, "", "string");
    cmd.add( output_arg );

    TCLAP::ValueArg<unsigned> verbosity_arg("v", "verbose", "level of verbosity [0 very low information, 1 much information]", false, 0, "string");
    cmd.add( verbosity_arg );

    cmd.parse( argc, argv );

    // read the number of multiplication to execute
    unsigned n_mults (n_mults_arg.getValue());
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



    delete custom_format;
    delete logogCout;
    delete logog_file;
    LOGOG_SHUTDOWN();

	return 0;
}
