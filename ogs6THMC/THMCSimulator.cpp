/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file THMCSimulator.cpp
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#include "THMCSimulator.h"

#include <iostream>

// external library
#include "logog/include/logog.hpp"
#include "tclap/CmdLine.h"
#ifdef USE_LIS
#include "lis.h"
#endif

// internal library
#include "BaseLib/CodingTools.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Options.h"
#include "BaseLib/OptionsXMLReader.h"
#include "BaseLib/FileTools.h"
#include "NumLib/TransientCoupling/TransientCouplingStructureBuilder.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "FemIO/ogs5/Ogs5FemIO.h"

// this module
#include "SimulationInfo.h"
#include "FormatterCustom.h"
#include "Ogs6FemData.h"
#include "Ogs5ToOgs6.h"
#include "MyConvergenceCheckerFactory.h"


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
: _sim_info(NULL), _cpl_system(NULL)
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
    BaseLib::releaseObject(_sim_info, _cpl_system);
    ogsExit();
}

bool THMCSimulator::checkInputFiles(const std::string& proj_path)
{
    // meanwhile OGS5 files are default
    if(!BaseLib::IsFileExisting(proj_path+".pcs"))
    {
        LOGOG_CERR << " Error: Cannot find a PCS file - " << proj_path << std::endl;
        return 1;
    }
    if(!BaseLib::IsFileExisting(proj_path+ ".prop"))
    {
        LOGOG_CERR << " Error: Cannot find a PCS file - " << proj_path << std::endl;
        return 1;
    }

    return true;
}

int THMCSimulator::execute()
{
    if (!_sim_info) return 0;

    BaseLib::Options op;
    Ogs6FemData* ogs6fem = Ogs6FemData::getInstance();
    const std::string proj_path = _sim_info->getProjectPath();

    //-------------------------------------------------------------------------
    // Read files
    //-------------------------------------------------------------------------
    LOGOG_CERR << "Reading input files..." << std::endl;
    // ogs5fem
    ogs5::Ogs5FemData ogs5femdata;
    ogs5femdata.read(proj_path);
    Ogs5ToOgs6::convert(ogs5femdata, *ogs6fem, op);

    // coupling
    BaseLib::addXMLtoOptions(proj_path+".prop", op);

    // ddc

    //-------------------------------------------------------------------------
    // Setup simulation
    //-------------------------------------------------------------------------
    LOGOG_CERR << "Setting up simulation..." << std::endl;

    // construct coupling system
    MyConvergenceCheckerFactory checkFac;
    NumLib::TransientCoulplingStrucutreBuilder cpl_builder;
    if (_cpl_system!=NULL) delete _cpl_system;
    _cpl_system = cpl_builder.build(&op, *GeoProcessBuilder::getInstance(), checkFac);
    if (!_cpl_system->check()) {
        LOGOG_CERR << " Error while checking coupled system " << std::endl;
        return 0;
    }

    // list up monolithic processes
    std::vector<NumLib::AbstractTransientMonolithicSystem*> &list_mono_system = cpl_builder.getListOfMonolithicSystem();
    std::vector<std::string> &list_mono_system_name = cpl_builder.getListOfMonolithicSystemName();
    for (size_t i=0; i<list_mono_system.size(); i++) {
        std::string &pcs_name = list_mono_system_name[i];
        ProcessLib::Process* pcs = list_mono_system[i];
        ogs6fem->list_pcs.insert(pcs_name, pcs);
        const BaseLib::Options* opPCS = op.getSubGroup("ProcessData")->getSubGroup(pcs_name);
        if (opPCS!=NULL) {
            pcs->initialize(*opPCS);
        } else {
            LOGOG_CERR << " Error: Cannot find Configuration for Process - " << pcs_name << std::endl;
            return 0;
        }
    }


    NumLib::TimeSteppingController timestepping;
    timestepping.addTransientSystem(*_cpl_system);

    double t_start = std::numeric_limits<double>::max();
    double t_end = std::numeric_limits<double>::min();

    for (size_t i=0; i<ogs6fem->list_tim.size(); i++) {
        t_start = std::min(t_start, ogs6fem->list_tim[i]->getBeginning());
        t_end = std::max(t_start, ogs6fem->list_tim[i]->getEnd());
    }

    //-------------------------------------------------------------------------
    // Run simulation
    //-------------------------------------------------------------------------
    LOGOG_CERR << "Start simulation...  start=" << t_start << ", end=" << t_end << std::endl;
    timestepping.setBeginning(t_start);
    timestepping.solve(t_end);

    LOGOG_CERR << "Finish simulation..." << std::endl;


    return 0;
}

} //end ogs6
