
#pragma once

#include <string>
#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"

#include "OutputTiming.h"

class IOutput
{
public:
    IOutput() : _output_timing(0) {};
    virtual ~IOutput()
    {
        BaseLib::releaseObject(_output_timing);
    };

    void setOutputPath(const std::string &dir, const std::string &base_name)
    {
        _output_dir_path = dir;
        _base_name = base_name;
    }

    std::string getOutputDir() const {return _output_dir_path;};
    std::string getOutputBaseName() const {return _base_name;};

    void addVariable(const std::string &var_name)
    {
        _list_var_name.push_back(var_name);
    };

    void addVariable(const std::vector<std::string> &vec_var_name)
    {
        _list_var_name.insert(_list_var_name.end(), vec_var_name.begin(), vec_var_name.end());
    };

    void setGeometry(const GeoLib::GeoObject* geo_obj) {_geo_obj = geo_obj;};
    const GeoLib::GeoObject* getGeometry() const {return _geo_obj;};

    void setOutputTiming(IOutputTiming* timing) {_output_timing = timing;};
    bool isActive(const NumLib::TimeStep &current_time) {return _output_timing->isActive(current_time);};

    virtual void write(const NumLib::TimeStep &current_time, std::map<std::string, NumLib::ITXFunction*> &data) = 0;

private:
    std::string _output_dir_path;
    std::string _base_name;
    IOutputTiming* _output_timing;
    std::vector<std::string> _list_var_name;
    const GeoLib::GeoObject* _geo_obj;
};

class PVDOutput : public IOutput
{
public:
    virtual void write(const NumLib::TimeStep &current_time, std::map<std::string, NumLib::ITXFunction*> &data)
    {

    }
};

class TecplotOutput : public IOutput
{
public:
    virtual void write(const NumLib::TimeStep &current_time, std::map<std::string, NumLib::ITXFunction*> &data)
    {

    }
};


class OutputBuilder
{
public:
    IOutput* create(const std::string &name)
    {
        if (name.compare("PVD")==0) {
            return new PVDOutput();
        } else if (name.compare("TECPLOT")==0) {
            return new TecplotOutput();
        }
        return NULL;
    }
};
