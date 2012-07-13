
#pragma once

#include <string>

namespace ProcessLib
{

/**
 *
 */
class ProcessInfo
{
public:
	~ProcessInfo() {};

    const char* getName() const { return _name.c_str(); }

private:
    std::string _name;
};

}
