
#include "OutputBuilder.h"

#include "PVDOutput.h"
#include "TecplotOutput.h"

IOutput* OutputBuilder::create(const std::string &name)
{
    if (name.compare("PVD")==0) {
        return new PVDOutput();
    } else if (name.compare("TECPLOT")==0) {
        return new TecplotOutput();
    }
    return NULL;
}
