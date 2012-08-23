#include "chemcomp.h"

namespace ogsChem
{

ChemComp::ChemComp(void):
    _name("NaN"), _composition("NaN"), _mole_mass(0.0),
    _charge(0.0), _G0(0.0)
{
    // set pointers to zero
    _mPhase    = 0;
    _mReaction = 0;
}

ChemComp::~ChemComp(void)
{
}

double ChemComp::get_mole_mass(void) const
{
return _mole_mass;
}

double ChemComp::get_mole_volume(void) const
{
return _mole_volume;
}

double ChemComp::get_charge(void) const
{
return _charge;
}

} // end of namespace
