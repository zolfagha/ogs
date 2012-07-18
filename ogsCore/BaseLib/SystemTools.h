
#pragma once


namespace BaseLib
{

inline bool IsLittleEndian()
{
    int x = 0x00000001;
    if (*(char*)&x)
        return true;              //am little
    else
        return false;             //am big
}

}
