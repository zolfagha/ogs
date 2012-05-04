
#pragma once

#include "Base/CodingTools.h"
#include "NumLib/TimeStepping/ITimeStepFunction.h"

namespace SolutionLib
{

class TimeSteppingProblem
{
public:
	TimeSteppingProblem() : _tim(0) {};

    virtual ~TimeSteppingProblem()
    {
        Base::releaseObject(_tim);
    }

	/// set  a time steppting function
	void setTimeSteppingFunction(NumLib::ITimeStepFunction &f)
	{
		_tim = f.clone();
	}

	/// get this time stepping function
	NumLib::ITimeStepFunction* getTimeSteppingFunction() const {return _tim;};

private:
	DISALLOW_COPY_AND_ASSIGN(TimeSteppingProblem);

    NumLib::ITimeStepFunction* _tim;
};

} //end
