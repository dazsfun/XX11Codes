#pragma once
#include "RCProcess.h"
#include <string>
using namespace std;

#pragma comment(lib, "RCProcess.lib")
class RadioMetricProcess
{
public:
	void RCInterface(string RCInput,int RCMode)
	{
		if (RCMode == 11)
		{
			int overlap[8] = {18,18,18,18,18,18,18,18};
			X11RCProcessed(RCInput.c_str(),"MeanStat",overlap);
		}

		if (RCMode == 14)
		{
			int overlap[7] = {0,0,0,0,0,0,0};
			X6RCProcessed(RCInput.c_str(),"MeanStat",overlap);
		}

		if (RCMode == 15)
		{
			int overlap[7] = {18,18,0,0,0,0,0};
			O2CRCProcessed(RCInput.c_str(),"MeanStat",overlap);
		}
	}
	RadioMetricProcess(void);
	~RadioMetricProcess(void);
};

