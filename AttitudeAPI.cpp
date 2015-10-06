#include "stdafx.h"
#include "AttitudeAPI.h"

/*
<AttType>Quater</AttType>
<AttTime>2014-04-14T22:45:59.261249</AttTime>
<Q1>0.829922</Q1>
<Q2>0.426848</Q2>
<Q3>-0.147407</Q3>
*/
void AttitudeAPI::SetPropertys()
{
	ClearPropertys();
	///序列化反序列化代码重构
	if (false)
	{

	}
	else
	{
		objectname = "姿态API";
		SetProperty(0 + 0, asDouble, &AttitudeUTC);
		SetProperty(1 + 0, asString, &AttitudeDateTime);
		SetProperty(2 + 0, asInt, &AttitudeType);
		SetProperty(3 + 0, asDouble, &AttitudeROLL);
		SetProperty(4 + 0, asDouble, &AttitudePITCH);
		SetProperty(5 + 0, asDouble, &AttitudeYAW);
		SetProperty(6 + 0, asDouble, &AttitudeQ0);
		SetProperty(7 + 0, asDouble, &AttitudeQ1);
		SetProperty(8 + 0, asDouble, &AttitudeQ2);
		SetProperty(9 + 0, asDouble, &AttitudeQ3);
		SetProperty(10 + 0, asDouble, &AttitudeROLLVec);
		SetProperty(11 + 0, asDouble, &AttitudePitchVec);
		SetProperty(12 + 0, asDouble, &AttitdueYAWVec);
		SetProperty(13 + 0, asString, &ReservedPara1);
		SetProperty(14 + 0, asString, &ReservedPara2);

	}
}

AttitudeAPI::AttitudeAPI(const AttitudeAPI& _getinfo)
{
	ClearPropertys();
	SetPropertys();
	utcapi = _getinfo.utcapi;
	this->AttitdueYAWVec = _getinfo.AttitdueYAWVec;
	this->AttitudePitchVec = _getinfo.AttitudePitchVec;
	this->AttitudeROLLVec = _getinfo.AttitudeROLLVec;
	this->AttitudeROLL = _getinfo.AttitudeROLL;
	this->AttitudePITCH = _getinfo.AttitudePITCH;
	this->AttitudeYAW = _getinfo.AttitudeYAW;
	this->AttitudeQ0 = _getinfo.AttitudeQ0;
	this->AttitudeQ1 = _getinfo.AttitudeQ1;
	this->AttitudeQ2 = _getinfo.AttitudeQ2;
	this->AttitudeQ3 = _getinfo.AttitudeQ3;
	this->AttitudeUTC = _getinfo.AttitudeUTC;
	this->AttitudeDateTime = _getinfo.AttitudeDateTime;
	this->AttitudeType = _getinfo.AttitudeType;
	this->ReservedPara1 = _getinfo.ReservedPara1;
	this->ReservedPara2 = _getinfo.ReservedPara2;
}

void AttitudeAPI::SetUTCAPIByTimestring()
{
	utcapi = UTCAPI(AttitudeUTC, Enum_TimeStandard_Standard);
}

void AttitudeAPI::NullInit()
{
	SetPropertys();
	AttitudeUTC = AttitudeType = AttitudeROLL = AttitudePITCH = AttitudeYAW
		= AttitudeQ0 = AttitudeQ1 = AttitudeQ2 = AttitudeQ3 = 0;
	AttitudeROLLVec = AttitudePitchVec = AttitdueYAWVec;
	ReservedPara1 = ReservedPara2 = "";
}

#include <math.h>
bool AttitudeAPI::IGiveMeFive(unsigned char* _content)
{
	unsigned char _attitudemode = _content[1];
	unsigned char _observemode = _content[0];

	if (_attitudemode == 0xC4)
	{
		unsigned char  reader[3][4];
		for (int i = 0; i < 3; i++)
		{
			reader[i][0] = *(_content + 10 + i * 4 + 2);
			reader[i][1] = *(_content + 10 + i * 4 + 3);
			reader[i][2] = *(_content + 10 + i * 4 + 0);
			reader[i][3] = *(_content + 10 + i * 4 + 1);
		}

		AttitudeQ1 = safestring::FormatDouble(reader[0], 4);
		//memcpy(&AttitudeQ1, reader[0], sizeof(char)* 4);
		AttitudeQ2 = safestring::FormatDouble(reader[1], 4);
		//memcpy(&AttitudeQ2, reader[1], sizeof(char)* 4);
		AttitudeQ3 = safestring::FormatDouble(reader[2], 4);
		//memcpy(&AttitudeQ3, reader[2], sizeof(char)* 4);

		AttitudeQ0 = sqrt(1 - pow(AttitudeQ1, 2.0) - pow(AttitudeQ2, 2) - pow(AttitudeQ3, 2.0));
	}

	else if (_attitudemode == 0x53)
	{
		int  shortreader[6];
		for (int i = 0; i < 6; i++)
		{
			shortreader[i] = safestring::FormatDouble
				(_content + 8 + i * 2, 2);
			//memcpy(&shortreader[i], _content + 8 + i * 2, 2);
		}

		AttitudeROLL = shortreader[0] * 0.0001;
		AttitudePITCH = shortreader[1] * 0.0001;
		AttitudeYAW = shortreader[2] * 0.0001;

		AttitudeROLLVec = shortreader[3] * 5.5 * 0.0000001;
		AttitudePitchVec = shortreader[4] * 5.5 * 0.0000001;
		AttitdueYAWVec = shortreader[5] * 5.5 * 0.0000001;
	}

	else
	{
		return false;
	}

	unsigned char timeinfo[6];
	memcpy(timeinfo, _content + 2, 6);
	unsigned int time1 = static_cast<unsigned int>(timeinfo[0] << 8) + timeinfo[1];
	unsigned int time2 = static_cast<unsigned int>(timeinfo[4] << 24) + static_cast<unsigned int>(timeinfo[5] << 16)
		+ static_cast<unsigned int>(timeinfo[2] << 8) + timeinfo[3];
	double currenttime = time1 * 0.000025 + time2;

	short  timeoffset = safestring::FormatShort(_content + 8, 2);
	//memcpy(&timeoffset, _content + 8,2);
	if (abs(timeoffset) >= 800)
	{
		cout << "warning:数管校时超过指定量程，相应值将不被记录" << "其值为:" << timeoffset << endl;
		IAddLog("Warning: 数管校时超过指定量程，相应值将不被 记录,其值为:%d", timeoffset);
		AttitudeUTC = currenttime;
	}
	else
		AttitudeUTC = currenttime + timeoffset * 0.001;
	UTCAPI api(AttitudeUTC, Enum_TimeStandard_ZY03);
	AttitudeDateTime = api.GetTimeString();
	AttitudeType = _attitudemode;
	utcapi = api;
	return true;
}