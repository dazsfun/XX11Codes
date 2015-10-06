#include "stdafx.h"
#include "LineTimeAPI.h"

/*
<LineTime>
<ImageMode>IMGMODE</ImageMode>
<MainFrameNum>555</MainFrameNum>
<ImageLine>41</ImageLine>
<ImageTime>2014-04-14T22:45:55.178775</ImageTime>
<GainLevel>2.000000</GainLevel>
<ImageCol>18676</ImageCol>
<HighTemperature>297.269582</HighTemperature>
<LowTemperature>276.577981</LowTemperature>
<FlatTemperature1>0.000000</FlatTemperature1>
<FlatTemperature2>0.000000</FlatTemperature2>
</LineTime>*/
void LineTimeAPI::SetPropertys()
{
	if (false)
	{
		//下列的初始化方式由于不够自动化已经被废弃
		/*
		SetProperty("ImageMode", asString, &ImageMode);
		SetProperty("MainFrameNum", asInt, &MainFrameNum);
		SetProperty("ImageLine", asInt, &ImageLine);
		SetProperty("ImageTime", asString, &ImageTime);
		SetProperty("GainLevel", asDouble, &GainLevel);
		SetProperty("ImageCol", asInt, &ImageCol);
		SetProperty("HighTemperature", asDouble, &HighTemperature);
		SetProperty("LowTemperature", asDouble, &LowTemperature);
		SetProperty("FlatTemperature1", asDouble, &FlatTemperature1);
		SetProperty("FlatTemperature2", asDouble, &FlatTemperature2);
		*/
	}
	else
	{
		SetProperty(0 +30, asDouble, &TIMEAGG);
		SetProperty(1 +30, asDouble, &TIMESEG);
		SetProperty(2 +30, asInt, &LINE);
		SetProperty(3 +30, asInt, &ActLine);
		SetProperty(4 +30, asDouble, &UTCTIME);
		SetProperty(5 +30, asString, &DATETIME);
		SetProperty(6 +30, asInt, &LINESEG);
		SetProperty(7 +30, asString, &ReservedPara0);
		SetProperty(8 +30, asString, &ReservedPara1);
		SetProperty(9 +30, asString, &ReservedPara2);
		SetProperty(10 +30, asString, &ReservedPara3);
		SetProperty(11 +30, asString, &ReservedPara4);
		SetProperty(12 +30, asString, &ReservedPara5);
		SetProperty(13 +30, asString, &ReservedPara6);
		SetProperty(14 + 30, asString, &ReservedPara7);
	}
}

LineTimeAPI::LineTimeAPI(const LineTimeAPI& _information)
{
	SetPropertys();
	TIMEAGG = _information.TIMEAGG;
	TIMESEG = _information.TIMESEG;
	LINE = _information.LINE;
	ActLine = _information.ActLine;
	UTCTIME = _information.UTCTIME;
	LINESEG = _information.LINESEG;
	DATETIME = _information.DATETIME;
	ReservedPara0 = _information.ReservedPara0;
	ReservedPara1 = _information.ReservedPara1;
	ReservedPara2 = _information.ReservedPara2;
	ReservedPara3 = _information.ReservedPara3;
	ReservedPara4 = _information.ReservedPara4;
	ReservedPara5 = _information.ReservedPara5;
	ReservedPara6 = _information.ReservedPara6;
	ReservedPara7 = _information.ReservedPara7;
	utcapi = _information.utcapi;
}


void LineTimeAPI::SetUTCAPIByTimestring()
{
	utcapi = UTCAPI(UTCTIME, Enum_TimeStandard_Standard);
}

bool LineTimeAPI::IGiveMeFive(unsigned char* fivecontent)
{
	NullInit();

	unsigned char infobyte[6];
	
	memcpy(infobyte, fivecontent, sizeof(unsigned char)* 6);

	unsigned int part1 = static_cast<unsigned int>(infobyte[0] << 8) + infobyte[1];
	unsigned int part2 = static_cast<unsigned int>(infobyte[4] << 24) + static_cast<unsigned int>(infobyte[5] << 16)
		+ static_cast<unsigned int>(infobyte[2] << 8) + infobyte[3];

	TIMEAGG = part1 * 0.000025 + part2;

	UTCAPI api(TIMEAGG, Enum_TimeStandard_ZY03);

	utcapi = api;

	UTCTIME = TIMEAGG;

	DATETIME = api.GetTimeString();

	return true;
}

void LineTimeAPI::NullInit()
{
	TIMEAGG = TIMESEG = LINE = ActLine = UTCTIME = 0;
	DATETIME = ReservedPara0 = ReservedPara1 = ReservedPara2 = ReservedPara3 = ReservedPara4 = ReservedPara5 = ReservedPara6 = ReservedPara7 = "";
	LINESEG = 0;
}