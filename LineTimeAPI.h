#pragma once
#include "jsonSerializer.h"

typedef struct SimplifiedLineTime
{
	int frameid;
	int lineid;
	int actLine;
	double utctime;
	double avgtime;
}SimpleTime;

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
</LineTime>
*/
class LineTimeAPI : public JSONOBJECTSerializer
{
public:
	//LineTimeAPI::LineTimeAPI(){}
	LineTimeAPI::LineTimeAPI(unsigned char* _heifive)
	{
		IGiveMeFive(_heifive);
	}
	LineTimeAPI::LineTimeAPI(string _input, CEnumDomType _type)
	{
		objectname = "LineTimeinfo";
		SetType(_type);
		DeSerialize(_input.c_str());
		if (_type == asXml)
			Serialize("d:\\testfactory\\jsonhost", true);
	}

public:

public:
	LineTimeAPI(const LineTimeAPI& _information);
public:
	int Gtesttest(int _none)
	{
		return 5;
	}

private:
	/*
	{ "TIME"
				, "TIMEAGG",  "TIMESEG",
				"LINE", "ActLine", "UTCTIME"
				, "DATETIME", "LINESEG", "ReservedPara1",
				"ReservedPara2", "ReservedPara3", "ReservedPara4","ReservedPara5",
				"ReservedPara6", "ReservedPara7" };
	
	*/
	double TIMEAGG;
	double TIMESEG;
	int LINE;
	int ActLine;
	double UTCTIME;
	string DATETIME;
	int LINESEG;
	string ReservedPara1;
	string ReservedPara2;
	string ReservedPara3;
	string ReservedPara4;
	string ReservedPara5;
	string ReservedPara6;
	string ReservedPara7;
	string ReservedPara0;

public:
	double IGetTime() const
	{
		return UTCTIME;
	}

	int Line()
	{
		return LINE;
	}

	void Line(const int& _lineinfo)
	{
		LINE = _lineinfo;
	}

	void ACTLine(const int& _actLineinfo)
	{
		ActLine = _actLineinfo;
	}

	int ACTLine()
	{
		return ActLine;
	}

public:
	//唯一的默认构造函数，只有一个默认值
	LineTimeAPI(string inputstring) :TIMEAGG(0), TIMESEG(0), LINE(0), ActLine(0), UTCTIME(0), DATETIME("TIMESTRING"), LINESEG(0)
	{
		SetPropertys();
		ReservedPara0 = ReservedPara1 = ReservedPara2 = ReservedPara3 = ReservedPara4 = ReservedPara5 = ReservedPara6 = ReservedPara7 = inputstring;
	}

public:
	virtual bool IGiveMeFive(unsigned char * _content);

public:
	/*
	string ImageMode;
	int MainFrameNum;
	int ImageLine;
	string ImageTime;
	double GainLevel;
	int ImageCol;
	double HighTemperature;
	double LowTemperature;
	double FlatTemperature1;
	double FlatTemperature2;
	*/
protected:
	void NullInit();
	virtual void SetPropertys();
	virtual void SetUTCAPIByTimestring();
};
