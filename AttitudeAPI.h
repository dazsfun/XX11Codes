#pragma once
#include "jsonSerializer.h"
using namespace std;


/*
<AttType>Quater</AttType>
<AttTime>2014-04-14T22:45:59.261249</AttTime>
<Q1>0.829922</Q1>
<Q2>0.426848</Q2>
<Q3>-0.147407</Q3>
*/
typedef class AttitudeAPI : public JSONOBJECTSerializer
{
public:
	AttitudeAPI(string _input, CEnumDomType _type)
	{
		objectname = "AttitudeInfo";
		SetType(_type);
		DeSerialize(_input.c_str());
		if (_type == asXml)
			Serialize("d:\\testfactory\\jsonhost", true);
	}

	AttitudeAPI(unsigned char* _information)
	{
		IGiveMeFive(_information);
	}

public:
	AttitudeAPI(const AttitudeAPI& _getinfo);

private:
	/*
	{ "AttitudeUTC"
	, "AttitudeDateTime", "AttitudeType", "AttitudeROLL",
	"AttitudePITCH", "AttitudeYAW", "AttitudeQ0"
	, "AttitudeQ1", "AttitudeQ2", "AttitudeQ3",
	"AttitudeRollVec", "AttitudePitchVec", "AttitudeYAWVec",
	"ReservedPara1", "ReservedPara2" };
	*/
	double  AttitudeUTC;
	string AttitudeDateTime;
	int AttitudeType;

public:
	double AttitudeROLL;
	double AttitudePITCH;
	double AttitudeYAW;
	double AttitudeQ0;
	double AttitudeQ1;
	double AttitudeQ2;
	double AttitudeQ3;
	double AttitudeROLLVec;
	double AttitudePitchVec;
	double AttitdueYAWVec;

private:
	string ReservedPara1;
	string ReservedPara2;

public:
	double IGetTime() const
	{
		return AttitudeUTC;
	}

public:
	AttitudeAPI(string inputstring) :AttitudeUTC(0), AttitudeDateTime("TIMESTRING"),
		AttitudeType(0), AttitudeROLL(0), AttitudePITCH(0), AttitudeYAW(0),
		AttitudeQ0(0), AttitudeQ1(0), AttitudeQ2(0), AttitudeQ3(0), AttitudeROLLVec(0), AttitudePitchVec(0),
		AttitdueYAWVec(0)
	{
		SetPropertys();
		ReservedPara1 = ReservedPara2 = inputstring;
	}
protected:
	virtual void SetPropertys();
	void NullInit();
	virtual void SetUTCAPIByTimestring();

public:
	virtual bool IGiveMeFive(unsigned char* _content);
}StrAttPoint;
