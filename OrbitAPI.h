#pragma once
#include "jsonSerializer.h"
#include <string>
#include <iostream>
#include <string>
using namespace std;


typedef class OrbitAPI : public JSONOBJECTSerializer
{
public:
	OrbitAPI(unsigned char* _information)
	{
		IGiveMeFive(_information);
	}
	OrbitAPI(string _input, CEnumDomType _inputType);
	virtual ~OrbitAPI();
public:

private:
	/*
	{ "GPSUTC"
	, "GPSDateTime", "GPSX",
	"GPSY", "GPSZ", "GPSXV"
	, "GPSYV", "GPSZV", "ReservedPara1",
	"ReservedPara2", "ReservedPara3", "ReservedPara4", "ReservedPara5",
	"ReservedPara6", "ReservedPara7" };
	*/
	double GPSUTC;
	string GPSDateTime;

public:
	double GPSX;
	double GPSY;
	double GPSZ;
	double GPSXV;
	double GPSYV;
	double GPSZV;

private:
	string ReservedPara1;
	string ReservedPara2;
	string ReservedPara3;
	string ReservedPara4;
	string ReservedPara5;
	string ReservedPara6;
	string ReservedPara7;
	double sideangle;


public:
	double IGetTime() const
	{
		return GPSUTC;
	}
public:
	OrbitAPI(string inputstring) : GPSUTC(0), GPSDateTime("TIMESTRING"), GPSX(0), GPSY(0), GPSZ(0),
		GPSXV(0), GPSYV(0), GPSZV(0)
	{
		SetPropertys();
		ReservedPara1 = ReservedPara2 = ReservedPara3 = ReservedPara4 = ReservedPara5 = ReservedPara6 = ReservedPara7 = inputstring;
	}


public:
	 void operator=( const OrbitAPI& hei)
	{
		 ClearPropertys();
		 SetPropertys();
		 GPSUTC = hei.GPSUTC;
		 GPSDateTime = hei.GPSDateTime; GPSX = hei.GPSX; GPSY = hei.GPSY; GPSZ = hei.GPSZ;
		 GPSXV = hei.GPSXV; GPSYV = hei.GPSYV; GPSZV = hei.GPSZV;
		 utcapi = hei.utcapi;
		 ReservedPara1 = hei.ReservedPara1;
		 ReservedPara2 = hei.ReservedPara2;
		 ReservedPara3 = hei.ReservedPara3;
		 ReservedPara4 = hei.ReservedPara4;
		 ReservedPara5 = hei.ReservedPara5;
		 ReservedPara6 = hei.ReservedPara6;
		 ReservedPara7 = hei.ReservedPara7;

	}
	
	OrbitAPI( const OrbitAPI& hei)
	{
		*this = hei;
	}


public:
	virtual bool IGiveMeFive(unsigned char* _content);
	//OrbitAPI(const OrbitAPI& _information);
protected:
	virtual void SetPropertys();
	void NullInit();
	///Set utc time by string 
	virtual void SetUTCAPIByTimestring();

}StrOrbitPoint;