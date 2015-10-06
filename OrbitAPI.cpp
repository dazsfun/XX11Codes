#include "stdafx.h"
#include "OrbitAPI.h"


void OrbitAPI::NullInit()
{
	SetPropertys();
	GPSUTC = GPSX = GPSY = GPSZ = GPSXV = GPSYV = GPSZV = 0;
	GPSDateTime = "";
	ReservedPara1 = ReservedPara2 = ReservedPara3 = ReservedPara4 = ReservedPara5
		= ReservedPara6 = ReservedPara7= "";
}

/*
<GpsData>
<GpsTime>2014-04-14T22:45:56.750369</GpsTime>
<PosX>3758458.750000</PosX>
<PosY>-12384.140000</PosY>
<PosZ>5717538.000000</PosZ>
<VelX>6198.200000</VelX>
<VelY>-2065.260000</VelY>
<VelZ>-4076.230000</VelZ>
</GpsData>
*/
void OrbitAPI::SetPropertys()
{
	ClearPropertys();
	///重构
	if (false)
	{
		/*  被废弃  这是最早版本的序列化尝试 因此被保留
		SetProperty("GpsTime", asString, &timestring);
		SetProperty("PosX", asDouble, &Xcomponent);
		SetProperty("PosY", asDouble, &Ycomponent);
		SetProperty("PosZ", asDouble, &Zcomponent);
		SetProperty("VelX", asDouble, &VecXcomponent);
		SetProperty("VelY", asDouble, &VecYcomponent);
		SetProperty("VelZ", asDouble, &VecZcomponent);
		*/
	}
	else
	{
		objectname = "轨道API";
		SetProperty(0 + 15, asDouble, &GPSUTC);
		SetProperty(1 + 15, asString, &GPSDateTime);
		SetProperty(2 + 15, asDouble, &GPSX);
		SetProperty(3 + 15, asDouble, &GPSY);
		SetProperty(4 + 15, asDouble, &GPSZ);;
		SetProperty(5 + 15, asDouble, &GPSXV);
		SetProperty(6 + 15, asDouble, &GPSYV);
		SetProperty(7 + 15, asDouble, &GPSZV);
		SetProperty(8 + 15, asString, &ReservedPara1);
		SetProperty(9 + 15, asString, &ReservedPara2);
		SetProperty(10 + 15, asString, &ReservedPara3);
		SetProperty(11 + 15, asString, &ReservedPara4);
		SetProperty(12 + 15, asString, &ReservedPara5);
		SetProperty(13 + 15, asString, &ReservedPara6);
		SetProperty(14 + 15, asString, &ReservedPara7);
	}
}

/*
OrbitAPI::OrbitAPI(const OrbitAPI& _information)
{
	SetPropertys();
	GPSUTC = _information.GPSUTC; 
	GPSDateTime = _information.GPSDateTime;
	GPSX = _information.GPSX;
	GPSY = _information.GPSY;
	GPSZ = _information.GPSZ;
	GPSXV = _information.GPSXV;
	GPSYV = _information.GPSYV;
	GPSZV = _information.GPSZV;
	ReservedPara1 = _information.ReservedPara1;
	ReservedPara2 = _information.ReservedPara2;
	ReservedPara3 = _information.ReservedPara3;
	ReservedPara4 = _information.ReservedPara4;
	ReservedPara5 = _information.ReservedPara5;
	ReservedPara6 = _information.ReservedPara6;
	ReservedPara7 = _information.ReservedPara7;
}
*/

OrbitAPI::OrbitAPI(string _input, CEnumDomType _type)
{
	SetPropertys();
	objectname = "OrbitInfo";
	SetType(_type);
	DeSerialize(_input.c_str());
	if (_type == asXml)
		Serialize("d:\\testfactory\\jsonhost", true);
}

void OrbitAPI::SetUTCAPIByTimestring()
{
	utcapi = UTCAPI(GPSDateTime, Enum_TimeStandard_Standard);
}

OrbitAPI::~OrbitAPI()
{

}

bool OrbitAPI::IGiveMeFive(unsigned char* _content)
{
	//unsigned char timesec[4];
	//memcpy(timesec, _content, 4);
	unsigned int time1 = safestring::FormatInt(_content, 4);
	//memcpy(&time1, _content, 4);
	unsigned int time2 = static_cast<unsigned int>(_content[4] << 16) + static_cast<unsigned int>(_content[5] << 8)
		+ _content[6];

	GPSUTC = time1 + time2 * 0.000001;

	UTCAPI api(GPSUTC, Enum_TimeStandard_ZY03);

	utcapi = api;

	GPSDateTime = api.GetTimeString();

	GPSX = safestring::FormatInt(_content + 8, 4) * 0.01;
	//memcpy(&GPSX, _content + 8, 4);
	GPSY = safestring::FormatInt(_content + 12, 4)*0.01;
	//memcpy(&GPSY, _content + 12, 4);
	GPSZ = safestring::FormatInt(_content + 16, 4)*0.01;
	//memcpy(&GPSZ, _content + 16, 4);
	GPSXV = safestring::FormatInt(_content + 20, 4)*0.01;
	//memcpy(&GPSXV, _content + 20, 4);
	GPSYV = safestring::FormatInt(_content + 24, 4)*0.01;
	//memcpy(&GPSYV, _content + 24, 4);
	GPSZV = safestring::FormatInt(_content + 28, 4)*0.01;
	//memcpy(&GPSZV, _content + 28, 4);

	unsigned int heitemp = 0;
	heitemp = safestring::FormatInt(_content + 32, 2);
	//memcpy(&heitemp, _content + 32, 2);
	sideangle = heitemp  * 0.0001;


	return true;
}