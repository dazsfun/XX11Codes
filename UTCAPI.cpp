#include "stdafx.h"
#include "UTCAPI.h"

#define BCDtoHEX(x)	((x)>>4)*10 + ((x)&0x0f)
#define HEXtoBCD(x) ((x)/10)*16+((x)%10)

int UTCAPI::_YDays[13] = { 0, 0, 31, 59, 90, 120, 151, 181, 212, 243,273 ,304, 334 };

bool UTCAPI::Is_Inited()
{
	return _inited;
}

bool UTCAPI::operator<=(const UTCAPI& utcinfo) const
{
	if (this->utc > utcinfo.GetUTC())
		return false;
	else return true;
}

bool UTCAPI::operator>=(const UTCAPI& utcinfo) const
{
	if (this->utc < utcinfo.GetUTC())
		return false;
	else return true;
}

bool UTCAPI::operator<(const UTCAPI& utcinfo) const
{
	return (this->utc < utcinfo.GetUTC());
}

bool UTCAPI::operator>(const UTCAPI& utcinfo) const
{
	return (this->utc > utcinfo.GetUTC());
}

bool UTCAPI::operator==(const UTCAPI& utcinfo) const
{
	return (this->utc == utcinfo.GetUTC());
}

UTCAPI::UTCAPI()
{}

UTCAPI::UTCAPI(double utctime, TimeStandard _standrd)
{
	SetTimeStandard();
	_inited = false;
	_errortype = Process_Success;
	map<TimeStandard, unsigned short>::iterator l_it;
	l_it = _dtimelookupTable.find(_standrd);
	if (l_it == _dtimelookupTable.end())
	{
		_errortype = Range_Out_UTCTYPE;
	}
	else
	{
		unsigned short mytimestring = l_it->second;
		CDateTime date_1(mytimestring, utctime);
		//startutc = utc;
		utc = utctime;
		//utc += 1;
		char tempreader[1024];
		sprintf_s(tempreader, 1024, "%04d-%02d-%02dT%02d-%02d-%02.8lf", date_1.year, date_1.month, date_1.day, date_1.hour, date_1.minute, date_1.second);
		timestring = tempreader;
		year = date_1.year;
		month = date_1.month;
		day = date_1.day;
		hour = date_1.hour;
		minute = date_1.minute;
		second = date_1.second;
		_inited = true;
	}

}

void UTCAPI::operator=(const UTCAPI& hei)
{
	this->year = hei.year;
	this->month = hei.month;
	this->day = hei.day;
	this->hour = hei.hour;
	this->minute = hei.minute;
	this->second = hei.second;
	this->utc = hei.GetUTC();
	this->timestring = hei.GetTimeString();
	this->_inited = true;
	//SetTimeStandard();
}

UTCAPI::UTCAPI(string _timestring, TimeStandard _standard)
{
	timestring = _timestring;
	SetTimeStandard();
	_inited = false;
	_errortype = Process_Success;
	map<TimeStandard, string>::iterator  l_it;
	l_it = _timelookupTable.find(_standard);
	if (l_it == _timelookupTable.end())
	{
		_errortype = Range_Out_UTCTYPE;
	}
	else
	{
		string mytimestring = l_it->second;
		ToUTC(_timestring);
		if (_errortype != Process_Success)
		{
			_inited = false;
			year = month = day = hour = minute = second = 0;
		}
		else
		{
			ToUTC(mytimestring);
			startutc = utc;
			ToUTC(_timestring);
			utc -= startutc;
			_inited = true;
		}
	}
}

/// year--year month--month day--day hour--hour minute-minute second--second utc--utc timestring--timestring


///算法函数，根据年月日时分秒，计算该时间信息
unsigned long UTCAPI::DateTimetoUTC(unsigned char* datetime)
{
	unsigned int _year;
	unsigned int _month, _day, _hour, _minute, _second;
	unsigned int _leaps = 0;
	unsigned long _days = 0;
	unsigned long _secs = 0;

	_year = BCDtoHEX(datetime[0]);
	if (_year <= 70) _year += 2000; else _year += 1900;
	_month = BCDtoHEX(datetime[1]);
	_day = BCDtoHEX(datetime[2]);
	_hour = BCDtoHEX(datetime[3]);
	_minute = BCDtoHEX(datetime[4]);
	_second = BCDtoHEX(datetime[5]);
	if (_year < 1970 || _year > 2099) return 0;
	_year -= 1970;
	_leaps = (_year + 2) / 4;
	if (!((_year + 2) & 3) && (_month < 3))
		--_leaps;

	_days = _year * 365L + _leaps + _YDays[_month] + _day - 1;
	_secs = _days * 86400L + _hour * 3600L + _minute * 60L + _second;
	return _secs;
}

void UTCAPI::UTCtoDateTime(unsigned long ulUTC, unsigned char* datetime)
{
	unsigned long ulHMS;
	unsigned int uiYMD, uiFYMD;
	unsigned int _year;
	unsigned char _month, _day, _hour, _minute, _second, _leap;

	uiYMD = ulUTC / 86400L;   //// 天数
	ulHMS = ulUTC % 86400L; ///   剩余总秒数 

	_second = (unsigned char)(ulHMS % 60L);
	ulHMS = ulHMS / 60L;
	_minute = (unsigned char)(ulHMS % 60L);
	ulHMS = ulHMS / 60L;
	_hour = (unsigned char)(ulHMS % 24L);

	_year = (uiYMD / 1461L) * 4;
	uiFYMD = uiYMD % 1461L;

	_leap = 0;
	if (uiFYMD > 1095)
	{
		uiFYMD -= 1096; _year += 3; _leap = 0;
	}
	else if (uiFYMD > 729)
	{
		uiFYMD -= 730; _year += 2; _leap = 1;
	}
	else if (uiFYMD > 364)
	{
		uiFYMD -= 365; _year += 1; _leap = 0;
	}

	_year += 1970;

	///月数判断
	for (int i = 0; i < 12; i++)
	{
		if (i == 11)
		{
			_month = 1;
			_day = uiFYMD + 1;
			break;
		}
		if (uiFYMD >(_YDays[12 - i] - 1 + _leap))
		{
			_day = uiFYMD - _YDays[12 - i] - _leap + 1;
			_month = 12 - i;
			break;
		}
	}

	
	datetime[0] = HEXtoBCD(_year % 100);
	datetime[1] = HEXtoBCD(_month);
	datetime[2] = HEXtoBCD(_day);
	datetime[3] = HEXtoBCD(_hour);
	datetime[4] = HEXtoBCD(_minute);
	datetime[5] = HEXtoBCD(_second);
	

}

void UTCAPI::ToTimeString(double myutc)
{
	unsigned long temputc = static_cast<unsigned long>(myutc);

	unsigned char * _temptimecharstring = new unsigned char[6];

	UTCtoDateTime(temputc, _temptimecharstring);
	year = (_temptimecharstring[0] > 70) ? _temptimecharstring[0] + 1900:_temptimecharstring[0] + 2000;
	month = _temptimecharstring[1];
	day = _temptimecharstring[2];
	hour = _temptimecharstring[3];
	minute = _temptimecharstring[4];
	second = _temptimecharstring[5] + (myutc - temputc);

	char tempreader[1024];
	sprintf_s(tempreader, 1024, "%04d-%02d-%02d-%02d-%02d-%02.8lf", year, month, day, hour, minute, second);
	timestring = tempreader;
	delete[] _temptimecharstring;
}

void UTCAPI::ToUTC(string mytimestring, unsigned short markid)
{
	if (safestring::Is_TimeFormated(mytimestring, year, month, day, hour, minute, second) != Process_Success)
	{
		_errortype = Format_Invalid_TimeString;
		return;
	}

	CDateTime date_1(markid, year, month, day, hour, minute, second);
	utc = date_1.countDate2Second();
	/*
	unsigned char* _tempin = new unsigned char[6];
	_tempin[0] = year % 100;
	_tempin[1] = month;
	_tempin[2] = day;
	_tempin[3] = hour;
	_tempin[4] = minute;
	_tempin[5] = static_cast<int>(second);
	
	utc = DateTimetoUTC(_tempin) + (second)-static_cast<int>(second);
	
	delete[] _tempin;
	*/
}

void UTCAPI::SetTimeStandard()
{
	_timelookupTable[Enum_TimeStandard_Standard] = "2009-01-01T00:00:00.000000";
	_dtimelookupTable[Enum_TimeStandard_Standard] = 0xC103;
	_dtimelookupTable[Enum_TimeStandard_ZY01] = 0xC101;
	_dtimelookupTable[Enum_TimeStandard_ZY03] = 0xC102;
}