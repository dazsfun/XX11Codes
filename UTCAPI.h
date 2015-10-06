#pragma once
#include "safestring.h"
#include <string>
#include <vector>
#include <map>
#include "DateTime.h"
using namespace std;

enum TimeStandard
{
	//XX–«        2000-01-01
	Enum_TimeStandard_Standard,
	///ZYŒ¿–«1∫≈  1998-01-01
	Enum_TimeStandard_ZY01,
	///ZYŒ¿–«»˝∫≈ 2012-01-09 2009-01-01
	Enum_TimeStandard_ZY03,
};

class UTCAPI 
{
public:
	UTCAPI();
	~UTCAPI()
	{
		_timelookupTable.clear();
		_dtimelookupTable.clear();
	}
	UTCAPI(string timestring, TimeStandard _standard);
	UTCAPI(double utcTime, TimeStandard _standard);

public:
	bool operator<=(const UTCAPI& utcinfo) const;
	bool operator>=(const UTCAPI& utcinfo) const;
	bool operator==(const UTCAPI& utcinfo) const;
	bool operator<(const UTCAPI& utcinfo) const;
	bool operator>(const UTCAPI& utcinfo) const;
	void operator=(const UTCAPI& utcinfo) ;


public:
	void ToUTC(string mystring,unsigned short = 0xC103);
	void ToTimeString(double myutc);
public:
	int year;
	int month;
	int day;
	int hour;
	int minute;
	double second;
	bool Is_Inited();
	double GetUTC() const
	{
		return utc;
	}

	string GetTimeString() const
	{
		return timestring;
	}
private:
	double utc;
	string timestring;
	void SetTimeStandard();
		
protected:
	bool _inited;
	ErrorType _errortype;
	double startutc;
	string startTimestring;
	TimeStandard mystandard;
	static int _YDays[13];

private:
	///Ω®¡¢≤È’“±Ì
	map<TimeStandard, string> _timelookupTable;
	map<TimeStandard, unsigned short> _dtimelookupTable;
	unsigned long DateTimetoUTC(unsigned char* datetime);
	void UTCtoDateTime(unsigned long ulUTC, unsigned char* datetime);
};