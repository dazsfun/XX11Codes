//========================
//Log.cpp
#include "stdafx.h"
#include "Log.h"
Log::Log()
:m_bEnabled(true)
{
}
Log::~Log()
{
}
bool Log::Open(string sFileName)
{
	m_tOLogFile.open(sFileName.c_str(), ios_base::out | ios_base::app);
	if (!m_tOLogFile)
	{
		return false;
	}
	return true;
}
void Log::Close()
{
	if (m_tOLogFile.is_open())
	{
		m_tOLogFile.close();
	}
}
bool Log::CommonLogInit(string apistring)
{
	time_t tNowTime;
	time(&tNowTime);
	tm tLocalTime ; 
	localtime_s(&tLocalTime,&tNowTime);
	//得到日期的字符串
	string sDateStr = ValueToStr(tLocalTime.tm_year + 1900) + "-" +
		ValueToStr(tLocalTime.tm_mon + 1) + "-" +
		ValueToStr(tLocalTime.tm_mday);
	////Modified
	///return Open("IRS" + sDateStr + ".log");
	return Open(apistring + "sDateStr" + ".log");
}
void Log::Enable()
{
	m_bEnabled = true;
}
void Log::Disable()
{
	m_bEnabled = false;
}
//得到当前时间的字符串
string Log::GetTimeStr()
{
	time_t tNowTime;
	time(&tNowTime);
	tm tLocalTime;
	localtime_s(&tLocalTime,&tNowTime);
	//得到日期的字符串
	string strDateTime = "[" + ValueToStr(tLocalTime.tm_year + 1900) + "-" +
		ValueToStr(tLocalTime.tm_mon + 1) + "-" +
		ValueToStr(tLocalTime.tm_mday) + " " +
		ValueToStr(tLocalTime.tm_hour) + ":" +
		ValueToStr(tLocalTime.tm_min) + ":" +
		ValueToStr(tLocalTime.tm_sec) + "] ";
	return strDateTime;
}