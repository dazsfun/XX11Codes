// Log.h
#pragma  once
#include <stdarg.h>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
using namespace std;


typedef enum ErrorType
{
	Process_Success,
	File_Invalid_Existence,
	File_Invalid_Content,
	File_Invalid_Existence_Fred,
	Fred_Cannot_Find_Header,
	Fred_Unexpected_End,
	Format_Invalid_TimeString,
	Range_Out_UTCTYPE,
}ReactionType;

class Log
{
public:
	Log();
	~Log();
	bool Open(string sFileName);
	void Close();
	bool CommonLogInit(string apistring = "logrecorder"); //打开默认的log 文件
	void Enable();
	void Disable();
	string GetTimeStr();
	template <typename T> void LogOut(const T& value)
	{
		if (m_bEnabled)
		{
			m_tOLogFile << value;
		}
	}
	template <typename T> void LogOutLn(const T& value)
	{
		if (m_bEnabled)
		{
			m_tOLogFile << value << endl;
		}
	}
	void LogOutLn()
	{
		if (m_bEnabled)
		{
			m_tOLogFile << endl;
		}
	}
	template <typename T> Log& operator<<(const T& value)
	{
		if (m_bEnabled)
		{
			m_tOLogFile << value;
		}
		return (*this);
	}
	Log& operator<<(ostream& (*_Pfn)(ostream&))
	{
		if (m_bEnabled)
		{
			(*_Pfn)(m_tOLogFile);
		}
		return (*this);
	}
private:
	template<typename T> string ValueToStr(T value)
	{
		ostringstream ost;
		ost << value;
		return ost.str();
	}
private:
	ofstream m_tOLogFile;
	bool m_bEnabled;
};

class ConfigAPI 
{
public:
	ConfigAPI(string typestring,string id, string cmdstring)
	{
		Log commandlog;
		commandlog.CommonLogInit(typestring);
		commandlog << id.c_str() << "---" << commandlog.GetTimeStr() << endl;
		commandlog << cmdstring.c_str() << endl;
		commandlog.Close();
	}
};

///日志类，承担程序的日志记录工作
class LogAPI
{
public:
	LogAPI(){}
	virtual ~LogAPI(){}

public:
	///使用方法，与 printf使用方式一样
	void IAddLog(const char* sz, ...)
	{

		char szData[1024];
		va_list args;
		va_start(args, sz);
		_vsnprintf_s(szData, 1024 - 2, sz, args);
		va_end(args);

		Log mainLog;
		mainLog.CommonLogInit();
		mainLog << mainLog.GetTimeStr() << szData << endl;
		mainLog.Close();
	}
};
