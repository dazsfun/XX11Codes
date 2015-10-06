// EquinoxBasedT2C.cpp: implementation of the CEquinoxBasedT2C class.
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/******************************************************************
/** 文件名：EquinoxBasedT2C.cpp
/** Copyright (c) 2003-2004 武汉武大吉奥信息技术有限公司
/** 创建人：张过
/** 日  期：2004/8/18
/** 修改人：
/** 日  期：
/** 描  述：
/** 算法说明：
/** 版  本：1.0
/**-----------------------------------------------------------------------------
/*******************************************************************/

#include "stdafx.h"

#include "EquinoxBasedT2C.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


#define Tref_year 2009
#define Tref_month 1
#define Tref_day 1
#define Tref_hour 0
#define Tref_minute 0
#define Tref_second 0.0





CEquinoxBasedT2C::CEquinoxBasedT2C()
{

}

CEquinoxBasedT2C::~CEquinoxBasedT2C()
{

}


BOOL CEquinoxBasedT2C::Initdir(CString path, CString ipath)
{
	CString dir = path + ipath;
	strncpy_s(m_dirconfig,LPCTSTR(dir),sizeof(m_dirconfig));
	return true;

}

void CEquinoxBasedT2C::Inittime(int year, int month, int day, int hour, int minute, double second)
{
	MJD_SOD(&year, &month, &day, &hour, &minute, &second, &mjd, &sod);
	UTC2TDT(&mjd, &sod,m_dirconfig, strlen(m_dirconfig));	
}



void CEquinoxBasedT2C::CeltoTer(double *R, double *dR)
{
	double mat_trans[3][3], dgmat[3][3];
	EARTHFIX2J2000(&mjd, &sod, &mat_trans[0][0], &dgmat[0][0], m_dirconfig, strlen(m_dirconfig));
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		{
			R[i*3+j]=mat_trans[i][j];
			dR[i*3+j] = dgmat[i][j];

		}


}

void CEquinoxBasedT2C::TertoCel(double *R, double *dR)
{
	double mat_trans[3][3], dgmat[3][3];
	EARTHFIX2J2000(&mjd, &sod, &mat_trans[0][0], &dgmat[0][0], m_dirconfig, strlen(m_dirconfig));
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		{
			R[j*3+i]=mat_trans[i][j];
			dR[j*3+i] = dgmat[i][j];

		}

}




double CEquinoxBasedT2C::FromYMDtoSecond(int year, int month, int day, int hour, int minute, double second)
{
	double reftime;
	int ref_year = Tref_year;
	int ref_month = Tref_month;
	int ref_day = Tref_day;
	int ref_hour = Tref_hour;
	int ref_minute= Tref_minute;
	double ref_second = Tref_second;
	reftime=DIFF_TIME(&year, &month, &day, &hour, &minute, &second, &ref_year,
						  &ref_month, &ref_day, &ref_hour, &ref_minute, &ref_second);
	return reftime;	

}

void CEquinoxBasedT2C::FromSecondtoYMD(double refsecond, int &year, int &month, int &day, int &hour, int &minute,double &second)
{

	int ref_year = Tref_year;
	int ref_month = Tref_month;
	int ref_day = Tref_day;
	int ref_hour = Tref_hour;
	int ref_minute= Tref_minute;
	double iref_second = Tref_second;
	SEC2YMDHMS(&ref_year, &ref_month, &ref_day, &ref_hour, &ref_minute,
                         &iref_second, &refsecond, &year, &month, &day, &hour,&minute,
	                     &second);	

}
