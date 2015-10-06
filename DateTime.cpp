


#include <math.h>
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
//#include <stdio.h>
#include <string>
#include "stdafx.h"
using namespace std;
#include "DateTime.h"
bool _isLeapYear ( int year )
{
// check if the year is a leap year

	if ( year % 4 != 0 )
		return false ;
	if ( (year%100==0) && (year%400!=0) )
		return false ;

	return true ;
}//ll add

CDateTime::CDateTime(void)
{
	
}
CDateTime::~CDateTime(void)
{
}

void CDateTime::InitDateRef(unsigned short Satellite_ID)
{
	m_Satellite_ID=Satellite_ID;
	if(Satellite_ID==0xC101)
	{		
		ref_year = 1998 ;
		ref_month = 1 ;
		ref_day = 1 ;
		ref_hour = 0 ;
		ref_minute = 0 ;
		ref_second = 0.0 ;
	}
	else if (Satellite_ID == 0xC102)
	{
		ref_year = 2009 ;
		ref_month = 1 ;
		ref_day = 1 ;
		ref_hour = 0 ;
		ref_minute = 0 ;
		ref_second = 0.0 ;
	}
	else if (Satellite_ID == 0xC103)
	{
		ref_year = 2009;
		ref_month = 1;
		ref_day = 1;
		ref_hour = 0;
		ref_minute = 0;
		ref_second = 0.0;
	}
}

CDateTime :: CDateTime (unsigned short Satellite_ID, int yr, int mon, int d, int hr, int min, double sec )
{
// constructor : format one
// the DateTime value must non-negative, otherwise, the value will be set with default 
	
	InitDateRef( Satellite_ID);
	setDateTime ( yr, mon, d, hr, min, sec ) ;
}

CDateTime :: CDateTime (unsigned short Satellite_ID, char *str ) 
{
// constructor : format two
// class constructor with a str : "yyyy:mm:dd:hh:mm:sss..."

	InitDateRef( Satellite_ID);
	setDateTime ( str ) ;
}

CDateTime :: CDateTime (unsigned short Satellite_ID, double inSeconds )
{
// constructor : format three
// using the timecode in seconds starting from the reference date
	
	InitDateRef( Satellite_ID);
	year = ref_year ;	month = ref_month ;		day = ref_day ;
	hour = ref_hour ;	minute = ref_minute ;	second = ref_second ;
		
	countSecond2Date ( this, inSeconds ) ;
}


void CDateTime :: setDateTime ( int yr, int mon, int d, int hr, int min, double sec )
{
	year	=  yr >= 1 ?  yr : 1 ;
	month	= mon >= 1 ? mon : 1 ; 
	day		=   d >= 1 ?   d : 1 ;
	hour	=  hr >= 0 ?  hr : 0 ;
	minute	= min >= 0 ? min : 0 ;
	second	= sec >= 0 ? sec : 0 ;
}

void CDateTime :: setDateTime ( char *str )
{
	int		i, j ;
	int		month_day[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 } ;
	char	buf[16] ;

// set the default value
	year = 1 ;	month = 1 ;		day = 1 ;
	hour = 0 ;	minute = 0 ;	second = 0 ;

	if ( (int)strlen(str) == 0 ) return ;

	i = j = 0 ;
	while ( (str[i]!=':') && (i<(int)strlen(str)) && (j<15) ) 
	{
		buf[j] = str[i] ;
		i++; j++ ;
	}
	buf[j] = '\0' ;
	year = atoi ( buf ) ;
	if ( year < 1 )	year = 1 ;
	if ( _isLeapYear(year) ) month_day[1] = 29 ;
	if ( i == (int)strlen(str) ) return ;

	i++ ;	j = 0 ;
	while ( (str[i]!=':') && (i<(int)strlen(str)) && (j<15) ) 
	{
		buf[j] = str[i] ;
		i++; j++ ;
	}
	buf[j] = '\0' ;
	month = atoi ( buf ) ;
	if ( (month<1) || (month>12) )	month = 1 ;
	if ( i == (int)strlen(str) ) return ;

	i++ ;	j = 0 ;
	while ( (str[i]!=':') && (i<(int)strlen(str)) && (j<15) ) 
	{
		buf[j] = str[i] ;
		i++; j++ ;
	}
	buf[j] = '\0' ;
	day = atoi ( buf ) ;
	if ( (day<1) || (day>month_day[month-1] ) )	day = 1 ;
	if ( i == (int)strlen(str) ) return ;

	i++ ;	j = 0 ;
	while ( (str[i]!=':') && (i<(int)strlen(str)) && (j<15) ) 
	{
		buf[j] = str[i] ;
		i++; j++ ;
	}
	buf[j] = '\0' ;
	hour = atoi ( buf ) ;
	if ( (hour<0) || (hour>23) )	hour = 0 ;
	if ( i == (int)strlen(str) ) return ;

	i++ ;	j = 0 ;
	while ( (str[i]!=':') && (i<(int)strlen(str)) && (j<15) ) 
	{
		buf[j] = str[i] ;
		i++; j++ ;
	}
	buf[j] = '\0' ;
	minute = atoi ( buf ) ;
	if ( (minute<0) || (minute>59) )	minute = 0 ;
	if ( i == (int)strlen(str) ) return ;

	i++ ;	j = 0 ;
	while ( (i<(int)strlen(str)) && (j<15) ) 
	{
		buf[j] = str[i] ;
		i++; j++ ;
	}
	buf[j] = '\0' ;
	second = atof ( buf ) ;
	if ( (second<0) || (second>60) )	second = 0 ;

	return ;
}

void CDateTime :: outputDate2String ( char *str ) 
{
// output the date is the format of "yyyy:mm:dd:hh:mm:ss.ssssss"

	sprintf_s( str,20, "%04d:%02d:%02d:%02d:%02d:%09.6f", year, month, day, hour, minute, second ) ;
	
}

double CDateTime :: countDate2Second ( ) 
{
// count the total second of the DateTime from the reference date which is set in design stage
	
	double		timeCode ;
	int			i ;
	const int	month_day[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 } ;
	
	timeCode = 0 ;

// if the date is after reference date, the return value is positive
// otherwise is negative

	if ( (*this) > CDateTime(m_Satellite_ID,ref_year,ref_month,ref_day,ref_hour,ref_minute,ref_second) )
	{
		if ( (*this) > CDateTime(m_Satellite_ID,ref_year+1,1,1) )
		{
			timeCode += month_day[ref_month-1] - ref_day + 1 ;
			for ( i=ref_month+1; i<=12; i++ )
			{
				timeCode += month_day[i-1] ;
				if ( _isLeapYear(ref_year) && (i==2) ) 
					timeCode ++ ;
			}

			for ( i=ref_year+1; i<year; i++ )
			{
				timeCode += 365 ;
				if ( _isLeapYear(i) ) timeCode++ ;
			}
			for ( i=1; i<month; i++ )
			{
				timeCode += month_day[i-1] ;
				if ( _isLeapYear(year) && (i==2) ) timeCode++ ;
			}
		}
		else
		{
			timeCode += month_day[ref_month-1] - ref_day + 1 ;
			for ( i=ref_month+1; i<month; i++ ) 
			{
				timeCode += month_day[i-1] ;
				if ( _isLeapYear(year) && (i==2) ) 
					timeCode ++ ;
			}
		}
		timeCode += day - 1 ;
		timeCode *= 3600 * 24 ;

		timeCode += ( hour - ref_hour ) * 3600 ;
		timeCode += ( minute - ref_minute ) * 60 ;
		timeCode += second - ref_second ;
	}

	if ( (*this) == CDateTime(m_Satellite_ID,ref_year,ref_month,ref_day,ref_hour,ref_minute,ref_second) )
		timeCode = 0 ;

	if ( (*this) < CDateTime(m_Satellite_ID,ref_year,ref_month,ref_day,ref_hour,ref_minute,ref_second) )
	{
		if ( CDateTime(m_Satellite_ID,ref_year,ref_month,ref_day) > CDateTime(m_Satellite_ID,year+1,1,1) )
		{
			timeCode += month_day[month-1] - day + 1 ;
			for ( i=month+1; i<=12; i++ )
			{
				timeCode += month_day[i-1] ;
				if ( _isLeapYear(year) && (i==2) ) 
					timeCode ++ ;
			}

			for ( i=year+1; i<ref_year; i++ )
			{
				timeCode += 365 ;
				if ( _isLeapYear(i) ) timeCode++ ;
			}
			for ( i=1; i<ref_month; i++ )
			{
				timeCode += month_day[i-1] ;
				if ( _isLeapYear(ref_year) && (i==2) ) timeCode++ ;
			}
		}
		else
		{
			timeCode += month_day[month-1] - day + 1 ;
			for ( i=month+1; i<ref_month; i++ ) 
			{
				timeCode += month_day[i-1] ;
				if ( _isLeapYear(ref_year) && (i==2) ) 
					timeCode ++ ;
			}
		}
		timeCode += ref_day - 1 ;
		timeCode *= 3600 * 24 ;

		timeCode += ( ref_hour - hour ) * 3600 ;
		timeCode += ( ref_minute - minute ) * 60 ;
		timeCode += ref_second - second ;

		timeCode *= (-1) ;
	}

	return timeCode ;
}

void CDateTime :: reset (unsigned short Satellite_ID )
{
// reset DateTime to the reference date	
	InitDateRef(Satellite_ID);
}


double	operator -  ( CDateTime a, CDateTime b ) 
{
	double asecond=a.countDate2Second();
	double bsecond=b.countDate2Second();
	return asecond-bsecond;

}
bool operator > ( CDateTime a, CDateTime b )
{
// friend function of DateTime, it will return true if date "a" is after date "b"
	
	if ( a.year > b.year ) return true ;
	if ( a.year < b.year ) return false ;
	if ( a.year == b.year )
	{
		if ( a.month > b.month ) return true ;
		if ( a.month < b.month ) return false ;
		if ( a.month == b.month ) 
		{
			if ( a.day > b.day ) return true ;
			if ( a.day < b.day ) return false ;
			if ( a.day == b.day )
			{
				if ( a.hour > b.hour ) return true ;
				if ( a.hour < b.hour ) return false ;
				if ( a.hour == b.hour ) 
				{
					if ( a.minute > b.minute ) return true ;
					if ( a.minute < b.minute ) return false ;
					if ( a.minute == b.minute ) 
					{
						if ( a.second > b.second + delta5 ) return true ;
						else return false ;
					}
				}
			}
		}
	}
	return true;
}

bool operator < ( CDateTime a, CDateTime b )
{
	if ( a.year < b.year ) return true ;
	if ( a.year > b.year ) return false ;
	if ( a.year == b.year )
	{
		if ( a.month < b.month ) return true ;
		if ( a.month > b.month ) return false ;
		if ( a.month == b.month ) 
		{
			if ( a.day < b.day ) return true ;
			if ( a.day > b.day ) return false ;
			if ( a.day == b.day )
			{
				if ( a.hour < b.hour ) return true ;
				if ( a.hour > b.hour ) return false ;
				if ( a.hour == b.hour ) 
				{
					if ( a.minute < b.minute ) return true ;
					if ( a.minute > b.minute ) return false ;
					if ( a.minute == b.minute )
					{
						if ( a.second < b.second - delta5 ) return true ;
						else return false ;
					}
				}
			}
		}
	}
	return true;
}

bool operator == ( CDateTime a, CDateTime b ) 
{
	if ( a.year != b.year ) return false ;
	else
	{
		if ( a.month != b.month ) return false ;
		else
		{
			if ( a.day != b.day ) return false ;
			else
			{
				if ( a.hour != b.hour ) return false ;
				else
				{
					if ( a.minute != b.minute ) return false ;
					else
					{
						if ( fabs(a.second-b.second) < delta5 ) return true ;
						else return false ;
					}
				}
			}
		}
	}
}

void countDay2Date ( CDateTime *theDay, long days )
{
// count the date from current date after/before specific days 
// after current date if days > 0, before current date if days < 0 
// suppose current time is 00:00:00

	const int	month_day[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 } ;
	long		daysCount ;
	int			dd ;	
	
	if ( days > 0 )
	{
		daysCount = 0 ;
		dd = month_day[theDay->month-1] - theDay->day + 1 ;
		if ( _isLeapYear(theDay->year) && (theDay->month==2) ) dd++ ;
		while ( daysCount + dd <= days )
		{
			daysCount += dd ;
			
			theDay->day = 1 ;
			theDay->month++ ;
			if ( theDay->month > 12 ) 
			{
				theDay->month = 1 ;
				theDay->year++ ;
			}
		
			dd = month_day[theDay->month-1] ;
			if ( _isLeapYear(theDay->year) && (theDay->month==2) ) dd = 29 ;
		}
		theDay->day += days - daysCount ;
	}
	if ( days < 0 ) 
	{
		daysCount = -days ; 
		theDay->day = theDay->day - 1 ;
		while ( daysCount - theDay->day > 0 )
		{
			daysCount -= theDay->day ;
			
			theDay->month-- ;
			if ( theDay->month < 1 ) 
			{
				theDay->month = 12 ;
				theDay->year-- ;
			}
			theDay->day = month_day[theDay->month-1] ;
			if ( _isLeapYear(theDay->year) && (theDay->month==2) ) theDay->day = 29 ;
		}
		theDay->day = theDay->day - daysCount + 1 ;
	}
}


void timestring2yearmonthdayhourminutessecond(string timestring,int &year,int &month,int &day,int &hour,int &minute,double &second)
{
	//20071031:122100.33
	year=atoi(timestring.substr(0,4).c_str());
	month=atoi(timestring.substr(4,2).c_str());
	day=atoi(timestring.substr(6,2).c_str());
	hour=atoi(timestring.substr(9,2).c_str());
	minute=atoi(timestring.substr(11,2).c_str());
	second=atof(timestring.substr(13).c_str());
	return;
}

void countSecond2Date ( CDateTime *theDay, double seconds )
{
// count the date after/before a certain seconds from current date 
// current date is specified by theDay
// if seconds < 0 --- before current date, seconds > 0 --- after current date
	
	long	ltime ;
	long	days ;
	int		carry ;

	if ( seconds > 0 )
	{
		ltime = (long)seconds ;
		days = ltime / ( 3600 * 24 ) ;

		theDay->second += ( ltime % 60 ) + ( seconds - ltime ) ;
		if ( theDay->second >= 60 )
		{
			carry = 1 ;
			theDay->second -= 60 ;
		}
		else
			carry = 0 ;

		theDay->minute += ( ltime / 60 ) % 60 + carry ;
		if ( theDay->minute > 59 )
		{
			carry = 1 ;
			theDay->minute -= 60 ;
		}
		else
			carry = 0 ;
	
		theDay->hour += ( ltime / 3600 ) % 24 + carry ;
		if ( theDay->hour > 23 )
		{
			carry = 1 ;
			theDay->hour -= 24 ;
		}
		else
			carry = 0 ;

		days += carry ;
		countDay2Date ( theDay, days ) ;
	}

	if ( seconds < 0 )
	{
		seconds = -seconds ;
		ltime = (long)(seconds) ;
		days = ltime / ( 3600 * 24 ) ;

		theDay->second = theDay->second - ( ltime % 60 ) - ( seconds - ltime ) ;
		if ( theDay->second < 0 )
		{
			carry = 1 ;
			theDay->second += 60 ;
		}
		else
			carry = 0 ;

		theDay->minute = theDay->minute - ( ltime / 60 ) % 60 - carry ;
		if ( theDay->minute < 0 )
		{
			carry = 1 ;
			theDay->minute += 60 ;
		}
		else
			carry = 0 ;
	
		theDay->hour = theDay->hour - ( ltime / 3600 ) % 24 - carry ;
		if ( theDay->hour < 0 )
		{
			carry = 1 ;
			theDay->hour += 24 ;
		}
		else
			carry = 0 ;

		days = ( days + carry ) * ( -1 ) ;
		countDay2Date ( theDay, days ) ;
	}
}


void ntohl(char *ptr,int n)
{
	char *ptr_o,ch;
	ptr_o = ptr+n-1;
	int i=0;
	for(i=0; i<n/2; i++)
	{
		ch = *ptr;
		*ptr = *ptr_o;
		*ptr_o = ch;
		ptr++;
		ptr_o--;			
	}
}
