#ifndef _DATETIME_H_
#define _DATETIME_H_
#include <string>
using namespace std;

#define  delta5 0.0001
class CDateTime  
{
private:
	int		ref_year ;
	int		ref_month ;
	int		ref_day ;
	int		ref_hour ;
	int		ref_minute ;
	double	ref_second ;	
public:
	int		year ;
	int		month ;
	int		day ;
	int		hour ;
	int		minute ;
	double	second ;
	unsigned short m_Satellite_ID;
private:
	
	void InitDateRef(unsigned short Satellite_ID);
	
public:
	CDateTime (void);
	
	CDateTime ( unsigned short Satellite_ID,int yr, int mon, int d, int hr = 0, int min = 0, double sec = 0 );
	CDateTime ( unsigned short Satellite_ID,char *str ) ;
	CDateTime ( unsigned short Satellite_ID,double inSeconds = 0 ) ;
	//virtual ~DateTime();
	
	void	outputDate2String ( char *str ) ;
	double	countDate2Second ( ) ;
	void	setDateTime ( int yr, int mon, int d, int hr = 0 , int min = 0, double sec = 0 ) ;
	void	setDateTime ( char *str ) ;
	void	reset ( unsigned short Satellite_ID) ;	

	
	virtual ~CDateTime();
};

bool	operator >  ( CDateTime a, CDateTime b ) ;
bool	operator <  ( CDateTime a, CDateTime b ) ;
bool	operator == ( CDateTime a, CDateTime b ) ;
double	operator -  ( CDateTime a, CDateTime b ) ;
void	countDay2Date ( CDateTime *theDay, long days ) ;
void	countSecond2Date ( CDateTime *theDay, double seconds ) ;
void timestring2yearmonthdayhourminutessecond(string timestring,int &year,int &month,int &day,int &hour,int &minute,double &second);
void ntohl(char *ptr,int n);

#endif /* _DATETIME_H_ */
