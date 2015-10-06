// EquinoxBasedT2C.h: interface for the CEquinoxBasedT2C class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_EQUINOXBASEDT2C_H__CDD852A8_A448_478B_802E_DD98A6245B4A__INCLUDED_)
#define AFX_EQUINOXBASEDT2C_H__CDD852A8_A448_478B_802E_DD98A6245B4A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <afx.h>
#include "Safeauxilary.h"

class  CEquinoxBasedT2C  
{
private:
	char m_dirconfig[1024];
	int mjd;
	double sod;

public:
	void FromSecondtoYMD(double refsecond,int &year,int &month,int &day,int &hour,int &minute,double &second);
	double FromYMDtoSecond(int year,int month,int day,int hour,int minute,double second);
	void TertoCel(double *R, double *dR);
	void CeltoTer(double *R,double *dR);
	void Inittime(int year,int month,int day,int hour,int minute ,double second);
	BOOL Initdir(CString path,CString ipath);


	CEquinoxBasedT2C();
	virtual ~CEquinoxBasedT2C();

};

#endif // !defined(AFX_EQUINOXBASEDT2C_H__CDD852A8_A448_478B_802E_DD98A6245B4A__INCLUDED_)
