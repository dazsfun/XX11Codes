// GridInterpolate.h: interface for the CGridInterpolate class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GRIDINTERPOLATE_H__2B050833_F23B_4CB0_A687_6D586C10848F__INCLUDED_)
#define AFX_GRIDINTERPOLATE_H__2B050833_F23B_4CB0_A687_6D586C10848F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "GeoRSDataStruct.h"

//##ModelId=4E4A89B30203
class CGridInterpolate  
{
private:
	//##ModelId=4E4A89B30205
	RPCGCPSTRUCT m_CurGRID;
public:
	//##ModelId=4E4A89B30209
	bool Init(RPCGCPSTRUCT CurGRID);
	//##ModelId=4E4A89B30214
	bool GetZ(double x,double y,double &ddx,double &ddy);
	//##ModelId=4E4A89B30219
	CGridInterpolate();
	//##ModelId=4E4A89B3021A
	virtual ~CGridInterpolate();

};

#endif // !defined(AFX_GRIDINTERPOLATE_H__2B050833_F23B_4CB0_A687_6D586C10848F__INCLUDED_)
