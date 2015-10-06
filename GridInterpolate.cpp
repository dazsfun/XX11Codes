// GridInterpolate.cpp: implementation of the CGridInterpolate class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "GridInterpolate.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//##ModelId=4E4A89B30219
CGridInterpolate::CGridInterpolate()
{

}

//##ModelId=4E4A89B3021A
CGridInterpolate::~CGridInterpolate()
{

}

//##ModelId=4E4A89B30209
bool CGridInterpolate::Init(RPCGCPSTRUCT CurGRID)
{
	m_CurGRID=CurGRID;
	return true;
}
//##ModelId=4E4A89B30214
bool CGridInterpolate::GetZ(double x,double y,double &ddx,double &ddy)
{
	short	lbGridx,lbGridy;
	long	lbOffset,ltOffset;
	double	dx,dy;
	double   z00,z10,z01,z11;

	if ( m_CurGRID.ddx == NULL ) 
		return  false;
	if ( m_CurGRID.ddy == NULL ) 
		return  false;	

	x-=m_CurGRID.x0;                                 // 这里应该还要减去初始偏移
	y-=m_CurGRID.y0;                                 // 这里也是一样的

	x/=m_CurGRID.dx;
	y/=m_CurGRID.dy;

	lbGridx =(short) (x);	dx =  x - lbGridx;
	lbGridy =(short) (y);	dy =  y - lbGridy;
	
	lbOffset = lbGridy * m_CurGRID.nx + lbGridx;
	ltOffset = lbOffset + m_CurGRID.nx;

	z00 = m_CurGRID.ddx[lbOffset];
	z01 = m_CurGRID.ddx[lbOffset+1];
	z10 = m_CurGRID.ddx[ltOffset];
	z11 = m_CurGRID.ddx[ltOffset+1];
	

	z00 += dx*(z01-z00);	/* px0 */
	z10 += dx*(z11-z10);	/* px1 */

	ddx=z00 + dy*(z10 - z00);


	z00 = m_CurGRID.ddy[lbOffset];
	z01 = m_CurGRID.ddy[lbOffset+1];
	z10 = m_CurGRID.ddy[ltOffset];
	z11 = m_CurGRID.ddy[ltOffset+1];
	

	z00 += dx*(z01-z00);	/* px0 */
	z10 += dx*(z11-z10);	/* px1 */

	ddy=z00 + dy*(z10 - z00);

	return true;

	

}
