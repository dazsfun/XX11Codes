// ResampleKernels.h: interface for the ResampleKernels class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RESAMPLEKERNELS_H__4098CA8E_AA11_4C45_976E_B4B4349705B9__INCLUDED_)
#define AFX_RESAMPLEKERNELS_H__4098CA8E_AA11_4C45_976E_B4B4349705B9__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Curve.h"
#include "SatOrbitNew.h"
#include "GeoRSDataStruct.h"


const int 	rs_rect		= 102;		// nearest neighbor, 1 or 2 points
const int 	rs_tri		= 202;		// piecewize linear, 2 point
const int 	rs_cc4p		= 304;		// cubic convolution, 4 point
const int 	rs_cc6p		= 306;		// cubic convolution, 6 point
const int 	rs_ts6p		= 406;		// truncated sinc, 6 point
const int 	rs_knab6p	= 506;		// knab window, 6 point
const int 	rs_rc6p  	= 606;		// raised cosine, 6 point
const int 	rs_kaizer24p = 724;		// raised cosine, 6 point





class ResampleKernels  
{
	CSatOrbit helper;
	static unsigned Ninterval;
	
	int m_ResampleMethod;
	int m_CurKernelSize;                      //in sample and line direcetion
	
	double **m_Kernel;
//	double *m_Kernel[128];
	double *m_CurKernel;
private:
	//获取采样方式
	bool ReSampleNumBasedReSampleMode(ReSampleMode ResamMode);

	//获取采样内核
	void MakeResampleKernelsTable(int method);
	
public:
	ResampleKernels();
	virtual ~ResampleKernels();

public:
	bool Init(ReSampleMode ResamMode);
	unsigned GetCurKernelSize() const;

	void GetCurPixelResampleKernels(double sample,double line,double* CurKernel);	
	void GetOneDemKernels(double sample,double* CurKernel);	

};

void rect(double *x,int N,double *y);
void tri(double *x,int N,double *y);
void cc4(double *x,int N,double *y);
void ts6(double *x,int N,double *y);
void cc6(double *x,int N,double *y);
void knab(double *x,int N,double *y);
void rc_kernel(double *x,int N,double *y);
void sc_kaiser_Sinc(double *x,int N,double *y);

#endif // !defined(AFX_RESAMPLEKERNELS_H__4098CA8E_AA11_4C45_976E_B4B4349705B9__INCLUDED_)
