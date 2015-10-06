// ResampleKernels.cpp: implementation of the ResampleKernels class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "ResampleKernels.h"
#ifndef PI 
#define PI  3.1415926535897932384626433832795
#endif
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double sqr(double x)
{
	return x*x;
}
double sinc(const double &x) 	
{
	return ((x==0) ? 1 : sin(PI*x)/(PI*x));
}
double rect(double x)
{
	double ans = 0.0;
	if (x<0.5 && x>-0.5) ans = 1;
	else if (x==0.5 || x==-0.5) ans = 0.5;
	return ans;		
}
double tri(double x)
{
	double ans = 0.0;
	if (x<1.0 && x>-1.0)
	{(x<0) ? ans=1+x : ans=1-x;}
	return ans;		
}
void rect(double *x,int N,double *y) 	
{	
	//double *y=new double[N];
	for (register int i=0;i<N;i++)
		y[i] = rect(x[i]);
//	return y;
}
void tri(double *x,int N,double *y)	
{
	//double *y=new double[N];		
	for (register int i=0;i<N;i++)
		y[i] = tri(x[i]);
//	return y;
}
void cc4(double *x,int N,double *y)
{
	double alpha = -1.0;
//	double *y=new double[N];
	for (register int i=0;i<N;i++)
	{
		double xx2 = sqr(x[i]);
		double xx  = sqrt(xx2);
		if      (xx < 1)
			y[i] = (alpha+2)*xx2*xx - (alpha+3)*xx2 + 1;
		else if (xx < 2)
			y[i] = alpha*xx2*xx - 5*alpha*xx2 + 8*alpha*xx - 4*alpha;
		else 
			y[i] = 0.0;
	}
//	return y;		
}
void  cc6(double *x,int N,double *y)
{
	float alpha = -.5;
	float beta  =  .5;
//	double *y=new double[N];
	for (register int i=0;i<N;i++)
	{
		double xx2 = sqr(x[i]);
		double xx  = sqrt(xx2);
		if      (xx < 1)
			y[i] = (alpha-beta+2)*xx2*xx - (alpha-beta+3)*xx2 + 1;			
		else if (xx < 2)
			y[i] =   alpha*xx2*xx - (5*alpha-beta)*xx2 
			+ (8*alpha-3*beta)*xx - (4*alpha-2*beta);
		else if (xx < 3)
			y[i] = beta*xx2*xx - 8*beta*xx2 + 21*beta*xx - 18*beta;
		else 
			y[i] = 0.0;
	}
//	return y;		
}
void  ts6(double *x,int N,double *y)
{
//	double *y=new double[N];
	for (register int i=0;i<N;i++)
	{
		y[i] = sinc(x[i]) * rect(x[i]/6.0);
	}
//	return y;
}
// ___ knab: oversampling factor of signal CHI, number of points N ___
void  knab(double *x,int N,double *y)
{
	double CHI=0.5;
	
//	double *y=new double[N];
	double v  = 1.0-1.0/CHI;
	double vv = PI*v*double(N)/2.0;
	for (register int i=0;i<N;i++)
	{
		y[i] = sinc(x[i])*cosh(vv*sqrt(1.0-sqr(2.0*x[i]/double(N))))/cosh(vv);
	}
//	return y;
	
}
// ___ raised cosine: oversampling factor of signal CHI, number of points N ___
void  rc_kernel(double *x,int N,double *y)
{
	
	double CHI=0.5;	
//	double *y=new double[N];
	double v  = 1.0-1.0/CHI;// alpha in paper cho05
	for (register int i=0;i<N;i++)
	{
		y[i] = sinc(x[i]) * rect(x[i]/double(N))*
			cos(v*PI*x[i]) / (1.0-sqr(2.0*v*x[i]));
	}
//	return y;		
}


double bsl1(int n,double x)
{ 
    int i,m;
    double t,y,p,b0,b1,q;
    static double a[7]={ 1.0,3.5156229,3.0899424,1.2067492,
		0.2659732,0.0360768,0.0045813};
    static double b[7]={ 0.5,0.87890594,0.51498869,
		0.15084934,0.02658773,0.00301532,0.00032411};
    static double c[9]={ 0.39894228,0.01328592,0.00225319,
		-0.00157565,0.00916281,-0.02057706,
		0.02635537,-0.01647633,0.00392377};
    static double d[9]={ 0.39894228,-0.03988024,-0.00362018,
		0.00163801,-0.01031555,0.02282967,
		-0.02895312,0.01787654,-0.00420059};
    if (n<0) n=-n;
    t=fabs(x);
    if (n!=1)
	{ if (t<3.75)
	{ 
		y=(x/3.75)*(x/3.75); p=a[6];
		for (i=5; i>=0; i--)
			p=p*y+a[i];
	}
	else
	{ 
		y=3.75/t; p=c[8];
		for (i=7; i>=0; i--)
			p=p*y+c[i];
		p=p*exp(t)/sqrt(t);
	}
	}
    if (n==0) return(p);
	q=p;
	if (t<3.75)
	{ 
		y=(x/3.75)*(x/3.75); p=b[6];
		for (i=5; i>=0; i--) p=p*y+b[i];
		p=p*t;
	}
    else
	{ 
		y=3.75/t; p=d[8];
   	    for (i=7; i>=0; i--) p=p*y+d[i];
	       p=p*exp(t)/sqrt(t);
	}
    if (x<0.0) p=-p;
    if (n==1) return(p);
    if (x==0.0) return(0.0);
    y=2.0/t; t=0.0; b1=1.0; b0=0.0;
    m=n+(int)sqrt(40.0*n);
    m=2*m;
    for (i=m; i>0; i--)
	{ 
		p=b0+i*y*b1; b0=b1; b1=p;
		if (fabs(b1)>1.0e+10)
		{
			t=t*1.0e-10; b0=b0*1.0e-10;
			b1=b1*1.0e-10;
		}
		if (i==n) t=b0;
	}
    p=t*q/b1;
    if ((x<0.0)&&(n%2==1)) p=-p;
    return(p);
}

double kaiser(double x,int N,double bx)
{
    //默认参数    
	//int bx = 3;	
	double temp0 = (x/N);
	double temp1 = temp0*temp0;
	double temp2 = bx*sqrt(fabs(1-temp1));
	double temp3 = bsl1(0,temp2)/bsl1(0,bx);
	return temp3;
}

void  sc_kaiser_Sinc(double *x,int N,double *y)
{	
	for (register int i=0;i<N;i++)
	{		
		y[i] = sinc(x[i])*kaiser(x[i],N/2,3);
	}
}



/*************************************************************************************
**************************************************************************************/


#include "SatOrbitNew.h"



unsigned ResampleKernels::Ninterval = 128;
//
ResampleKernels::ResampleKernels()
{
	m_Kernel = new double*[Ninterval];
	m_CurKernel = NULL;
}

ResampleKernels::~ResampleKernels()
{
	if ( NULL != m_CurKernel )	delete []m_CurKernel;

	for (int i=0; i<Ninterval; ++i)
	{
		if ( m_Kernel[i] != NULL)
		 delete []m_Kernel[i];
	}
	delete []m_Kernel;
}
//
//由采样方法确定相关采样信息，应该放在初始化代码中
bool ResampleKernels::ReSampleNumBasedReSampleMode(ReSampleMode ResamMode)
{
	switch( ResamMode)
	{
	case Nearest_Neighbor:
		m_ResampleMethod=rs_rect;
		break;
	case Piecewise_Linear:
		m_ResampleMethod=rs_tri;
		break;
	case Cubic_Convolution_4P:
		m_ResampleMethod=rs_cc4p;
		break;
	case Truncated_Sinc_6P:
		m_ResampleMethod=rs_ts6p;
		break;
	case Cubic_Convolution_6P:
		m_ResampleMethod=rs_cc6p;
		break;
	case Knab_6P:
		m_ResampleMethod=rs_knab6p;
		break;
	case Raised_Cosine_6P:
		m_ResampleMethod=rs_rc6p;
		break;
//	case Kaiser_SINC_24P:
//		m_ResampleMethod=rs_kaizer24p;
//		break;
	default:
		m_ResampleMethod=rs_tri;
		return FALSE;
	}
	
	return TRUE;
}

//建立查询表，更新变量m_Kernel
void ResampleKernels::MakeResampleKernelsTable(int method)
{
	const int Npoints      = method%100; 
	const int Npointsd2    = Npoints/2;
	const int Npointsd2m1  = Npointsd2-1;
	register int i,j;
	double *x_axis=new double[Npoints];
	for (i=0; i<Npoints; ++i)
		x_axis[i] = (1.0f - Npointsd2 + i);              // start at [-1 0 1 2]

	const int INTERVAL  = Ninterval - 1;                 // precision: 1./INTERVAL [pixel]
	const double dx        = 1.0/INTERVAL;               // interval look up table

	for (i=0; i<Ninterval; ++i)
	{
		m_Kernel[i]=new double[Npoints];

		switch(method)
		{
			// --- Extremely simple kernels (not good, but fast) ---
		case rs_rect:
			rect(x_axis,Npoints,m_Kernel[i])  ;       
			break;
		case rs_tri:
			tri(x_axis,Npoints,m_Kernel[i])  ; 		
			break;
			// --- Cubic Convolution kernel: theoretical better than truncated sinc. ---	
		case rs_cc4p:
			cc4(x_axis,Npoints,m_Kernel[i])  ; 		
			break;
			// --- Truncated sinc ---
		case rs_ts6p:
			ts6(x_axis,Npoints,m_Kernel[i])  ;         
			break;
			// --- Cubic Convolution kernel: theoretical better than truncated sinc. ---		
		case rs_cc6p:
			cc6(x_axis,Npoints,m_Kernel[i])  ;         
			break;
			// --- KNAB kernel: theoretical better than cubic conv. ---	
		case rs_knab6p:
			knab(x_axis,Npoints,m_Kernel[i])  ; 		
			break;	
			// --- Raised cosine: theoretical best ---
		case rs_rc6p:
			rc_kernel(x_axis,Npoints,m_Kernel[i])  ;         
			break;	
		case rs_kaizer24p:
			sc_kaiser_Sinc(x_axis,Npoints,m_Kernel[i])  ;         
			break;	
		default:
			
			break;
		}

		for( j=0;j<Npoints;j++)
		{
			x_axis[j]  -= dx;    // Note: 'wrong' way (mirrored)
		}
	}
	
	for (i=0; i<Ninterval; i++)
	{	
		
		double sum = 0.0;		
		
		for ( j=0; j<Npoints; j++)
		{			
			sum += m_Kernel[i][j];			
		}
		
		for ( j=0; j<Npoints; j++)
		{
			m_Kernel[i][j] /= sum;			
		}
	}
	delete[] x_axis;
}

bool ResampleKernels::Init(ReSampleMode ResamMode)
{
	if ( !ReSampleNumBasedReSampleMode(ResamMode)) return FALSE;
	MakeResampleKernelsTable(m_ResampleMethod);
	
	//要保证m_CurKernel只分配一次内存，而且在最后会被析构掉
	m_CurKernelSize = m_ResampleMethod%100;
	m_CurKernel  = new double[m_CurKernelSize*m_CurKernelSize];
	
	return TRUE;
}


unsigned ResampleKernels::GetCurKernelSize() const
{
	return m_CurKernelSize;
}

//通过输入当前的点的坐标，返回当前的采样内核，注意，这里获取的是二维内核
//而输入的点位应当位于影像内，这是应当自我保证的，那么就需要告知外部程序，我们采样的窗口的大小.故有GetCurKernelSize()函数
void ResampleKernels::GetCurPixelResampleKernels(double sample,double line,double* CurKernel)
{
	const int Npoints      = m_ResampleMethod%100;  // #pnts interpolator
	const int Npointsd2    = Npoints/2;
	const int Npointsd2m1  = Npointsd2-1;
	
	const int INTERVAL  = Ninterval - 1;                 // precision: 1./INTERVAL [pixel]
	
	const int fl_line = int(line);
	const int fl_sample = int(sample);
	const double linedec = line - fl_line;            // e.g. .35432
	const double sampledec = sample - fl_sample;            // e.g. .5232
	
	const int kernelnoLine   = int(linedec*INTERVAL+0.5); // lookup table index
	const int kernelnoSample   = int(sampledec*INTERVAL+0.5); // lookup table index
	double *kernelline =   (m_Kernel[kernelnoLine]);           // local copy 
	double *kernelsample = (m_Kernel[kernelnoSample]);           // local copy		
	helper.mult(kernelline,kernelsample,m_CurKernel,Npoints,1,Npoints);
	memcpy(CurKernel, m_CurKernel, sizeof(double)*(m_CurKernelSize*m_CurKernelSize));
}


//传入的指针应当有足够的空间
//而且传入的sample应当自己保证位于影像范围内
void ResampleKernels::GetOneDemKernels(double sample,double* CurKernel)
{
	int Npointsd2m1    = m_CurKernelSize/2 -1;
	const int INTERVAL = Ninterval - 1;                 // precision: 1./INTERVAL [pixel]
	const double dx    = 1.0/INTERVAL;               // interval look up table
	
	int fl_sample    = int(sample);
	int firstsample  = fl_sample - Npointsd2m1;            
	double sampledec = sample - fl_sample;
	
	const int kernelnoSample   = int(sampledec*INTERVAL+0.5);  // lookup table index
	memcpy(CurKernel, m_Kernel[kernelnoSample],sizeof(double)*m_CurKernelSize);
}
