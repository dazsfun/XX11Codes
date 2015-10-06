#include "stdAfx.h"
#include "Curve.h"
#include <algorithm>
#include "SatOrbitNew.h"

#ifndef PI 
#define PI  3.1415926535897932384626433832795
#endif




const double ESP = 1e-5;




bool IsInRect(const DRect& WinSize,const DPoint2D& opt) 
{
	return opt.x <= WinSize.right && opt.x >= WinSize.left && opt.y >= WinSize.top && opt.y <= WinSize.bottom;
}


bool IsInRect(const DRect& WinSize,double x, double y) 
{
	return x <= WinSize.right && x >= WinSize.left && y >= WinSize.top && y <= WinSize.bottom;
}



//////////////////////////////////////////////////////////////////////////

//##ModelId=4E4A89B30332
DPolygon::DPolygon()
{
}

//##ModelId=4E4A89B30333
DPolygon::DPolygon(DPoint2D* pPo, unsigned nNum)
{
	vSides.clear();
	vPoints.clear();
	for (unsigned i=0; i<nNum; ++i) 
	{
		LineSegment pSi;
		
		pSi.pt1 = pPo[i];
		pSi.pt2 = pPo[(i+1)%nNum];
		vSides.push_back( pSi);
		vPoints.push_back( pPo[i]);
	}
}

//##ModelId=4E4A89B3033E
DPolygon::~DPolygon()
{
}

//##ModelId=4E4A89B3033C
DPolygon::DPolygon(const std::vector<DPoint2D>& vPo )
{
	vSides.clear();
	vPoints.clear();
	vPoints = vPo;
	size_t nSize = vPo.size();
	for (unsigned i=0; i<nSize; ++i) 
	{
		LineSegment pSi;
		pSi.pt1 = vPo[i];
		pSi.pt2 = vPo[(i+1)%nSize];
		vSides.push_back( pSi);
	}
}

//##ModelId=4E4A89B30340
DPolygon::DPolygon(const DPolygon& rhs)
{
	vPoints.assign( rhs.vPoints.begin(), rhs.vPoints.end() );
	vSides.assign( rhs.vSides.begin(), rhs.vSides.end() );
}

//##ModelId=4E4A89B30342
DPolygon& DPolygon::operator =(const DPolygon& rhs)
{
	if( this == &rhs) return *this;
	else 
	{
		vSides = rhs.vSides;
		vPoints = rhs.vPoints;
	}
	return *this;
}

bool IsInPoygon(const DPolygon& WinSize,const DPoint2D& opt)
{
	return WinSize.InPolygon( opt);
}

double Multiply(DPoint2D p1, DPoint2D p2, DPoint2D p0)
{
	return ( (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y) );
}

// 判断线段是否包含点point
bool IsOnline(DPoint2D point, LineSegment line)
{
	return( ( fabs(Multiply(line.pt1, line.pt2, point)) < ESP ) &&
		( ( point.x - line.pt1.x ) * ( point.x - line.pt2.x ) <= 0 ) &&
		( ( point.y - line.pt1.y ) * ( point.y - line.pt2.y ) <= 0 ) );
}

// 判断线段相交
bool Intersect(LineSegment L1, LineSegment L2)
{
	return( (max(L1.pt1.x, L1.pt2.x) >= min(L2.pt1.x, L2.pt2.x)) &&
		(max(L2.pt1.x, L2.pt2.x) >= min(L1.pt1.x, L1.pt2.x)) &&
		(max(L1.pt1.y, L1.pt2.y) >= min(L2.pt1.y, L2.pt2.y)) &&
		(max(L2.pt1.y, L2.pt2.y) >= min(L1.pt1.y, L1.pt2.y)) &&
		(Multiply(L2.pt1, L1.pt2, L1.pt1) * Multiply(L1.pt2, L2.pt2, L1.pt1) >= 0) &&
		(Multiply(L1.pt1, L2.pt2, L2.pt1) * Multiply(L2.pt2, L1.pt2, L2.pt1) >= 0)
		);
}

// 判断点在多边形内
//##ModelId=4E4A89B3034B
bool DPolygon::InPolygon(const DPoint2D& point) const
{
	unsigned nPointsNum = vPoints.size();
    if (nPointsNum == 1)
	{
		return ( (fabs(vPoints[0].x - point.x) < ESP) && (fabs(vPoints[0].y - point.y) < ESP) );
	}
	else if (nPointsNum == 2) 
	{
		return IsOnline(point, vSides[0]);
	}
	
	int count = 0;
	LineSegment line;
	line.pt1 = point;
	line.pt2.y = point.y;
	line.pt2.x = - INFINITY;
	
	for( unsigned i = 0; i < nPointsNum; i++ ) 
	{
		if( IsOnline(point, vSides[i]) ) return TRUE;
		
		// 如果side平行x轴则不作考虑
		if( fabs(vSides[i].pt1.y - vSides[i].pt2.y) < ESP )
		{
			continue;
		}
		
		if( IsOnline(vSides[i].pt1, line) )
		{
			if( vSides[i].pt1.y > vSides[i].pt2.y ) count++;
        } 
		else if( IsOnline(vSides[i].pt2, line) ) 
		{
			if( vSides[i].pt2.y > vSides[i].pt1.y ) count++;
		} 
		else if( Intersect(line, vSides[i]) )  count++;
    }
	
    return ( count % 2 == 1 );
}



//
#include "Polygon.h"
unsigned PolygonIntersection(int n1,DPoint2D *p1,int n2,DPoint2D *p2,DPoint2D* p) 
{
		int i =0;
//	const unsigned& n1 = nPointsNum;
//	const unsigned& n2 = otherpolygon.nPointsNum;
//	const DPoint2D *p1 = pPoints;
//	const DPoint2D *p2 = otherpolygon.pPoints;
	unsigned n;

	gpc_polygon  pp1,pp2,pp3;
	pp1.num_contours=1;
	pp1.contour=new gpc_vertex_list[1];
	pp1.hole=new int[1];
	pp1.hole[0]=FALSE;
	pp1.contour->num_vertices=n1;
	pp1.contour->vertex=new gpc_vertex[n1];
	for(int i=0;i<n1;i++)
	{
		pp1.contour->vertex[i].x=p1[i].x;
		pp1.contour->vertex[i].y=p1[i].y;
	}
	
	pp2.num_contours=1;
	pp2.contour=new gpc_vertex_list[1];
	pp2.hole=new int[1];
	pp2.hole[0]=FALSE;
	pp2.contour->num_vertices=n2;
	pp2.contour->vertex=new gpc_vertex[n2];
	for( i=0;i<n2;i++)
	{
		pp2.contour->vertex[i].x=p2[i].x;
		pp2.contour->vertex[i].y=p2[i].y;
	}
	
	gpc_polygon_clip(GPC_INT,&pp1,&pp2,&pp3);
	
	delete pp2.contour->vertex;pp2.contour->vertex=NULL;
	delete pp1.contour->vertex;pp1.contour->vertex=NULL;
	
	delete pp1.contour;pp1.contour=NULL;
	delete pp2.contour;pp2.contour=NULL;
	
	
	if(pp3.contour->num_vertices==0)
	{
		n=0;
		return NULL;
	}
//	DPoint2D *p;
//	p=new DPoint2D[pp3.contour->num_vertices];
	for(i=0;i<pp3.contour->num_vertices;i++)
	{
		p[i].x=pp3.contour->vertex[i].x;
		p[i].y=pp3.contour->vertex[i].y;
	}
	n=pp3.contour->num_vertices;
	
	delete pp3.contour->vertex;pp3.contour->vertex=NULL;
	
	delete pp3.contour;pp3.contour=NULL;
	
	
	return n;
}


//##ModelId=4E4A89B3034E
unsigned DPolygon::GetPolygon(DPoint2D* p) const
{
	size_t num = vPoints.size();
	for (unsigned i=0; i<num; ++i )
	{
		p[i] = vPoints[i];
	}

	return num;
}



//计算直线的系数
//##ModelId=4E4A89B303BC
bool DirectLine::Init(const std::vector<DPoint2D>&  p_pts)
{
	size_t num = p_pts.size();
	if ( num < 2 )	return FALSE;
	
	double NE[4],L[2];
	memset(NE, 0, sizeof(double)*4);
	memset(L , 0, sizeof(double)*2);
	double A[2];
	int i;

	//计算所有点的x y的最大最小值
	double xmax = p_pts[0].x;
	double xmin = p_pts[0].x;
	double ymax = p_pts[0].y;
	double ymin = p_pts[0].y;
	for( i=0;i<num;++i)
	{
		xmax = p_pts[i].x > xmax ? p_pts[i].x : xmax;
		xmin = p_pts[i].x < xmin ? p_pts[i].x : xmin;	
		ymax = p_pts[i].y > ymax ? p_pts[i].y : ymax;
		ymin = p_pts[i].y < ymin ? p_pts[i].y : ymin;		
	}
	
	//若是x的差值大于y上面的差值，则选y = k*x + b
	if ( (xmax - xmin) >= (ymax - ymin) )
	{
		for( i=0;i<num;++i)
		{
			A[0]=p_pts[i].x;
			A[1]=1;
			
			helper.pNormal(A,2,p_pts[i].y,NE,L,1.0);		
		}
		helper.Gauss(NE, L, 2);
		a_ = L[0];
		b_ = -1;
		c_ = L[1];	
	}
	//否则，选择 x = k*y + b
	else
	{
		for( i=0;i<num;++i)
		{
			A[0]=p_pts[i].y;
			A[1]=1;
			
			helper.pNormal(A,2,p_pts[i].x,NE,L,1.0);		
		}
		helper.Gauss(NE, L, 2);
		a_ = -1;
		b_ = L[0];
		c_ = L[1];		
	}
/*
	for( i=0;i<num;++i)
	{
		A[0]=p_pts[i].x;
		A[1]=1;
		
		pNormal(A,2,p_pts[i].y,NE,L,1.0);		
	}

	//需要考虑到无解的情况，就是当k为∞的时候
	if ( fabs(NE[2]*NE[1] - NE[0]*NE[3]) < 0.00000001)
	{
		a_ = -1;
		b_ = 0;
		c_ = 0;
		for ( unsigned i=0; i<p_pts.size(); ++i)
		{
			c_+= p_pts[i].x;
		}
		c_ /= num;
	}
	else
	{
		Gauss(NE, L, 2);
		a_ = L[0];
		b_ = -1;
		c_ = L[1];
	}
*/
	return TRUE;
}

// y = k*x + b 参数为矩形窗口
//这里返回的两个端点并没有顺序性
//##ModelId=4E4A89B303CD
std::vector<DPoint2D> DirectLine::GetEndpoint(DRect& WinSize) const
{
	std::vector<DPoint2D> vec_p(2);
	if (a_ == 0)
	{
		if ( b_==0)
		{
			throw "无最佳直线";                             //抛出异常
		}
		else
		{
			vec_p[0].x = 0;
			vec_p[0].y = vec_p[1].y = c_/(-b_);
			vec_p[1].x = WinSize.right;
		}
	}
	else if ( b_==0)
	{
		vec_p[0].y = 0;
		vec_p[0].x = vec_p[1].x = c_/(-a_);
		vec_p[1].y = WinSize.bottom;
	}
	else
	{
		//计算出四个点，然后判断是否位于矩形内，需要小心的是可能存在两个点相同
		//即存在(0,0)、(0,0)、(135,135)，那么需要判断点是否相同，即保证相同的点只在程序中出现一次
		std::vector<DPoint2D> vp;
		DPoint2D temp1(WinSize.left, (a_*WinSize.left+c_)/(-b_) );
		if ( IsInRect(WinSize, temp1) )
		{
			vp.push_back(temp1);
		}
		DPoint2D temp2(WinSize.right, (a_*WinSize.right+c_)/(-b_));
		if ( IsInRect(WinSize, temp2)&& std::count(vp.begin(), vp.end(), temp2)==0 )
		{
			vp.push_back(temp2);
		}
		DPoint2D temp3( (b_*WinSize.top+c_)/(-a_) ,WinSize.top);
		if( IsInRect(WinSize, temp3) && std::count(vp.begin(), vp.end(), temp3)==0 )
		{
			vp.push_back(temp3);
		}
		DPoint2D temp4((b_*WinSize.bottom+c_)/(-a_) ,WinSize.bottom);
		if( IsInRect(WinSize, temp4) && std::count(vp.begin(), vp.end(), temp4)==0 )
		{
			vp.push_back(temp4);
		}
		vec_p = vp;
	}

	return vec_p;
}

//返回系数，对于a*x+b*y+c=0模型
//##ModelId=4E4A89B303D2
void DirectLine::GetCoefficient(double* coe) const
{
	coe[0] = a_;
	coe[1] = b_;
	coe[2] = c_;
}

//返回值类型是角度
//##ModelId=4E4A89B303D9
double DirectLine::GetAngle() const
{
	if ( b_ == 0)	return 90;
	else	return 180.0*atan(-a_/b_)/PI;
}

//返回斜率
//##ModelId=4E4A89B303D0
double DirectLine::GetSlot() const
{
	if ( b_ == 0)	return INT_MAX;
	else	return -a_/b_;
}



//##ModelId=4E4A89B303CA
double DirectLine::GetDistance(const DPoint2D& po) const
{
	if ( a_ == 0 && b_ ==0 )	return DBL_MAX;
	else	return fabs(a_*po.x + b_*po.y + c_)/sqrt(a_*a_ + b_*b_);
//	else	return (a_*po.x + b_*po.y + c_)/sqrt(a_*a_ + b_*b_);
}

//##ModelId=4E4A89B303DB
DirectLine& DirectLine::operator =(const DirectLine& rhs)
{
	a_ = rhs.a_;
	b_ = rhs.b_;
	c_ = rhs.c_;
	
	return *this;
	
}



//简化传感器模型a_*x+b_*x*y+c_*y+d_=0
//只能先假定y的系数不为o了，那么定义为1，求解其他值
//##ModelId=4E4A89B40036
bool Conics::Init(const std::vector<DPoint2D>& p_pts)
{
	size_t num = p_pts.size();
	if (num < 3 ) return FALSE;
	double NE[9],L[3];
	memset(NE, 0, sizeof(double)*9);
	memset(L,  0, sizeof(double)*3);

	double A[3];
	int i;
	//计算所有点的x y的最大最小值
	double xmax = p_pts[0].x;
	double xmin = p_pts[0].x;
	double ymax = p_pts[0].y;
	double ymin = p_pts[0].y;
	for( i=0;i<num;++i)
	{
		xmax = p_pts[i].x > xmax ? p_pts[i].x : xmax;
		xmin = p_pts[i].x < xmin ? p_pts[i].x : xmin;	
		ymax = p_pts[i].y > ymax ? p_pts[i].y : ymax;
		ymin = p_pts[i].y < ymin ? p_pts[i].y : ymin;		
	}


	if ( (xmax - xmin) >= (ymax - ymin) )
	{
		for( i=0;i<num;++i)
		{
			A[0]=p_pts[i].x;
			A[1]=p_pts[i].x * p_pts[i].y;
			A[2]=1;
			
			helper.pNormal(A,3,p_pts[i].y,NE,L,1.0);		
		}	
		int ll = helper.Gauss(NE, L, 3);
		a_ = L[0];
		b_ = L[1];
		c_ = -1;
		d_ = L[2];	
	}
	else
	{
		for( i=0;i<num;++i)
		{
			A[0]=p_pts[i].x * p_pts[i].y;
			A[1]=p_pts[i].y;
			A[2]=1;
			
			helper.pNormal(A,3,p_pts[i].x,NE,L,1.0);		
		}	

//		std::ofstream out("D:\\AA.txt");
//		for (int xx = 0; xx < 3; ++xx)
//		{
//			for (int yy = 0; yy <3; ++yy)
//			{
//				out << NE[xx*3+yy] << " \t";
//			}
//			out << L[xx] << "\n";
//		}
//		out.close();
		
		int ll = helper.Gauss(NE, L, 3);

		a_ = -1;
		b_ = L[0];
		c_ = L[1];
		d_ = L[2];	
	}



	return TRUE;
}

/*
double Conics::GetX(double y) const
{
	double denominator = a_ + b_ * y;
	if ( denominator == 0)    return DBL_MAX;
	else	return -(c_*y+d_)/denominator;
}

double Conics::GetY(double x) const
{
	double denominator = c_ + b_ * x;
	if ( denominator == 0)    return DBL_MAX;
	else	return -(a_*x+d_)/denominator;
}
*/
//ss是近似解   dis是严密值
//##ModelId=4E4A89B40045
double Conics::GetDistance(const DPoint2D& po) const
{
	double x = GetX( po.y);
	double y = GetY( po.x);
/*
	double k = (y - po.y)/(x - po.x);
	double cc = (a_*c_+c_*c_*k)/b_ - d_;
	double aa = b_*k;
	double bb = 2*c_*k;

	double temp = (-bb + sqrt(bb*bb - 4*aa*cc) ) / (2*aa);
	double xx;
	if ( (temp-x)*(temp-po.x) <= 0) 
		xx = temp;
	else
		xx = (-bb - sqrt(bb*bb - 4*aa*cc) ) / (2*aa);

	double yy = GetY( xx);
	double dis = sqrt(  (xx-po.x)*(xx-po.x) + (yy-po.y)*(yy-po.y) );
*/
	double ss = fabs( (x-po.x)*(y-po.y) )/sqrt( (x-po.x)*(x-po.x) + (y-po.y)*(y-po.y) );
	
	return ss;
}

//返回系数
//##ModelId=4E4A89B4004F
void Conics::GetCoefficient(double* coe) const
{
	coe[0] = a_;
	coe[1] = b_;
	coe[2] = c_;
	coe[3] = d_;
}

