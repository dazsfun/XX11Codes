#pragma once
#ifndef CURVE_H
#define CURVE_H   
#include <float.h>
#include <cmath>
#include <vector>
#include "SatOrbitNew.h"

// 应该建立成模板方式

//##ModelId=4E4A89B30290
class   DPoint2D
{
public:

	//##ModelId=4E4A89B30291
	double x;
	//##ModelId=4E4A89B30292
	double y;
public:
	//##ModelId=4E4A89B30293
	DPoint2D():x(0),y(0)
	{
	}
	//##ModelId=4E4A89B30294
	DPoint2D(double sample, double line):x(sample),y(line)
	{
	}
	//##ModelId=4E4A89B30297
	DPoint2D(const DPoint2D& rhs):x(rhs.x), y(rhs.y) {}

	//##ModelId=4E4A89B3029F
	bool operator==(const DPoint2D& rhs) const
	{
		if ( fabs(x - rhs.x) < 0.00001 && fabs(y - rhs.y) < 0.00001 )  	return true;
		return false;
	}
};


//##ModelId=4E4A89B302A2
class  __declspec(dllexport)  DPoint3D
{
public:
	//##ModelId=4E4A89B302A3
	double lat;
	//##ModelId=4E4A89B302A4
	double lon;
	//##ModelId=4E4A89B302AF
	double h;
public:
	//##ModelId=4E4A89B302B0
	DPoint3D():lat(0),lon(0),h(0)
	{
	}
	//##ModelId=4E4A89B302B1
	DPoint3D(double z):lat(0),lon(0),h(z) {}
	//##ModelId=4E4A89B302B3
	DPoint3D(double x, double y, double z):lat(x),lon(y),h(z)
	{
	}
};



//##ModelId=4E4A89B302B7
struct  __declspec(dllexport)  LineSegment {
	DPoint2D pt1, pt2;
};


//##ModelId=4E4A89B302C0
class  __declspec(dllexport)  DRect
{
public:
	//##ModelId=4E4A89B302CE
	DRect():top(0),bottom(0),left(0),right(0) {}
	//##ModelId=4E4A89B302CF
	~DRect() {}

	//##ModelId=4E4A89B302D0
	DRect(double l, double t, double r, double b):left(l), top(t),right(r) , bottom(b) {}

	//##ModelId=4E4A89B302D5
	DRect(DPoint2D LU, DPoint2D RD):left(LU.x), top(LU.y), right(RD.x), bottom(RD.y) {}

public:
	//##ModelId=4E4A89B302D8
	double top;
	//##ModelId=4E4A89B302D9
	double bottom;
	//##ModelId=4E4A89B302DE
	double left;
	//##ModelId=4E4A89B302DF
	double right;
};


//##ModelId=4E4A89B302E0
struct  __declspec(dllexport)   ConjugatePoint
{
	//##ModelId=4E4A89B302E2
	char  ID[24];
	//##ModelId=4E4A89B302EF
	DPoint2D LeftPoint;
	//##ModelId=4E4A89B302F4
	DPoint2D RightPoint;
	//##ModelId=4E4A89B302F9
	DPoint3D GroundPoint;
};


//##ModelId=4E4A89B302FD
class  __declspec(dllexport)  DPolygon
{
private:
	//##ModelId=4E4A89B3031D
	std::vector<DPoint2D> vPoints;
	//##ModelId=4E4A89B3032E
	std::vector<LineSegment> vSides;
public:
	//##ModelId=4E4A89B30332
	DPolygon();
	//##ModelId=4E4A89B30333
	DPolygon(DPoint2D* pPo, unsigned nNum);
	//##ModelId=4E4A89B3033C
	DPolygon(const std::vector<DPoint2D>& vPo );
	//##ModelId=4E4A89B3033E
	virtual ~DPolygon();
	
	//##ModelId=4E4A89B30340
	DPolygon(const DPolygon& rhs);
	//##ModelId=4E4A89B30342
	DPolygon& operator=(const DPolygon& rhs);
public:
	//##ModelId=4E4A89B3034B
	bool InPolygon(const DPoint2D& point) const;
	//##ModelId=4E4A89B3034E
	unsigned GetPolygon(DPoint2D* p) const;	
};

bool  __declspec(dllexport)  IsInRect(const DRect& WinSize,const DPoint2D& opt);
bool __declspec(dllexport) IsInRect(const DRect& WinSize,double x, double y);
double Multiply(DPoint2D p1, DPoint2D p2, DPoint2D p0);
bool IsOnline(DPoint2D point, LineSegment line);
bool Intersect(LineSegment L1, LineSegment L2);

bool  __declspec(dllexport)  IsInPoygon(const DPolygon& WinSize,const DPoint2D& opt);
unsigned PolygonIntersection(int n1,DPoint2D *p1,int n2,DPoint2D *p2,DPoint2D* p);


//////////////////////////////////////////////////////////////////////////
// 曲线类型
// 直线和双曲线
//##ModelId=4E4A89B3035B
enum CurveType
{
	CT_None, CT_DirectLine, CT_Conics,CT_AUTO
};



/////////////////////////////////////////////////////////////////////////////
// 曲线应该包含以下几个部分：
// ① 初始化，或求解过程
// ② 求解对应x和y的值，应当存在多个解的情况，怎么处理呢？ 可以给出解的范围
// ③ 求解点到曲线的距离，应当对于任何曲线都可以
/////////////////////////////////////////////////////////////////////////////
//##ModelId=4E4A89B3036B
class   Curve
{
public:
	CSatOrbit helper;
	//##ModelId=4E4A89B3036C
	virtual ~Curve()
	{}
	//##ModelId=4E4A89B3036E
	virtual bool   Init(const std::vector<DPoint2D>&  p_pts) = 0;
	//##ModelId=4E4A89B3037A
	virtual double GetX(double y) const = 0;
	//##ModelId=4E4A89B3037D
	virtual double GetY(double x) const = 0;
	//##ModelId=4E4A89B30380
	virtual double GetDistance(const DPoint2D& po) const = 0;
	//##ModelId=4E4A89B30383
	virtual void   GetCoefficient(double* coe) const = 0;
};

//##ModelId=4E4A89B30399
class  __declspec(dllexport)  DirectLine:public Curve
{
private:
	//##ModelId=4E4A89B3039B
	double a_;
	//##ModelId=4E4A89B3039C
	double b_;
	//##ModelId=4E4A89B303A9
	double c_;
public:
	//##ModelId=4E4A89B303AA
	DirectLine():a_(0),b_(0),c_(0) {}                     
	//##ModelId=4E4A89B303AB
	DirectLine( double a, double b, double c) 
	{
		a_ = a;
		b_ = b;
		c_ = c;
	}
	//##ModelId=4E4A89B303BA
	virtual ~DirectLine() {}
	
	//计算未知数即直线系数和常数 y=k*x+b (中的k和b)
	//##ModelId=4E4A89B303BC
	bool Init(const std::vector<DPoint2D>&  p_pts);  
	
	//##ModelId=4E4A89B303BE
	double DirectLine::GetY(double x) const
	{
		if ( b_ == 0)	return DBL_MAX;
		else	return (a_*x+c_)/(-b_);
	}
	
	//给定y，返回相应x
	//##ModelId=4E4A89B303C1
	double DirectLine::GetX(double y) const
	{
		if ( a_ == 0)	return DBL_MAX;
		else	return (b_*y+c_)/(-a_);
	}


	//计算点到直线的距离
	//##ModelId=4E4A89B303CA
	double GetDistance(const DPoint2D& po) const;
	
	//利用系数计算直线的两个端点,输入参数窗口大小
	//##ModelId=4E4A89B303CD
	std::vector<DPoint2D> GetEndpoint(DRect& WinSize) const;   
	
	//返回左端点坐标和直线斜率，参数为端点和窗口大小和斜率
	//##ModelId=4E4A89B303D0
	double GetSlot() const;   
	
	//返回系数，必须保证其大小能包含所有数据
	//##ModelId=4E4A89B303D2
	void GetCoefficient(double* coe) const;
	
	//返回角度值，范围(-PI/2, PI/2]
	//##ModelId=4E4A89B303D9
	double GetAngle() const;

	//重载=
	//##ModelId=4E4A89B303DB
	DirectLine& operator =(const DirectLine& rhs);
};


//简化传感器模型a_*x+b_*x*y+c_*y+d_=0
//##ModelId=4E4A89B40000
class  __declspec(dllexport)  Conics:public Curve
{
private:
	//##ModelId=4E4A89B40010
	double a_;
	//##ModelId=4E4A89B40011
	double b_;
	//##ModelId=4E4A89B4001F
	double c_;
	//##ModelId=4E4A89B40020
	double d_;
public:
	//##ModelId=4E4A89B4002E
	Conics():a_(0),b_(0),c_(0),d_(0){}
	//##ModelId=4E4A89B4002F
	Conics(double a,double b,double c,double d):a_(a),b_(b),c_(c),d_(d) {}
	//##ModelId=4E4A89B40034
	virtual ~Conics() {}
public:
	//##ModelId=4E4A89B40036
	bool Init(const std::vector<DPoint2D>&  p_pts);
	//##ModelId=4E4A89B4003F
	double GetX(double y) const
	{
		double denominator = a_ + b_ * y;
		if ( denominator == 0)    return DBL_MAX;
		else	return -(c_*y+d_)/denominator;
	}
	//##ModelId=4E4A89B40042
	double GetY(double x) const
	{
		double denominator = c_ + b_ * x;
		if ( denominator == 0)    return DBL_MAX;
		else	return -(a_*x+d_)/denominator;
	}		
	//##ModelId=4E4A89B40045
	double GetDistance(const DPoint2D& po) const;
	//##ModelId=4E4A89B4004F
	void GetCoefficient(double* coe) const;
};



#endif /* CURVE_H */
