#ifndef _CLASS_NONESTATIC_AGLLIBHELPER_H__
#define _CLASS_NONESTATIC_AGLLIBHELPER_H__

#include <string>
#include <fstream>
#include "ap.h"
using namespace alglib;
using namespace std;



class AGLLIBHELPER
{
public:
	AGLLIBHELPER(void)
	{
		_P_SingleCheckActivated = false;
	};


	///最主要的两个接口函数
	double GetXFromYValue(double y);
	double GetYFromXValue(double x);

	///新的接口函数
	///主要目的是获取在任意一点任意类型曲线的法向


	///获取最合适的内插方法
	void InitFromImgPoints(vector<ImgS> imgpoints,int squence,bool left = true);

	void DebugSingle(double x,double y);
	void CloseDebugSingle();

private:
	double _imgHeight;
	double GetHeight()
	{
		return _imgHeight;
	}
	double _imgWidth;
	double GetWidth()
	{
		return _imgWidth;
	}

	void SetHeight(double height)
	{
		_imgHeight = height;
	}

	void SetWidth(double width)
	{
		_imgWidth = width;
	}

public:
	int  GetMax(double&x,double& y)
	{
		x = _endxvalue;
		y = _endyvalue;
		return 0;
	}
	int  GetMin(double&x,double& y)
	{
		x = _startxvalue;
		y = _startyvalue;

		return -1;
	}

	int GetMaxAll(double& x,double & y)
	{
		x = _maxxvalue;
		y = _maxyvalue;
		return -1;
	}

	int GetMinAll(double & x,double & y)
	{
		x = _minxvalue;
		y = _minyvalue;
		return -2;
	}
private:
	///记录核线在影像部分内得范围
	double _maxxvalue,_minxvalue;
	double _maxyvalue,_minyvalue;
	int _valuetype;
	double _startxvalue,_startyvalue;
	double _endxvalue,_endyvalue;

private:
	bool _SingleCheck(CurveInfo curvenow);

	bool _f_GetLineType();

	bool _f_GetPOW();

protected:
	///一系列检查函数 
	///单点检查
	double _P_SinglePoint_X;
	double _P_SinglePoint_Y;
	bool _P_SingleCheckActivated;

	bool _P_F_GetNormalLine(barycentricinterpolant p,double x,DirectLine& normalline);
	TheCurve _theCurve;

	///检查工具，检查展开的正确性
	bool _P_F_CheckSpread(barycentricinterpolant p,real_1d_array a,double value);

	//求交实现函数，给出两条直线的barycentric表达式，然后给出交点，
	void _P_F_GetIntersectionPoint(const ImgS& p1,const ImgS& p2,const ImgS& p3,const ImgS& p4,ImgS& intersectionPoint);



public:
	///用于求交，也可以用于求核线的边界
	bool InterSectWithDLine(const DirectLine&,double startx,bool bdecide = false,double* intersectx=NULL,double* intersecty=NULL);

	double GetLineLength()
	{
		return sqrt(pow(_maxxvalue-_minxvalue,2) + pow(_maxyvalue-_minyvalue,2));
	}

public:
	~AGLLIBHELPER(void);
};
#endif
