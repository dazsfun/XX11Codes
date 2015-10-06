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


	///����Ҫ�������ӿں���
	double GetXFromYValue(double y);
	double GetYFromXValue(double x);

	///�µĽӿں���
	///��ҪĿ���ǻ�ȡ������һ�������������ߵķ���


	///��ȡ����ʵ��ڲ巽��
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
	///��¼������Ӱ�񲿷��ڵ÷�Χ
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
	///һϵ�м�麯�� 
	///������
	double _P_SinglePoint_X;
	double _P_SinglePoint_Y;
	bool _P_SingleCheckActivated;

	bool _P_F_GetNormalLine(barycentricinterpolant p,double x,DirectLine& normalline);
	TheCurve _theCurve;

	///��鹤�ߣ����չ������ȷ��
	bool _P_F_CheckSpread(barycentricinterpolant p,real_1d_array a,double value);

	//��ʵ�ֺ�������������ֱ�ߵ�barycentric���ʽ��Ȼ��������㣬
	void _P_F_GetIntersectionPoint(const ImgS& p1,const ImgS& p2,const ImgS& p3,const ImgS& p4,ImgS& intersectionPoint);



public:
	///�����󽻣�Ҳ������������ߵı߽�
	bool InterSectWithDLine(const DirectLine&,double startx,bool bdecide = false,double* intersectx=NULL,double* intersecty=NULL);

	double GetLineLength()
	{
		return sqrt(pow(_maxxvalue-_minxvalue,2) + pow(_maxyvalue-_minyvalue,2));
	}

public:
	~AGLLIBHELPER(void);
};
#endif
