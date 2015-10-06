
#ifndef _GEOBASE
#define _GEOBASE
#pragma once
#include "AttitudeAPI.h"
#include "OrbitAPI.h"
#include "UTCAPI.h"
#include "LineTimeAPI.h"
#pragma comment(lib,"CEOPDLL.lib")

#include "Safeauxilary.h"
#include "log.h"
#include "OrbitAPI.h"
#include "AttitudeAPI.h"
#include "LineTimeAPI.h"
#include "SatOrbitNew.h"
#include <iostream>
using namespace std;

#define PI 3.1415926535897932384626433832795

enum AttitdueType
{
	AttitudeAsOrbit2Body,
	AttitudeAsJ20002Body,
};

// 基准：建立参考椭球
//////////////////////
struct StrDATUM
{
	string DatumName;	// 参考椭球名称
	double a;			// 半长轴(单位：m)
	double f;			// 偏心率的倒数
	double b;			// 短半轴
	double a2;			// 半长轴的平方
	double b2;			// 半短轴的平方
	double a2_b2;		// 长短半轴比的平方
	double e2;			// 偏心率的平方
	// 重载构造函数,默认用WGS84
	StrDATUM()
	{
		DatumName = "";
		a = 6378137.0;	f = 298.257223563;	
		e2 = (2*f-1)/pow(f, 2);	// e = (2*f-1)/f
		b = (f-1)/f*a;
		a2 = a*a;	b2 = b*b;
		a2_b2 = pow(a/b, 2);
	}
	// 重载操作符=
	StrDATUM& operator=(const StrDATUM &s)
	{
		this->DatumName = s.DatumName;
		this->a = s.a;	this->f = s.f;
		this->b = s.b;	this->a2 = s.a2;
		this->b2 = s.b2;	this->a2_b2 = s.a2_b2;
		this->e2 = s.e2;
		return *this;
	}
};
//轨道的离散点信息
/*
struct StrOrbitPoint
{
	double X[6];
	int year,month,day,hour,minute;
	double second;
	double UT;
	//重载构造函数
	StrOrbitPoint()
	{
		memset(X,0,sizeof(double)*6);
		UT=second=0.0;
		year=month=day=hour=minute=0;
	}
};
*/



//指向角的离散点信息
struct StrInner
{
	  long m_Innum;	 // 垂轨向像素个数
	  StrInner()
	  {
		  m_Innum=0;
	  }
};




class GeoBase :public LogAPI
{
public:
	GeoBase(void);
	AttitudeAPI attleft, attright;
	virtual ~GeoBase(void);
public:
	// 对向量进行归一化
	void NormVector(double *R, int num);
	// 三阶向量叉乘
	void CrossMult(double *u, double *v, double *w);
	// 求矩阵相乘,A矩阵为[m,p],B矩阵为[p,n],C为[m,n] 
	void Multi(double *A,double *B,double *C ,int m,int p,int n);
	

	void transpose(double *m1, double *m2, int m, int n);
	void RotationX(double angle, double *R);
	void RotationY(double angle, double *R);
	void RotationZ(double angle, double *R);
	void rot(double fai, double omega, double kappa, double *R);
	int invers_matrix(double *m1, int n);
	void matrix2quat(double *R, double &q1, double &q2, double &q3, double &q4);
	
public:
	// 将格里高利历与UT时转化为用户规定的累计秒
	void FromYMDtoSecond(double refMJD, int year, int month, int day, int hour, 
						 int minute, double second, double& refsecond);
public:
	// 从四元数获得旋转矩阵
	void Quat2Matrix(double q1, double q2, double q3, double q4, double *R);
public:
	// 轨道拉格朗日内插
	void LagrangianInterpolation(const vector<OrbitAPI>& m_EphWGS84, double UT, int m_Ephnum, OrbitAPI &m_point);
	//姿态四元数内插
	void QuatInterpolation(const vector<AttitudeAPI>& m_Body2Orbit, double UT, int m_Attnum, AttitudeAPI &m_att);

public:
	void Rect2Geograph(StrDATUM datum, double X, double Y, double Z,double &B, double &L, double &H);
	CSatOrbit orbitHelper;
public:
	void ConvertOrbit2WGS84(double* attRot, double* j2000284rot, double* OrbitComponent, double* Body2WGS84Rot,bool  isj2000);

protected:
	const unsigned int orbitShift = 20;
	const unsigned int attshift = 20;
};
//外部封装的库函数
// 从格里高利历到约化儒略日的转化
extern "C" int  _stdcall Cal2JD(int year, int month, int day,double fracday, double *jd0, double *mjd);

#endif