
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

// ��׼�������ο�����
//////////////////////
struct StrDATUM
{
	string DatumName;	// �ο���������
	double a;			// �볤��(��λ��m)
	double f;			// ƫ���ʵĵ���
	double b;			// �̰���
	double a2;			// �볤���ƽ��
	double b2;			// ������ƽ��
	double a2_b2;		// ���̰���ȵ�ƽ��
	double e2;			// ƫ���ʵ�ƽ��
	// ���ع��캯��,Ĭ����WGS84
	StrDATUM()
	{
		DatumName = "";
		a = 6378137.0;	f = 298.257223563;	
		e2 = (2*f-1)/pow(f, 2);	// e = (2*f-1)/f
		b = (f-1)/f*a;
		a2 = a*a;	b2 = b*b;
		a2_b2 = pow(a/b, 2);
	}
	// ���ز�����=
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
//�������ɢ����Ϣ
/*
struct StrOrbitPoint
{
	double X[6];
	int year,month,day,hour,minute;
	double second;
	double UT;
	//���ع��캯��
	StrOrbitPoint()
	{
		memset(X,0,sizeof(double)*6);
		UT=second=0.0;
		year=month=day=hour=minute=0;
	}
};
*/



//ָ��ǵ���ɢ����Ϣ
struct StrInner
{
	  long m_Innum;	 // ���������ظ���
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
	// ���������й�һ��
	void NormVector(double *R, int num);
	// �����������
	void CrossMult(double *u, double *v, double *w);
	// ��������,A����Ϊ[m,p],B����Ϊ[p,n],CΪ[m,n] 
	void Multi(double *A,double *B,double *C ,int m,int p,int n);
	

	void transpose(double *m1, double *m2, int m, int n);
	void RotationX(double angle, double *R);
	void RotationY(double angle, double *R);
	void RotationZ(double angle, double *R);
	void rot(double fai, double omega, double kappa, double *R);
	int invers_matrix(double *m1, int n);
	void matrix2quat(double *R, double &q1, double &q2, double &q3, double &q4);
	
public:
	// �������������UTʱת��Ϊ�û��涨���ۼ���
	void FromYMDtoSecond(double refMJD, int year, int month, int day, int hour, 
						 int minute, double second, double& refsecond);
public:
	// ����Ԫ�������ת����
	void Quat2Matrix(double q1, double q2, double q3, double q4, double *R);
public:
	// ������������ڲ�
	void LagrangianInterpolation(const vector<OrbitAPI>& m_EphWGS84, double UT, int m_Ephnum, OrbitAPI &m_point);
	//��̬��Ԫ���ڲ�
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
//�ⲿ��װ�Ŀ⺯��
// �Ӹ����������Լ�������յ�ת��
extern "C" int  _stdcall Cal2JD(int year, int month, int day,double fracday, double *jd0, double *mjd);

#endif