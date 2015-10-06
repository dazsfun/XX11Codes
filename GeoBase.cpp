#include "stdafx.h"
#include "GeoBase.h"

//构造/析构函数

GeoBase::GeoBase(void) :attleft(AttitudeAPI("hei")), attright(AttitudeAPI("hei")){}
GeoBase::~GeoBase(void) {}


void GeoBase::FromYMDtoSecond(double refMJD, int year, int month, int day, int hour, int minute, double second, double& refsecond)
{
	double jd0, mjd;
	Cal2JD(year, month, day, 0, &jd0, &mjd);
	refsecond = (mjd - refMJD) * 86400 + hour * 3600 + minute * 60 + second;
}

// 三阶向量叉乘
void GeoBase::CrossMult(double *u, double *v, double *w)
{
	w[0] = u[1] * v[2] - u[2] * v[1];
	w[1] = u[2] * v[0] - u[0] * v[2];
	w[2] = u[0] * v[1] - u[1] * v[0];
}

// 对向量进行归一化
void GeoBase::NormVector(double *R, int num)
{
	double retVal = 0.0;
	for (int i = 0; i < num; i++)
		retVal += pow(R[i], 2);
	retVal = sqrt(retVal);
	for (int i = 0; i < num; i++)
		R[i] /= retVal;
}

// 求矩阵相乘,A矩阵为[m,p],B矩阵为[p,n],C为[m,n] 
void GeoBase::Multi(double *A, double *B, double *C, int m, int p, int n)
{
	for (int i = 0; i < m; i++)
	for (int j = 0; j < n; j++)
	{
		double sum = 0;
		for (int k = 0; k < p; k++)
			sum = sum + A[i*p + k] * B[k*n + j];
		C[i*n + j] = sum;
	}
}


//姿态四元数内插
void GeoBase::QuatInterpolation(const vector<AttitudeAPI>& m_Body2Orbit, double UT, int m_Attnum, AttitudeAPI &m_att)
{
	//memset(&m_att, 0, sizeof(AttitudeAPI));
	m_att.utcapi = UTCAPI(UT, Enum_TimeStandard_Standard);
	if (m_Attnum < 2)
	{
		printf("m_num is less than 2, can't interpolation!\n");
		IAddLog("姿态数目过少，个数为:%d.最低允许个数为:2", m_Attnum);
		return;
	}
	// 寻找临近的两个点(对分查找)
	attleft = (m_Body2Orbit[0]);  attright = (m_Body2Orbit[0]);
	long posstart, posend, pos;
	posstart = 0, posend = m_Attnum - 1, pos = 0;
	while (posstart<posend)
	{
		pos = (posstart + posend) / 2;
		if (pos == posstart)  break;
		if ((m_Body2Orbit[pos] <= m_att) && (m_Body2Orbit[pos + 1] > m_att))
			break;
		if (m_Body2Orbit[pos] <= m_att)
			posstart = pos;
		else
			posend = pos;
	}
	if (pos < 0)	pos = 0;
	if (pos >= m_Attnum - 1)		pos = m_Attnum - 2;
	attleft = m_Body2Orbit[pos];		attright = m_Body2Orbit[pos + 1];

	// 进行内插
	///如果姿态时间过小，不允许进行内插，最大限度，根据基础文件定义进行
	double sp, sq;
	double t = (UT - attleft.utcapi.GetUTC()) / (attright.utcapi.GetUTC() - attleft.utcapi.GetUTC());
	if (abs(static_cast<int>(t)) > attshift)
	{
		cout << "姿态数据内插出错！" << endl;
		IAddLog("姿态数据内插出错，要求时间不在以后姿态所覆盖的时间范围内\n,要求时间:%lf,最小姿态时间", UT, m_Body2Orbit[0].utcapi.GetUTC());
		IAddLog("姿态时间外推上限为:%d个",attshift);
		return;
	}
	double cosa = pow(attleft.AttitudeQ0, 2) + pow(attleft.AttitudeQ1, 2) + pow(attleft.AttitudeQ2, 2) + pow(attleft.AttitudeQ3, 2);
	// 这个错误需要注意了,防止邻近两个值互为反号的情况,需要确保length>0
	if (cosa<0)
	{
		cosa = -cosa;
		attright.AttitudeQ0 = -attright.AttitudeQ0;	attright.AttitudeQ1 = -attright.AttitudeQ1;	attright.AttitudeQ2 = -attright.AttitudeQ2;	attright.AttitudeQ3 = -attright.AttitudeQ3;
	}
	if (cosa>0.9999f)
	{
		sp = 1.0 - t;	sq = t;
	}
	else
	{
		double sina = sqrt(1.0 - pow(cosa, 2));	double a = atan2(sina, cosa);	double invSina = 1.0 / sina;
		sp = sin((1.0 - t)*a)*invSina;			sq = sin(t*a)*invSina;
	}
	m_att.AttitudeQ0 = sp*attleft.AttitudeQ0 + sq*attright.AttitudeQ0;	m_att.AttitudeQ1 = sp*attleft.AttitudeQ1 + sq*attright.AttitudeQ1;
	m_att.AttitudeQ2 = sp*attleft.AttitudeQ2 + sq*attright.AttitudeQ2;	m_att.AttitudeQ3 = sp*attleft.AttitudeQ3 + sq*attright.AttitudeQ3;
	double matrixtemp[9];
	Quat2Matrix(m_att.AttitudeQ0, m_att.AttitudeQ1, m_att.AttitudeQ2, m_att.AttitudeQ3, matrixtemp);
	orbitHelper.rot2eulor(matrixtemp, m_att.AttitudeROLL, m_att.AttitudePITCH, m_att.AttitudeYAW);
}

// 功能：从四元数获得旋转矩阵
//////////////////////////////////////
void GeoBase::Quat2Matrix(double q1, double q2, double q3, double q4, double *R)
{
	// 求取欧式空间长度
	double length2 = pow(q1, 2) + pow(q2, 2) + pow(q3, 2) + pow(q4, 2);
	// 非单位四元数
	if (fabs(length2 - 1.0) >= DBL_MIN * 10)
	{
		if (fabs(length2) <= DBL_MIN * 2)
		{
			printf("Quat2Matrix Error!\n");	// 输入不是单位四元数,并且欧式距离太短
			memset(R, 0, sizeof(double)* 9);
			return;
		}
		double length = sqrt(length2);
		q1 /= length;  q2 /= length;  q3 /= length;  q4 /= length;
	}
	// 计算旋转矩阵
	R[0] = 1.0 - 2.0*(q2*q2 + q3*q3);	R[1] = 2.0 * (q1*q2 + q3*q4);	R[2] = 2.0 * (q1*q3 - q2*q4);
	R[3] = 2.0 * (q1*q2 - q3*q4);	R[4] = 1.0 - 2.0*(q1*q1 + q3*q3);	R[5] = 2.0 * (q2*q3 + q1*q4);
	R[6] = 2.0 * (q1*q3 + q2*q4);	R[7] = 2.0 * (q2*q3 - q1*q4);	R[8] = 1.0 - 2.0*(q1*q1 + q2*q2);
}

//轨道拉格朗日内插
void GeoBase::LagrangianInterpolation(const vector<OrbitAPI>& m_EphWGS84, double UT, int m_Ephnum, OrbitAPI &m_point)
{
	int order = 7;
	UTCAPI HEI(UT, Enum_TimeStandard_Standard);
	m_point.utcapi = HEI;
	double up = 1, down = 1;
	long posstart, posend, pos;
	if (m_EphWGS84.size() <= 7)
	{
		cout << "轨道数据内插错误" << endl;
		IAddLog("轨道数据个数过少 ，不满足内插要求，轨道数据内插错误，实际轨道个数为:%d",m_EphWGS84.size());
		return;
	}
	posstart = 0, posend = m_Ephnum - 1, pos = 0;
	while (posstart < posend)
	{
		pos = (posstart + posend) / 2;
		if (pos == posstart)  break;
		if ((m_EphWGS84[pos] <= m_point) && (m_point <= m_EphWGS84[pos + 1])) break;
		if (m_EphWGS84[pos] < m_point)
			posstart = pos;
		else
			posend = pos;
	}
	if (pos - order / 2 < 0)   posstart = 0;
	else                posstart = pos - order / 2;
	if (pos + order / 2 >= m_Ephnum - 1) posend = m_Ephnum - 1;
	else                     posend = pos + order / 2;
	int i, j;
	int shiftleft =fabs( (UT - m_EphWGS84[posstart].utcapi.GetUTC()) / (m_EphWGS84[posstart+1].utcapi.GetUTC() - m_EphWGS84[posstart].utcapi.GetUTC()));
	int shiftright = fabs((UT - m_EphWGS84[posend].utcapi.GetUTC()) / (m_EphWGS84[posend].utcapi.GetUTC() - m_EphWGS84[posend -1].utcapi.GetUTC()));
	if (min(shiftleft, shiftright) > orbitShift)
	{
		cout << "轨道数据内插出错" << endl;
		IAddLog("轨道数据内插出错，要求时间不在轨道时间所覆盖区域内\n,要求时间：%lf，最小时间：%lf",UT,m_EphWGS84[0].utcapi.GetUTC());
		IAddLog("最大允许的轨道外推格式为:%d",orbitShift);
		return;
	}
	/*
				if(i!=j)
				{
				up *=(UT-m_EphWGS84[i].UT);
				down *=(m_EphWGS84[j].UT-m_EphWGS84[i].UT);
				}
				*/
	m_point.GPSX = m_point.GPSXV = m_point.GPSY = m_point.GPSYV = m_point.GPSZ = m_point.GPSZV = 0;
	for (j = posstart; j <= posend; j++)
	{
		up = 1, down = 1;
		for (i = posstart; i <= posend; i++)
		if (i != j)
		{
			double temp = m_EphWGS84[i].utcapi.GetUTC();
			up *= (UT - m_EphWGS84[i].utcapi.GetUTC());
			down *= (m_EphWGS84[j].utcapi.GetUTC() - m_EphWGS84[i].utcapi.GetUTC());
		}
		/*
		for(i=0;i<6;i++)
		m_point.X[i] += m_EphWGS84[j].X[i]*up/down;
		*/
		m_point.GPSXV += m_EphWGS84[j].GPSXV * up / down;
		m_point.GPSYV += m_EphWGS84[j].GPSYV * up / down;
		m_point.GPSZV += m_EphWGS84[j].GPSZV * up / down;
		m_point.GPSX += m_EphWGS84[j].GPSX * up / down;
		m_point.GPSY += m_EphWGS84[j].GPSY * up / down;
		m_point.GPSZ += m_EphWGS84[j].GPSZ * up / down;
	}
}

void GeoBase::rot(double fai, double omega, double kappa, double *R)
{
	memset(R, 0, 9 * sizeof(double));
	double RX[9], RY[9], RZ[9], t[9];
	RotationX(omega, RX);
	RotationX(fai, RY);
	RotationX(kappa, RZ);
	Multi(RY, RX, t, 3, 3, 3);
	Multi(t, RZ, R, 3, 3, 3);
}

void GeoBase::Rect2Geograph(StrDATUM datum, double X, double Y, double Z, double &B, double &L, double &H)
{
	// 对于球的特殊处理
	if (datum.f == 0)
	{
		// 对于椭球情况的处理
		double XY = sqrt(X*X + Y*Y);				// XY投影面上的长度
		// 先求出经度
		L = atan2(Y, X);
		// 再求出纬度
		B = atan2(Z, XY);
		// 再求出高程
		H = sqrt(XY*XY + Z*Z) - datum.a;
		return;
	}
	// 对于椭球情况的处理
	double XY = sqrt(X*X + Y*Y);				// XY投影面上的长度
	// 先求出经度
	L = atan2(Y, X);
	////////////////////////////////////////
	// 接着求纬度
	////////////////////////////////////////
	double z1, z2, dz, B0, B1, B2, delta;
	B0 = atan2(Z, XY);							// 地心纬度
	delta = 0.0001;
	int num = 0;
	do
	{
		B1 = B0;		B2 = B0 - delta;
		z1 = XY*datum.a2_b2*tan(B1) - (datum.a2_b2 - 1)*sin(B1)*datum.b / sqrt(1 - datum.e2*pow(cos(B1), 2));
		z2 = XY*datum.a2_b2*tan(B2) - (datum.a2_b2 - 1)*sin(B2)*datum.b / sqrt(1 - datum.e2*pow(cos(B2), 2));
		dz = z1 - z2;		B0 += (Z - z1) / dz*delta;
		if (num++ > 1000) break;
	} while (fabs(Z - z1) > fabs(dz));
	B = atan(tan(B0)*datum.a2_b2);
	////////////////////////////////////////
	// 最后求高程
	////////////////////////////////////////
	double re = datum.b / sqrt(1.0 - datum.e2*pow(cos(B0), 2));
	H = sqrt(pow(Z - re*sin(B0), 2) + pow(XY - re*cos(B0), 2));
	if ((XY*XY + Z*Z) < pow(re, 2))
		H *= -1.0;
}

void GeoBase::transpose(double *m1, double *m2, int m, int n)
{
	int i, j;

	for (i = 0; i < m; i++)
	for (j = 0; j < n; j++)
		m2[j*m + i] = m1[i*n + j];
}

void GeoBase::matrix2quat(double *R, double &q1, double &q2, double &q3, double &q4)
{
	double q1q2q3 = sqrt((R[1] + R[3])*(R[2] + R[6])*(R[5] + R[7]));
	q1 = q1q2q3 / (R[5] + R[7]) / 2;
	q2 = q1q2q3 / (R[2] + R[6]) / 2;
	q3 = q1q2q3 / (R[1] + R[3]) / 2;
	q4 = sqrt((R[0] + R[4] + R[8] + 1) / 4);

	if (fabs(2 * (q1*q2 + q3*q4) - R[1]) > 0.001)
		q4 = -q4;

}
void GeoBase::RotationX(double angle, double *R)
{
	memset(R, 0, 9 * sizeof(double));
	R[0] = 1;
	R[4] = cos(angle); R[5] = -sin(angle);
	R[7] = sin(angle); R[8] = cos(angle);
}

void GeoBase::RotationY(double angle, double *R)
{
	memset(R, 0, 9 * sizeof(double));
	R[0] = cos(angle);        R[2] = -sin(angle);
	R[4] = 1;
	R[6] = sin(angle);        R[8] = cos(angle);
}
void GeoBase::RotationZ(double angle, double *R)
{
	memset(R, 0, 9 * sizeof(double));
	R[0] = cos(angle); R[1] = -sin(angle);
	R[3] = sin(angle); R[4] = cos(angle);
	R[8] = 1;
}

int GeoBase::invers_matrix(double *m1, int n)
{
	int *is, *js;

	int i, j, k, l, u, v;

	double temp, max_v;

	is = (int *)malloc(n*sizeof(int));

	js = (int *)malloc(n*sizeof(int));

	if (is == NULL || js == NULL)
	{

		printf("out of memory!\n");

		return(0);

	}

	for (k = 0; k < n; k++){
		max_v = 0.0;
		for (i = k; i < n; i++)
		for (j = k; j<n; j++){
			temp = fabs(m1[i*n + j]);
			if (temp>max_v){
				max_v = temp; is[k] = i; js[k] = j;
			}
		}
		if (max_v == 0.0){
			free(is); free(js);
			printf("invers is not availble!\n");
			return(0);
		}
		if (is[k] != k)
		for (j = 0; j < n; j++){
			u = k*n + j; v = is[k] * n + j;
			temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
		}
		if (js[k] != k)
		for (i = 0; i < n; i++){
			u = i*n + k; v = i*n + js[k];
			temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
		}
		l = k*n + k;
		m1[l] = 1.0 / m1[l];
		for (j = 0; j < n; j++)
		if (j != k){
			u = k*n + j;
			m1[u] *= m1[l];
		}
		for (i = 0; i < n; i++)
		if (i != k)
		for (j = 0; j < n; j++)
		if (j != k){
			u = i*n + j;
			m1[u] -= m1[i*n + k] * m1[k*n + j];
		}
		for (i = 0; i < n; i++)
		if (i != k){
			u = i*n + k;
			m1[u] *= -m1[l];
		}
	}
	for (k = n - 1; k >= 0; k--){
		if (js[k] != k)
		for (j = 0; j < n; j++){
			u = k*n + j; v = js[k] * n + j;
			temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
		}
		if (is[k] != k)
		for (i = 0; i < n; i++){
			u = i*n + k; v = i*n + is[k];
			temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
		}
	}
	free(is); free(js);
	return(1);
}

void GeoBase::ConvertOrbit2WGS84(double* attRot, double* j2000284rot, double* position, double* Body2WGS84Rot, bool  isj2000)
{
	if (!isj2000)
	{
		double XYZs[3]; double XYZvs[3];
		for (int i = 0; i < 3; i++)
		{
			XYZs[i] = position[i];
			XYZvs[i] = position[3 + i];
		}
		double Orbit2WGS84[9] = { 0 };
		double X2[3], Y2[3], Z2[3];
		NormVector(XYZs, 3);
		Z2[0] = -XYZs[0];
		Z2[1] = -XYZs[1];
		Z2[2] = -XYZs[2];
		CrossMult(Z2, XYZvs, Y2);
		NormVector(Y2, 3);
		CrossMult(Y2, Z2, X2);
		NormVector(X2, 3);
		Orbit2WGS84[0] = X2[0]; Orbit2WGS84[1] = Y2[0]; Orbit2WGS84[2] = Z2[0];
		Orbit2WGS84[3] = X2[1]; Orbit2WGS84[4] = Y2[1]; Orbit2WGS84[5] = Z2[1];
		Orbit2WGS84[6] = X2[2]; Orbit2WGS84[7] = Y2[2]; Orbit2WGS84[8] = Z2[2];
		//invers_matrix(Orbit2WGS84, 3);
		Multi(Orbit2WGS84, attRot, Body2WGS84Rot, 3, 3, 3);
	}
	else
	{
		invers_matrix(attRot, 3);
		/*
		cout << "attis" << endl;
		for (int i = 0; i < 9; i++)
			cout << attRot[i] << ",";
		cout << endl;
		cout << "after invers attis" << endl;
		invers_matrix(attRot, 3);
		for (int i = 0; i < 9; i++)
			cout << attRot[i] << ",";
		cout << endl;*/
		invers_matrix(j2000284rot, 3);
		Multi(j2000284rot, attRot, Body2WGS84Rot, 3, 3, 3);
		//invers_matrix(Body2WGS84Rot, 3);
	}
}
