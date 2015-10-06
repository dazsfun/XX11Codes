// SatOrbit.h: interface for the CSatOrbit class.
//
//////////////////////////////////////////////////////////////////////
//#include "SatOrbit.h"
#pragma once
#include "GeoRSDataStruct.h"
#include <string>
#include <vector>
using namespace std;
#ifndef ClassInDLL_H
#define ClassInDLL_H
#ifdef _CLASSINDLL
#define CLASSINDLL_CLASS_DECL __declspec(dllexport)
#else
#define CLASSINDLL_CLASS_DECL __declspec(dllimport)
#endif
#endif // ClassInDLL_H



typedef double Quat[4];
class  CLASSINDLL_CLASS_DECL CSatOrbit
{
public:
	//直角坐标和经纬度坐标之间的换算关系
	static void geograph2rect(const double& lat, const double& lon,const double& zn, 
						  DATUM *datum, double *x, double *y, double *z);
	static const bool rect2geograph(const double& x,const double& y, const double& z, DATUM *datum,
						  double *lat, double *lon, double *zn);
	static void wgs842geoPseudo(const double& x,const double& y,const double& z,double *psi, double *lambda, double *geoDisNadir);
	static void geoPseudo2wgs84(const double& psi,const double& lambda,const double& geoDisNadir,
						  double *x, double *y, double *z);
	//地心纬度变参心纬度
	static const double lat2psi(double lat, DATUM *datum);
	//wgs84到极坐标 km
	void wgs842geoPseudo(double x, double y, double z,double *psi, double *lambda, double *geoDisNadir);
	void geoPseudo2wgs84(double psi, double lambda, double geoDisNadir,double *x, double *y, double *z);
	static const double MapXY2LatLon(double X1,double Y1,double& lon,double& lat);
private:
	
	static const double lat2psi(const double& lat, DATUM *datum);
	static const double Z_XbarPsi(const double& xBar, const double& psi, DATUM *datum);
public:
	///如果不希望函数的参数在计算过程中被改变（对象（指针所指的对象），指针），就应该加上限定
	///此外，尽量的传引用而不是值
	///该类的函数当然都应该是const和static的，声明一下也不费事，避免日后代码升级时本类发生其不应当具有的一些变化
	static const double distance(const double& x1,const double& y1,const double& x2,const double &y2);
	//static const double singlepixelsize(const GCPSTRUCT& gcp1,const GCPSTRUCT& gcp2);
	static const void norn(double* result);
	static const double pointmult(const double *a,const double *b);
	static const void ephemeris2OrbitEle(const POINT_Ep& ep,ORBIT& orbit);
	static const void OrbitEle2ephemeris_M(const ORBIT& orbit,POINT_Ep &ep);
	static const void OrbitEle2ephemeris_E(const ORBIT& orbit,POINT_Ep &ep);

	//计算偏流角
	static const void FromORBIT2biasangle(const POINT_Ep& ep,double &biasangle);

	static const void mult(const double*m1,const double* m2,double* result,const int& i_1,const int& j_12,const int j_2);

	static const void pNormal(const double* a,const int& n,const double& b,double* aa,double* ab,const double p);

	static const int Gauss(double *A,double *b,int n);

	static const void FromtotalSecond2hourminutesecond(const double& totalsecond,int &hour,int &minute,double &second); 

	static const void normalvect(const double *x,double* y);

	static const void crossmultnorm(const double* x,const double* y,double* z);

	static const void transpose(const double* m1,double *m2,const int&m,const int& n);

	static const int invers_matrix(double* const m1,const int& n);

	static const void GaussExt(double* const aa, double* ab,double* const x,const int&n);

	static const void GaussL(double* const aa,  double* const ab,double* const x,const int& n,const double &btb);

	static const void RotationX(const double& angle,double* const R);

	static const void GaussExtQ(double *aa,double *ab,double *x,int n,double *Q);

	static const void RotationY(const double& angle,double* const R);

	static const void RotationZ(const double& angle,double* const R);

	static const void PredictbasedAffine(const RPCMODEL &rpcmodel,const double& x0,const double& y0,const double &h0,double &P0,double &L0);

	//公历和儒略日之间的关系
	static const double JulDate(struct ut_instant* const date);
	static const long CalDate(struct ut_instant* const date);

	//JD为UTC时间
	static const void JD2ChristianYear(const double& JD,int& year,int &month,int& day,int& hour,int& minute,double& second);
	static const void ChristianYear2JD(const int& year,const int& month,const int& day,const int& hour,const int& minute,const double second,double& JD);

	//
	static const void prbub(GCPSTRUCT* const p, const int& n);
	static const BOOL larger(const GCPSTRUCT& m,const GCPSTRUCT& n);
	static const void rsplit(GCPSTRUCT* const p,const int& n,int* const m);
	static const void prqck(GCPSTRUCT* const p,const int& n);

	static const void statistic(GCPSTRUCT * const p,const int & num,double &avx,double &avy,double &avxy,double &dx,double &dy,double &dxy);
	static const void statistic(Attitude * const p,const int& num,double &avroll,double &avpitch,double &avyaw,double &droll,double &dpitch,double &dyaw);
	static const void statistic(GCPSTRUCTEX* const p,const int & num,
		double &avx,double &avy,double &avxy,double &avsita,
		double &avrx,double &avry,double &avrxy,double &avrsita,
		double &dx,double &dy,double &dxy,double &dsita,
		double &drx,double &dry,double &drxy,double &drsita);

	static const void statistic(GCPSTRUCTEXEX * const p,const int& num,double &avC,double &avrC,double &avsi,
		double &avCx,double &avrCx,double &avCy,double &avrCy,
		double &dCx,double &drCx,
		double &dCy,double &drCy,
		double &dC,double &drC,
		double &dsi,
		double &avsita,double &dsita,
		double &avrsita,double &drsita);

	static const void ConvertTimeString2TimeString(const char* str,CString& timestring);

	static const void ConvertTimeString2YearMonthDayHourMinuteSecond(CString timestring,int &year,int &month,int &day,int &hour,int &minute,double &second);

	static int* const Combination(const int& n,const int& r,int&k); 

	static const double lag(const double& t,const double& t0,const double& dt,const double& e0,const double& e1,const double& e2,const double& e3);

	static const void Rotate3D(const int& direction,const double& theta,double a[3][3]);

	static const void MatrixMulitVec(double a[3][3],double f[3]);
	static const void MatrixMulitMatrix(double a[3][3],double b[3][3]);
	static const double ehpP(const double& t,const double& t0,const double& dt,const double& ep0,const double& ep1,const double& ev0, const double& ev1); 
	static const double ehpV(const double& t,const double& t0,const double& dt,const double& ep0,const double& ep1,const double& ev0, const double& ev1);
	static const int multiplyQ(double qbe[4],double qaf[4]);
	static const int interpolationQ(const double& t,const double& t0,const double& t1,double qpre[4],double qpost[4]);
	static const bool MatrixInvert(double a[3][3],double ainv[3][3]);
	static const double MatrixDeterminant(double a[3][3]);
	//向量操作
	static const bool Vector3Normalize(double a[3]);
	static const double Vector3Length(double a[3]);
	static const bool Vector3CrossProduct(double a[3],double b[3],double result[3]);
	static const double Vector3InnerProduct(double a[3],double b[3]);
	//高斯消元法
	static const int guass_metrix(const int& n,double a[4][4],double b[4],const double& eps);
	static const double AngleFromsinAndcos(const double& sinangle,const double& cosangle);
	static const void LagrangianInterpolation(POINT_Ep * const Ep_Array,const int& epNum,const double& CurUT,POINT_Ep &point);
	static const void crossmult(double * const a,double * const b,double * const result);
	static const double dabs(const double& a);
	static const void solve33(double * const a1,double* const A);
	static const void GetMaxMin(double * const xxxx,const int& num,double& xmax,double &xmin);
	static const void GetMaxMinVDx(double * const xxxxx,const int & num,double& xmax,double& xmin,double &v,double& dx);
	//static const void RPCcompensation(const RPCMODEL& rpc,double & sample,double & line,const double& H);
	static const double RPCPLH(const double *  rpcof,const double& P,const double& L,const double& H);
	static const double RPCPLHDynamic(const double* rpcof, const double& P, const double& L, const double& H, int * mark);
	static const void GetMaxMinRMS(double * const xxxxx,const int& num,double & xmax,double& xmin,double& RMS);
	static const void LineSample2LATLONGHEIGHT(const double& line,const double& sample,const double& height,const RPCMODEL& rpcmodel,double& latitude,double &longitude);
	static const void LATLONGHEIGHT2LineSample(RPCMODEL rpcmodel,const double& latitude,const double& longitude,const double& height,double& x,double& y);
	static const BOOL GetRPCMODEL(CString rpcname,RPCMODEL& rpcmodel);
	static const BOOL GetAffineMODEL(CString affinename,RPCMODEL & rpcmodel);
	static const double SpaceIntersection(const long& modelnum,int* const modelid,RPCMODEL* const rpcmodel,double* const x,double* const y,double& lat,double& lon,double& height); 
	static const void Compute_avAnddx(const double*a,const int& num,double& av,double& dx);
	static const double PartialDerivativeRPC2P(const double *rpccof,const double& P,const double& L,const double& H);
	static const double PartialDerivativeRPC2L(const double *rpccof,const double& P,const double& L,const double& H);
	static const double PartialDerivativeRPC2H(const double *rpccof,const double& P,const double& L,const double& H);

	static const void PartialDerivativesampleline2latlonheight(const double& P,const double &L,const double& H,const RPCMODEL & rpcmodel,
		double &s2lat,double &s2lon,double &s2height,
		double &l2lat,double &l2lon,double &l2height,double &s,double &l);


	static const void rot(const double& fai,const double& omega,const double& kappa,double * const R) ;

	static const void crossmult(const double& ax,const double& ay,const double& az,const double& bx,const double& by,const double& bz,double * const result);

	static const void rotx(const double& angle,double * const R);

	static const void roty(const double& angle,double * const R);

	static const void rotz(const double& angle,double * const R);

	static const void quat2matrix(double m[9], const stQuaternion& quat);

	static void const quat2matrix(const double& q1,const double& q2,const double& q3,const double& q4,double * const R);

	static void const jcbj(double * const a,const int& n,const double& keps);

	static const void diag(double * const a, const int& m,const int& n,double * const b);

	static const void  trmul(double * consta , double * const b, const int& m,const  int& n,const int& k,double * const c);

	static const double vectornorm(double * const a,const int& m);

	static const void lamda2etarho(const double& lambda,double * const ata,double * const atl,const int& n,const double& ltl,double &eta,double &rho);

	static const void  lamda2detadrho(const double& lambda,double * const ata,double * const atl,const int& n,const double& ltl,double &deta,double &drho);

	static const void lamda2ddetaddrho(const double& lambda,double * const ata,double * const atl,const int& n, const double& ltl,double &ddeta,double &ddrho);

	static const void lcfun(double * const lambda,const int& nump,double * const ata,double * const atl,const int& n ,const double& ltl,double *g);

	static const void  minvector(double * const vector,const int& m,int &mini);

	static const double fminbnd(const double& reg_min,const double& reg_max,double * const ata,double * const atl,const int& n ,const double& ltl,const double& eps);

	static const double L_curve(double * const ata,double * const atl,const int& n,const double& ltl);

	static const void crossmultnorn(const double& ax,const double& ay,const double& az,const double& bx,const double& by,const double& bz,double * const result);

	static const double pointmult(const double& ax,const double& ay,const double& az,const double& bx,const double& by,const double& bz);

	static const double  Etof(const double& e,const double& E);

	static const double  ftoE(const double& e,const double& E);

	static const double EtoM(const double& e,const double& E);

	static const double MtoE(const double& e,const double& M);

	static const void OrbitEle2ephemeris_f(const ORBIT& orbit,POINT_Ep &ep);

	static const void OrbitEle2ephemeris_t0(const ORBIT& orbit,POINT_Ep &ep);

	static const int findlm(double *  A,const int& ic,const int& n,double * const maxx);

	static const void exlij(double *  A,double *  b,const int& li,const int& lj,const int& n);

	static const void CSatOrbit::eliminate(double *  A,double *  b,const int& ii,const int& n);

	//static const void exlij(double * const A,double * const b,const int& li,const int& lj,const int& n);

	static const void GaussExtQ(double * const aa,double * const ab,double * const x,const int& n,double * const Q);

	static const __int64 factorial(const int& n);

//	static int * const Combination(const int& n,const int& r,int &k);

//	static const void PredictbasedAffine(const RPCMODEL& rpcmodel,const double& x0,const double& y0,const double& h0,double &P0,double &L0);

	static const void PreciseBasedAffine(const RPCMODEL& rpcmodel,double &latitude,double &longitude,const double& dlat,const double& dlon,const double& x0,const double& y0,const double& h0);

	static const double distanceOnePixel(const RPCMODEL& rpcmodel,const double& latitude,const double& longitude);

	static const void RPCSpaceIntersectionErrorEquation(const double& x,const double& y,const double& lat,const double& lon,const double& height,const RPCMODEL& rpcmodel,double * const bx,double * const by,double &lx,double &ly);

	static const void RPCBlockAdjustmentErrorEquation(const double& x,const double& y,const double& lat,const double& lon,const double& height,const RPCMODEL& rpcmodel,double * const ax,double * const ay,double * const bx,double * const by,double &lx,double &ly);

	static const void RPCSpaceResectionErrorEquation(const double& x,const double& y,const double& lat,const double& lon,const double& height,const RPCMODEL& rpcmodel,double * const ax,double *const ay,double &lx,double &ly);

	static const void matrix2quat(double * const R,double &q1,double &q2,double &q3,double &q4);
	//static void RPCSpaceResectionErrorEquation(const double& x,const double& y,const double& lat,const  double& lon,const double& height,const RPCMODEL& rpcmodel,double * const ax,double * const ay,double &lx,double &ly);

	static const double delta;

	static const BOOL rot2eulor(double * const R,double &Phi,double &Omega,double &Kappa);
	// static const void RPCBlockAdjustmentErrorEquation(const double& x,const double& y,const double& lat,const double& lon,const double& height,const RPCMODEL& rpcmodel,double * const ax,double * const ay,double * const bx,double * const by,double &lx,double &ly);
//	static const void PreciseBasedAffine(const RPCMODEL& rpcmodel,double &latitude,double &longitude,const double& dlat,const double& dlon,const double& x0,const double& y0,const double& h0);
//	static const void OrbitEle2ephemeris_f(const ORBIT& orbit,POINT_Ep &ep);

	//static const void OrbitEle2ephemeris_t0(const ORBIT& orbit,POINT_Ep &ep)

	static const int TruncateEphemerisForInterpolation(std::vector<POINT_Ep>& vEp, const double& dStarScanTime,const  double& dEndScanTime);

	static const int TruncateAttitudeForInterpolation(std::vector<Attitude>& vAtt,const double& dStarScanTime, const double& dEndScanTime);

	static const bool Compute_OrbitAttitudeModel_Regression(std::vector<Attitude> vAtt,std::vector<POINT_Ep> vEp, const double& dStarScanTime,const double& dEndScanTime,
										   const double& dCenterTime, OrbitAttitudeRegressionModel& Rmdl,const int & nEpOrder,const int& nAttOrder);

	static const void PolynominalFitting(double *x,double *y,double *p,const int& n,const int& order);

	static const void PolynominalFittingError(double *x,double *y,double *p,int n,const int& order,double *error);

	static const double PolyValue(const double *p,const double& UT);

	static const CString GetProcessPath();

	static const POINT_Ep EpLagrangianInterpolation(const std::vector<POINT_Ep>& vEp, const double& CurUT);

	static const void AttQuaternionInterpolation(const std::vector<Attitude>& vAtt, const double& CurUT, double* R);

	static const POINT_Ep EpPolyInterpolation(const OrbitAttitudeRegressionModel&  OARM, const double& CenterTimeUT, const double& CurUT);

	static const void quat_slerp(Quat c, Quat a, Quat b, const double& t);
};

//void bisectroot(double a, double b, double h, int n, double eps, int* m, double x[4]);

//double __declspec(dllexport)  SpaceIntersection(long modelnum,int *modelid,RPCMODEL *rpcmodel,double *x,double *y,double &lat,double &lon,double &height);

//void RPCSpaceIntersectionErrorEquation(double x,double y,double lat,double lon,double height,RPCMODEL rpcmodel,double *bx,double *by,double &lx,double &ly);







