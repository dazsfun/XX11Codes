#pragma once
#include "OrbitAPI.h"
#include "LineTimeAPI.h"
#include "AttitudeAPI.h"
#include <opencv2\opencv.hpp>

#pragma comment(lib,"iers2003.lib")

#ifdef __cplusplus
extern "C" {
#endif

	void _stdcall DIR_TABLES(char* file_table, int len);

	void _stdcall MJD_SOD(int* year, int* month, int* day, int* hour, int* minute, double* second, int* mjd, double* sod);

	double _stdcall DIFF_TIME(int* year, int* month, int* day, int* hour, int* minute, double* second, int* ref_year,
		int* ref_month, int* ref_day, int* ref_hour, int* ref_minute, double* ref_second);

	double _stdcall SEC2YMDHMS(int* ref_year, int* ref_month, int* ref_day, int* ref_hour, int* ref_minute,
		double* ref_second, double* seclen, int* year, int* month, int* day, int* hour, int* minute,
		double* second);

	void _stdcall UTC2TDT(int* mjd, double* sod);

	void _stdcall GPST2TDT(int* mjd, double* sod);

	void _stdcall ECEF2ECIF(int* mjd, double* sod, double* trans, double* dot_trans);
	void _stdcall ECIF2ECEF(int* mjd, double* sod, double* trans, double* dot_trans);

#ifdef __cplusplus
}
#endif


namespace Safeauxilary
{
	double InterPolation(const double& seed, const vector<double>& _seeds);
	bool GetJ2000ToWGS84(int year, int month, int day, int hour, int  minute, double second, double* rotMatrix,double* rotMatrixDot);
	//void SelfBuild(MatchToolType& _lefttool, VirtualMachine* _leftModel);
};

