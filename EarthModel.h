#pragma once
#include "SatOrbitNew.h"
#include <iostream>
#include "log.h"
using namespace std;
#include <opencv2\opencv.hpp>
using namespace cv;
#pragma comment(lib,"SatLightAndOarth.lib")
class EarthModel : public LogAPI
{
public:
	EarthModel(double avlue,double bvalue,double cvalue,double dvalue,double evalue,double fvalue,double gvalue);
	~EarthModel();

	CSatOrbit orbithelper;
public:
	DATUM* IDatum()
	{
		DATUM WGS84[] =
		{ _a, _b, _c, _d, _e, _f, _g };
		return WGS84;
	}

public:
	void Rect2Geo(const double& xvalue, const double& yvalue, const double& zvalue, double& lon, double& lat, double& height)
	{
		DATUM WGS84[] =
		{ _a, _b, _c, _d, _e, _f, _g };
		CSatOrbit orbithelper;
		double test1, test2, test3;
		orbithelper.rect2geograph(xvalue, yvalue, zvalue, WGS84, &lat, &lon , &height);
		lon = -lon / 3.1415926535897932384626433832795 * 180.0;
		lat = lat / 3.1415926535897932384626433832795 * 180.0;
	}

	void Geo2Rect(const double& lon, const double& lat, const double& height, double& xvalue, double& yvalue, double& zvalue)
	{
		DATUM WGS84[] =
		{ _a, _b, _c, _d, _e, _f, _g };
		CSatOrbit helper;
		helper.geograph2rect(lat, lon, height, WGS84, &xvalue, &yvalue, &zvalue);
	}

private:
	double _a;
	double _b;
	double _c;
	double _d;
	double _e;
	double _f;
	double _g;

public:
	bool IPosition( double* _satposition,  double* _lightvector, const double& height,  double* rectPosition,  double* earthposition);
};

