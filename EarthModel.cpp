#include "stdafx.h"
#include "EarthModel.h"
#include <math.h>

EarthModel::EarthModel(double avalue, double bvalue,double cvalue,double dvalue,double evalue,
	double fvalue, double gvalue) :_a(avalue), _b(bvalue), _c(cvalue), _d(dvalue), _e(evalue), _f(fvalue),
	_g(gvalue)
{

}


EarthModel::~EarthModel()
{
}

///当把原来一段100行的代码分散成4个60行左右的函数
///是为了增强可读性，并提高查错的效率，当然并不是针对每个长函数都这么做
bool EarthModel::IPosition( double* _satposition,  double* lightvector, const double& height,  double* rectPosition,   double* earthposition)
{

	double A = _a + height;
	double B = _b + height;

	double sumtest = 0;
	for (int i = 0; i < 3; i++)
	{
		//cout << lightvector[i] << "\t";
		sumtest += lightvector[i] * lightvector[i];
	}

	//cout << endl;
	sumtest = sqrt(sumtest);

	if (sumtest != 1)
	{
		for (int i = 0; i < 3; i++)
			lightvector[i] = lightvector[i] / sumtest;
	}



	double rotcompenst[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 }; //{ 0.6767922681376681, -0.5409248699687399, 0.4993520910523681,
	//0.7242393318599709, 0.6108741491698663, -0.3198596005484512,
	//-0.1320212709228849, 0.5781289293274413, 0.8051939673759436 };

	double _lightvector[3];


	double roll;
	
	double pitch;
	double yaw;
	roll = pitch = yaw = 0;
	roll = -1.645487798741025047 - 0.27187423317486969 + 0.0014685925645 + 0.01186503466776 + 0.049231786731971233
		+ 0.0034971484082110659 - 0.00033083061214364377
		+ 0.0061503490127519673
		- 0.0095721650392792818
		+ 0.0048305350932080701
		+ 0.013334701983824913
		-0.01071544575250853
		//-0.012780945066091901
		;// +0.0075236204997510547;// -0.00014188594818104333;// -1.075125297987262107;
	pitch = -0.75948931401617531 + 0.0420430719632504 + 0.0148339286 - 0.0010836224 + 0.030592053542564165
		+ -0.001097921209459401 - 0.0010067393937711834
		+ 0.0063446391300211588
		- 0.0032938231305866673
		+ 0.0056215699394276666
		+ 0.007416771062703764
		-0.006327730102806307
		//- 0.0086598314788523534
		;//	+0.0025349011363896587;// +0.0058690842523870029;// -0.24828320355427008;
	
	yaw = -2.705626471380425914 + 0.042342165339081667 + 0.059684815497401282 + 0.003482870997925182 + 0.059289455776508071
		+ -0.0026948754918504044 - 0.0014221279183076598
		+ 0.012147362676649969
		- 0.0048444526234804163
		+ 0.0077773807892022769
		+ 0.016194812304292105
		-0.014090265101906628
		//- 0.010286895668387488
		;//	+0.0051545229779616231;// -0.0024537623203155908;// -1.576824892399777614 + 0.01;
	double rotvecvalue[3] = { roll / 180.0 * CV_PI, pitch / 180.0 * CV_PI, yaw / 180.0 * CV_PI };
	Mat rotatevector(3, 1, CV_64F, &rotvecvalue);

	double rotvalue[9];
	CvMat value1 = rotatevector;
	Mat rotationMatrix(3, 3, CV_64F, &rotvalue);
	CvMat value2 = rotationMatrix;
	cvRodrigues2(&value1, &value2);

	

	orbithelper.mult(lightvector, rotvalue, _lightvector, 1, 3, 3);
	
	sumtest = 0;
	for (int i = 0; i < 3; i++)
		sumtest += _lightvector[i] * _lightvector[i];

	sumtest = sqrt(sumtest);

	if (sumtest != 1)
	{
		//cout << "sumtest is :" << sumtest << endl;
		for (int i = 0; i < 3; i++)
			_lightvector[i] = _lightvector[i] / sumtest;
	}

	//_lightvector[0] = 0; _lightvector[1] = 0; _lightvector[2] = -1;
	/*
	double AA=(PosInCT[0]*PosInCT[0]+PosInCT[1]*PosInCT[1])/A/A+PosInCT[2]*PosInCT[2]/B/B;
	double BB =2*(XYZs[0]*PosInCT[0]+XYZs[1]*PosInCT[1])/A/A+2*XYZs[2]*PosInCT[2]/B/B;
	double CC=(XYZs[0]*XYZs[0]+XYZs[1]*XYZs[1])/A/A+XYZs[2]*XYZs[2]/B/B-1;
	double scale=(-BB-sqrt(BB*BB-4*AA*CC))/2.0/AA;
	*/

	double AA = (_lightvector[0] * _lightvector[0] + _lightvector[1] * _lightvector[1]) / A / A + _lightvector[2] * _lightvector[2] / B / B;
	double BB = 2 * (_satposition[0] * _lightvector[0] + _satposition[1] * _lightvector[1]) / A / A + 2 * _satposition[2] * _lightvector[2] / B / B;
	double CC = (_satposition[0] * _satposition[0] + _satposition[1] * _satposition[1]) / A / A + _satposition[2] * _satposition[2] / B / B - 1;

	double middlevalue = 0;
	middlevalue = (BB*BB - 4 * AA * CC);
	if (((BB*BB - 4 * AA * CC)) < 0)
	{
		cout << "姿态数据错误或者转交发生变化，导致卫星相机光束与地面没有交集，请检查日志" << endl;
		IAddLog("姿态光束与地球椭球不相交，请检查转角模式与星上数据是否相符");
		return false;
	}

	double scale = (-BB - sqrt(BB*BB - 4 * AA*CC)) / 2.0 / AA;




	DATUM WGS84[] =
	{ _a, _b, _c, _d, _e, _f, _g };
	double X, Y, Z;
	rectPosition[0] = _satposition[0] + scale*_lightvector[0];
	rectPosition[1] = _satposition[1] + scale*_lightvector[1];
	rectPosition[2] = _satposition[2] + scale*_lightvector[2];

	double test1, test2, test3;
	orbithelper.rect2geograph(_satposition[0], _satposition[1], _satposition[2], WGS84, &test1, &test2, &test3);
	test1 = test1 / 3.1415926535897932384626433832795 * 180.0;
	test2 = test2 / 3.1415926535897932384626433832795 * 180.0;
	orbithelper.rect2geograph(rectPosition[0], rectPosition[1], rectPosition[2], WGS84, &earthposition[0], &earthposition[1], &earthposition[2]);
	earthposition[1] = -earthposition[1];
	rotatevector.~Mat();
	rotationMatrix.~Mat();
	return true;
}
