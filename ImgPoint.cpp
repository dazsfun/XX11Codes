#include "stdafx.h"
#include "ImgPoint.h"
#include "Safeauxilary.h"

ImgPoint::ImgPoint(vector<double> swingangleTable, vector<double> lookangleXtable, vector<double> lookangleYtable, double* utctime, double* rotationmatrix)
: _utctime(utctime)
{
	SetPropertys();
	safestring::CopyVector(_lookangleXtable, lookangleXtable);
	safestring::CopyVector(_lookangleYtable, lookangleYtable);
	safestring::CopyVector(_swingangletable, swingangleTable);
	//_rotationmatrix = rotationmatrix;
	_rotationmatrix = new double[9];

	memcpy(_rotationmatrix, rotationmatrix, sizeof(double)* 9);
}

void ImgPoint::ISetB2W(double* body2WGS84)
{
	memcpy(_rotationB2W, body2WGS84, sizeof(double)* 9);

	_rotation1 = _rotationB2W[0];
	_rotation2 = _rotationB2W[1];
	_rotation3 = _rotationB2W[2];
	_rotation4 = _rotationB2W[3];
	_rotation5 = _rotationB2W[4];
	_rotation6 = _rotationB2W[5];
	_rotation7 = _rotationB2W[6];
	_rotation8 = _rotationB2W[7];
	_rotation9 = _rotationB2W[8];
}


void ImgPoint::ILightVector(double imgx, double imgy, double* lightvector)
{



	double lookanglex = Safeauxilary::InterPolation(imgx + 200.0, _lookangleXtable);
	double lookangley = Safeauxilary::InterPolation(imgx + 200.0, _lookangleYtable);

	lightvector[0] = (lookanglex);
	lightvector[1] =  (lookangley);
	lightvector[2] = 1.0;

	_lookangleX = lookanglex;
	_lookangleY = lookangley;

	_imgx = imgx;
	_imgy = imgy;

	double modevalue = 0;
	for (int i = 0; i < 3; i++)
		modevalue += pow(lightvector[i], 2.0);

	for (int i = 0; i < 3; i++)
		lightvector[i] /= sqrt( modevalue);
}

void ImgPoint::IGetB2W(double* position, double* _body2wgs84,DATUM* WGS84)
{
	if (_X == -99 && _longitude == -99)
	{
		cout << "����Ӱ���Ϊ��ʼ����Ͼͱ�����ƽ����㣬���飡" << endl;
		IAddLog("Ӱ���δ��ʼ��������ƽ�����");
	}

	if (_X == -99)
	{
		CSatOrbit sat;
		double height = 100;
		double Xvalue, Yvalue, Zvalue;
		sat.geograph2rect(_latitude, _longitude, height,WGS84, &Xvalue, &Yvalue, &Zvalue);
	}
}

void ImgPoint::SetPropertys()
{
	/*
	double _imgx;   ///�������л�
	double _imgy;   ///ͬ��
	double _lookangle; //ͬ��
	double _swingangleX;//ͬ��
	double _swingangleY; //ͬ��
	double _longitude;  //ͬ��
	double _latitude; //ͬ��
	double _centerlon; //ͬ��
	double _X; //ͬ��
	double _Y; //ͬ��
	double _Z;//ͬ��
	double _rotation1;
	double _rotation2;
	double _rotation3;
	double _rotation4;
	double _rotation5;
	double _rotation6;
	double _rotation7;
	double _ratation8;
	double _ratation9;
	*/
	SetProperty("GroundX", asDouble, &_imgx);
	SetProperty("GroundY", asDouble, &_imgy);
	SetProperty("LookAngleX", asDouble, &_lookangleX);
	SetProperty("LookAngleY", asDouble, &_lookangleY);
	SetProperty("SwingAngle", asDouble, &_swingangle);
	SetProperty("longitude", asDouble, &_longitude);
	SetProperty("latitude", asDouble, &_latitude);
	SetProperty("centerlon", asDouble, &_centerlon);
	SetProperty("ImgX", asDouble, &_X);
	SetProperty("ImgY", asDouble, &_Y);
	SetProperty("GroundZ", asDouble, &_Z);
	SetProperty("RotationMatrix1", asDouble, &_rotation1);
	SetProperty("RotationMatrix2", asDouble, &_rotation2);
	SetProperty("RotationMatrix3", asDouble, &_rotation3);
	SetProperty("RotationMatrix4", asDouble, &_rotation4);
	SetProperty("RotationMatrix5", asDouble, &_rotation5);
	SetProperty("RotationMatrix6", asDouble, &_rotation6);
	SetProperty("RotationMatrix7", asDouble, &_rotation7);
	SetProperty("RotationMatrix8", asDouble, &_rotation8);
	SetProperty("RotationMatrix9", asDouble, &_rotation9);
}