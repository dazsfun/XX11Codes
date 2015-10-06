#pragma once
#include <vector>
#include "jsonSerializer.h"
#include "SatOrbitNew.h"
#include <math.h>
using namespace std;
///这个类关心的有4个因素：
///1. 成像时间的导入
///2. 成像时刻指向角（必须是严格意义上的指向角）的导入，不包含外畸变的改正
///3. 成像时刻偏置（外畸变改正的导入），static变量，对一景影像的处理来说，它肯定是一个定值
///4. 
///在被初始化之后，对于imgPoint来说就没有成像模型的概念了
class ImgPoint :
	public CJson_Ojbect_Base_Serializer
{
public:
	ImgPoint(vector<double> lookangletable, vector<double> swingangletable,vector<double> swingangle1table, double* utctime,double* rotationmatrix);
	~ImgPoint()
	{
		if (_rotationmatrix != nullptr)
		{
			//delete[] _rotationmatrix;
			_rotationmatrix = nullptr;
		}

		if (_utctime != nullptr)
		{
			//delete[] _utctime;
			_utctime = nullptr;
		}
	}

	void SetAngle(int index, double _xinfo, double _yinfo)
	{
		ofstream ofs_angle;
		ofs_angle.open("c:\\Temp\\match\\look_error_test.txt", ios::app);
		ofs_angle << index << "\t" << (_xinfo - _lookangleXtable[index]) / _lookangleXtable[index]
			<< "\t" << (_yinfo - _lookangleYtable[index]) / _lookangleYtable[index] << endl;
		ofs_angle.close();
		_lookangleXtable[index+ 35] = _xinfo;
		_lookangleYtable[index+ 35] = _yinfo;
	}

	ImgPoint(string info, CEnumDomType domtype)
	{
		objectname = "ImgPoint";
		SetType(domtype);
		
		_rotationmatrix = new double[9];
		double tempvaluerot[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

		memcpy(_rotationmatrix, tempvaluerot, sizeof(double)* 9);
		DeSerialize(info.c_str());
	}

	void ActualVector(double* satposition, double*  calcpostion, double* actlook)
	{
		double linevector[3] = { satposition[0] - calcpostion[0], satposition[1] - calcpostion[1], satposition[2] - calcpostion[2] };
		double sumvalue = 0;
		for (int i = 0; i < 3; i++)
		{
			sumvalue += pow(linevector[i], 2.0);
		}
		for (int i = 0; i < 3; i++)
		{
			actlook[i] = -linevector[i] / sqrt(sumvalue);
		}
	}

	void TrueVector(double* satposition, double* lookvector)
	{
		double linevector[3] = { satposition[0] - _X, satposition[1] - _Y, satposition[2] - _Z };
		double sumvalue = 0;
		for (int i = 0; i < 3; i++)
		{
			sumvalue += pow(linevector[i], 2.0);
		}
		for (int i = 0; i < 3; i++)
		{
			lookvector[i] = -linevector[i]/ sqrt(sumvalue);
		}

	}

	ImgPoint operator+(const ImgPoint& _point)
	{
		safestring::CopyVector(this->_swingangletable, _point._swingangletable);
		safestring::CopyVector(this->_lookangleXtable, _point._lookangleXtable);
		safestring::CopyVector(this->_lookangleYtable, _point._lookangleYtable);
		
		return *this;
	}

	void operator/(const ImgPoint& _point)
	{
		_utctime = _point._utctime;
		//memcpy(_rotationB2W, _point._rotationB2W, sizeof(double)* 9);
		_imgx = _point._imgx;
		_imgy = _point._imgy;

		_swingangle = _point._swingangle;
		_lookangleX = _point._lookangleX;
		_lookangleY = _point._lookangleY;
		_longitude = _point._longitude;
		_latitude = _point._latitude;
		_centerlon = _point._centerlon;
		_X = _point._X;
		_Y = _point._Y;
		_Z = _point._Z;
		_rotation1 = _point._rotation1;
		_rotation2 = _point._rotation2;
		_rotation3 = _point._rotation3;
		_rotation4 = _point._rotation4;
		_rotation5 = _point._rotation5;
		_rotation6 = _point._rotation6;
		_rotation7 = _point._rotation7;
		_rotation8 = _point._rotation8;
		_rotation9 = _point._rotation9;

		utcapi = _point.utcapi;
	}

	ImgPoint(const ImgPoint& _pointnow, bool iftrue)
	{
		SetPropertys();
		if (!iftrue)
		{
			(*this) / _pointnow;
		}
		else
			*this = _pointnow;
	}

	ImgPoint(const ImgPoint& _point)
	{
		SetPropertys();
		_X = _point._X;
		_swingangle = _point._swingangle;
		//_swingangletable = _point._swingangletable;
		//_lookangleXtable = _point._lookangleXtable;
		//_lookangleYtable = _point._lookangleYtable;
		if (_point._swingangletable.size() >0)
		safestring::CopyVector(_swingangletable, _point._swingangletable);
		if (_point._lookangleXtable.size() >0)
		safestring::CopyVector(_lookangleXtable, _point._lookangleXtable);
		if (_point._lookangleYtable.size() > 0)
		safestring::CopyVector(_lookangleYtable, _point._lookangleYtable);

		//memcpy(_rotationmatrix, _point._rotationmatrix, sizeof(double)* 3);
		_utctime = _point._utctime;
		//memcpy(_rotationB2W, _point._rotationB2W, sizeof(double)* 9);
		_imgx = _point._imgx;
		_imgy = _point._imgy;

		_swingangle = _point._swingangle;
		_lookangleX = _point._lookangleX;
		_lookangleY = _point._lookangleY;
		_longitude = _point._longitude;
		_latitude = _point._latitude;
		_centerlon = _point._centerlon;
		_X = _point._X;
		_Y = _point._Y;
		_Z = _point._Z;
		_rotation1 = _point._rotation1;
		_rotation2 = _point._rotation2;
		_rotation3 = _point._rotation3;
		_rotation4 = _point._rotation4;
		_rotation5 = _point._rotation5;
		_rotation6 = _point._rotation6;
		_rotation7 = _point._rotation7;
		_rotation8 = _point._rotation8;
		_rotation9 = _point._rotation9;

		utcapi = _point.utcapi;
	}
	virtual void SetUTCAPIByTimestring() {};

	const double& OriX()
	{
		return _rotation1;
	}

	const double& OriY()
	{
		return _rotation2;
	}

	const double& OriZ()
	{
		return _rotation3;
	}

	void OriX(const double& valueinput)
	{
		_rotation1 = valueinput;
	}
	
	void OriY(const double& valueinput)
	{
		_rotation2 = valueinput;
	}

	void OriZ(const double& valueinput)
	{
		_rotation3 = valueinput;
	}

	const double& ImgX()
	{
		return _imgx;
	}

	void ImgX(const double& value)
	{
		_imgx = value;
	}

	void ImgY(const double& value)
	{
		_imgy = value;
	}

	const double& ImgY()
	{
		return _imgy;
	}

	const double& Lon()
	{
		return _longitude;
	}

	const double& Lat()
	{
		return _latitude;
	}

	const double& X()
	{
		return _X;
	}

	const double& Y()
	{
		return _Y;
	}

	const double& Z()
	{
		return _Z;
	}

	void X(const double& xvalue)
	{
		_X = xvalue;
	}

	void Y(const double& yvalue)
	{
		_Y = yvalue;
	}

	void Z(const double& zvalue)
	{
		_Z = zvalue;
	}

public:
	void ILightVector(double imgx = 0, double imgy = 0, double* lightvector = nullptr);
	void ISetB2W(double* _body2wgs84);
	void IGetB2W(double* position,double* _body2wgs84,DATUM* WGS84);
private:
	  vector<double> _swingangletable; //不进行序列化
	  vector<double> _lookangleXtable;//同上
	 vector<double> _lookangleYtable;
	double* _utctime; //同上
	 double* _rotationmatrix; // 同上
	double _rotationB2W[9];
	double _imgx;   ///参与序列化
	double _imgy;   ///同上

public:
	double _swingangle; //同上
	double _lookangleX;//同上
	double _lookangleY; //同上
	double _longitude;  //同上
	double _latitude; //同上
	double _centerlon; //同上
	double _X; //同上    
	double _Y; //同上 
	double _Z;//同上 
	double _rotation1;
	double _rotation2;
	double _rotation3;
	double _rotation4;
	double _rotation5;
	double _rotation6;
	double _rotation7;
	double _rotation8;
	double _rotation9;

private:
	virtual void SetPropertys();
};

