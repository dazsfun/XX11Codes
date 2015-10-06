#ifndef _VIRTUALMATHICN_HEADER_H_INFO_
#define _VIRTUALMATHICN_HEADER_H_INFO_
#pragma once
#include "log.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "Scanner.h"
#include "jsonSerializer.h"
#include "LineTimeAPI.h"
#include "OrbitAPI.h"
#include "AttitudeAPI.h"
#include "UTCAPI.h"

#include "RPCEvaluationHelper.h"
#include "ImageModel.h"
#include "ImgPoint.h"
#include "EarthModel.h"
#include "CalibratedModel.h"
#include "RadioMetricProcess.h"
#include "RPCModel.h"

#pragma comment(lib,"RCProcess.lib")
using namespace std;

enum CEnumMachineCompent
{
	MachineComponent_asOrbit,
	MachineComponent_asAttitude,
	MachineComponent_asLienTime,
};

///本类的作用是定义了Fred文件的规格后，实施FRED内容的读取或者XML的读取或者JSON文件的读取
///并且建立最基础的成像模型，在这个模型中，不考虑内方位元素和偏置，摆扫角度等
///它存在的好处：1 。继承的利用：是让其子类（更高阶成像类）设计时不需要再考虑辅助数据的组建
///             2 。多态的利用：不用再去关心文件读取代码的实现，将格式保存为极简单的配置文件
///                 交给同一套代码去处理
/*与成像模型有关的两个主要接口
	virtual IGetPosition(double _time);
	virtual IDirction(double _x = 0,double _y = 0);
*/

class RSRPCModel;
class VirtualMachine : public LogAPI
{
public:
	UTCAPI GetTime(CEnumMachineCompent _utcapi, int _index);
private:
	bool GetAuxFile(string auxxml,CEnumMachineCompent _componenttype);
	bool GetJSonFile(string _filepath,CEnumMachineCompent _componenttype);
	map<CEnumMachineCompent, string> _compent_lookuptable;
	void SetLookUpTable();
	bool BuildRpcModel(int index);
	vector<RSRPCModel> _rpcInfo;

	double comphi, comome, comkappa;

	void Form_RPC(int index_id);

	double GetCompensate()
	{
		ifstream ifs("C:\\Temp\\match\\compensate.txt");
		char reader[1024];
		ifs.getline(reader, 1024);
		sscanf_s(reader, "%lf %lf %lf", &comphi, &comome, &comkappa);
		ifs.close();
		return 0;
	}



public:
	double GiveTime(double time)
	{
		string hei = IGetPosition(time);
		OrbitAPI apiforOrbit(hei, asJSon);
		return apiforOrbit.utcapi.GetUTC();
	}

	 

	void FromXY2LonLatTest(const double& imgx, const double& imgy, const double& height, double& lon, double& lat,double& rot1,double& rot2,double& rot3,double* posinfo, int imageid)
	{
		imageid -= 41;
		//currentImageID = imageid;
		if (imageid == -1) imageid = currentImageID;

		///根据行号和景ID确定时间 
		int _lineindex = imageid * _singleframeHeight + _singleframeHeight - 1 - imgy;
		int _oriline = imageid * _singleframeHeight;
		double time = _simpleTimes[_lineindex].utctime;  //时间的设置 
		double addon1 = 10 - fabs(imgy - 5392.5) / 5392.5 * 10.0;
		time += addon1 / 480.0 * 10786 * 3.5122304903644342682694142339248e-5;
		//time -= 10 * 3.5122304903644342682694142339248e-5;
		cout.precision(20);
		//cout << "----current utc time is: " << time<<"----";

		//cout << "whole time is" << _simpleTimes[(imageid + 1) * _singleframeHeight].utctime - _simpleTimes[imageid * _singleframeHeight].utctime << endl;
		//cout << "starttimesi:" << timestartapi.year << "-" << timestartapi.month << "-" << timestartapi.day << "-" << timestartapi.hour << "-" << timestartapi.minute << "-" << timestartapi.second << endl;

		//cout << "Average time is" << _simpleTimes[(imageid + 1) * _singleframeHeight].avgtime<<endl;
		double body2WGS84Rotation[9] = { 0 };
		UTCAPI timecheck(time, Enum_TimeStandard_ZY03);
		double position[3] = { 0 };
		m_pImgModel->IBody2WGS84Test(time, body2WGS84Rotation, position);
		//cout << "body 2 wgs84 finished" << endl;
		memcpy(posinfo, position, sizeof(double)* 3);

		double lightvector[3] = { 0 };
		base_point.ILightVector(imgx, imgy, lightvector);
		double earth_LonLatposition[3] = { 0 };
		double earth_RectPosition[3] = { 0 };
		double lookvector[3] = { 0 };

		CSatOrbit helper;
		//helper.invers_matrix(body2WGS84Rotation, 3);
		double rotmirror[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
		double rotcamera[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
		//{ 0.0141, 89.9914, 90.0111,
		//90.0887, 0.0292, 89.9721,
		//89.9888, 90.0227, 0.0253 };

		double rotaftercamera[9]; m_base.Multi(rotcamera, body2WGS84Rotation, rotaftercamera, 3, 3, 3);
		double rotaftermirror[9]; m_base.Multi(rotmirror, rotaftercamera, rotaftermirror, 3, 3, 3);


		//memcpy(body2WGS84Rotation, rotaftermirror, sizeof(double)* 9);




		double roll = 0;
		double pitch = 0;
		double yaw = 0;// -1.4 + 2.38594403;
		double rotvecvalue[3] = { roll / 180.0 * CV_PI, pitch / 180.0 * CV_PI, yaw / 180.0 * CV_PI };
		Mat rotatevector(3, 1, CV_64F, &rotvecvalue);



		double rotvalue[9];
		CvMat value1 = rotatevector;
		Mat rotationMatrix(3, 3, CV_64F, &rotvalue);
		CvMat value2 = rotationMatrix;
		cvRodrigues2(&value1, &value2);

		double swingrot[9];
		double startpos = (-2.0875 * CV_PI / 180.0);
		double swingoffset = 0;// (0.5* 0.0381971863420548805845321032094 / (9226.0 - / 9226.0 / 9226.0) * CV_PI / 180.0 * imgy*imgy*imgy;
		double wholevalue = +fabs(startpos) * 2 * 1.00549388739566009032939488764 * 1.000250944846080260950021170825 / 1.0057495352300516595331044322155 * 1.000931158 / 1.0001140023647047061221734279613;
		;
		double imgyshift = 0;


		double imgx_yshift = 0;
		if (p_TimeCalib_X != nullptr)
		{
			imgx_yshift = GetAlongAngle(imgx, true);
		}


		double valueget = 0;
		if (_gofuckyourselfSwingError.size() > 0)
		{
			//cout << "hei" << endl;
			if (imgy > _gofuckyourselfSwingError[0]._start && imgy < endswingvalue)
			{


				if (imgy < _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1]._start)
				{
					double heibegin = _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1]._start;
					int step = (heibegin - imgy) / 5.0;


					for (int i = 0; i < step; i++)
					{
						valueget += 5.0 * _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1 - 1 - i]._avg;
					}

					valueget += _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 2 - step]._avg * fabs(heibegin - 5.0 * step - imgy);

					valueget += _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1]._avg * (endswingvalue - heibegin);
				}

				else
				{
					valueget += _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1]._avg * fabs(endswingvalue - imgy);
				}
			}
		}

		double swingangle1 = 0;
		if (true/*fabs(valueget - 0) < 0.0001*/)
		{
			double addon = 6 - fabs(imgy - 5392.5) / 5392.5 * 6;
			swingangle1 = (startpos + (10785 - imgy) / 10785.0 * wholevalue) * 2 + (-5.0 - addon) / 10785.0 * wholevalue * 2;
		}

		else
		{
			//cout << "hei" << endl;
			swingangle1 = startpos * 2 + (10785 - endswingvalue) / 10785 * wholevalue * 2 + valueget / 10785 * wholevalue * 2;
		}
		//swingangle1 -= swingoffset;
		//double realswan = -swingangle1 * 2 + 0;
		double swingan[3] = { swingangle1, 0, 0 };
		//cout << "swing angle is " << swingangle1 << "iMAGEY IS" << imgy * 180.0 / CV_PI << endl;
		//m_base.rot(0, swingan,0, swingrot);
		Mat swingvector(3, 1, CV_64F, &swingan);
		CvMat values1 = swingvector;
		Mat rotationMatrixS(3, 3, CV_64F, &swingrot);
		CvMat values2 = rotationMatrixS;
		cvRodrigues2(&values1, &values2);
		double rotafterswing[9];
		//helper.rot(0, swingangle1, 0, swingrot);
		//m_base.invers_matrix(swingrot, 3);

		double rotafterCompensate[9];
		helper.mult(rotaftermirror, rotvalue, rotafterCompensate, 3, 3, 3);
		helper.mult(rotafterCompensate, swingrot, body2WGS84Rotation, 3, 3, 3);

		


		double rotcom1value[9];
		//cout << "comphi" << comphi << "\t" << comome << "\t" << comkappa << endl;
		safeopencv::Eulor2Rot(comphi, comome, comkappa, rotcom1value);		
		helper.invers_matrix(rotcom1value, 3);
		double wgs84camera[9];
		double wgs84cameravector[3];
		helper.mult(rotcom1value, body2WGS84Rotation, wgs84camera, 3, 3, 3);
		safeopencv::Rot2Eulor(wgs84camera, rot1, rot2, rot3);
	
		//double lightvectorbeforehand[3] = { 0 };
		helper.mult(body2WGS84Rotation, lightvector, lookvector, 3, 3, 1);
		//helper.mult(lightvectorbeforehand, swingrot, lookvector, 1, 3, 3);
		pEarth->IPosition(position, lookvector, height, earth_RectPosition, earth_LonLatposition);
		earth_LonLatposition[0] = earth_LonLatposition[0] / 3.1415926535897932384626433832795 * 180.0;
		earth_LonLatposition[1] = earth_LonLatposition[1] / 3.1415926535897932384626433832795 * 180.0;


		lon = earth_LonLatposition[1];
		lat = earth_LonLatposition[0];
	}

	double GetRangeForOneStep(int ccdindex);
private:
	vector<RPCEvaluationHelper> rpc_container;

	void GetControlPoints(string matchfile, vector<ImgPoint>& allinfo,bool bwithoutground = false)
	{
///38.061753410420017	90.181367931438203	43.33130719125591	84.431024003666764		-0.93824658957998175	-0.81863206856179105	 0.33130719125591046	2.4310240036667667	0.9647651438266448
		ifstream heimatch(matchfile.c_str());
		char reader[1024];
		ImgPoint newpoints(base_point);
		ofstream constream("c:\\Temp\\match\\compensatelist.txt");
		ifstream groundpoint("c:\\Temp\\match\\groundlist.txt");
		ofstream heiground;
		bool groundexist = groundpoint.is_open();
		if (!groundexist)
		{
			groundpoint.close();
			heiground.precision(20);
			heiground.open("c:\\Temp\\match\\groundlist.txt");
		}
		constream.precision(20);
		int ncounter = 0; 
		do
		{
			cout << ncounter << "\r";
			ncounter++;
			heimatch.getline(reader, 1024);
			double a[9];
			sscanf_s(reader, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &a[0], &a[1],
				&a[2], &a[3], &a[4], &a[5], &a[6], &a[7], &a[8]);

			//cout << a[0] << "\t" << a[1] << "\t" << a[2] << "\t" << a[3] << endl;
			///确定时间
			int _currentImageID = 21;
			newpoints.ImgX(a[0]);
			newpoints.ImgY(a[1]);
			newpoints._rotation4 = a[2];
			newpoints._rotation5 = a[3];
			int _lineindex = _currentImageID * _singleframeHeight + _singleframeHeight - newpoints.ImgY();
			double time = _simpleTimes[_lineindex].utctime;  //时间的设置 
			double b[5];
			char readground[1024];
			if (groundexist)
			{
				groundpoint.getline(readground, 1024);
				sscanf_s(readground, "%lf %lf %lf %lf %lf", &b[0], &b[1], &b[2], &b[3], &b[4]);
			}

			newpoints.utcapi = UTCAPI(time, Enum_TimeStandard_ZY03);
			//newpoints = newpoints + base_point;
			double xvalue, yvalue, zvalue;
			double heilon, heilat;
			double heilonnext, heilatnext;
			////Amend Here!!!  20150512 1800
			////在这一过程中需要知道整个大矩阵的值，使用一个新的函数在FromXY2LonLatTest基础上完成任务
			///Add Function:——FromXY2Rot(newpoints.ImgX(), newpoints.ImgY(), 100,rot1,rot2,rot3,62);
			///Its interface should be like above 
			///预计完成时间：  20分钟
			///完成
			double posinfo[3] = { 0 };
			FromXY2LonLatTest(newpoints.ImgX(), newpoints.ImgY(), 100, heilon, heilat, newpoints._rotation7, newpoints._rotation8, newpoints._rotation9, posinfo,
				 66);
			//cout << newpoints._rotation7 << "\t" << newpoints._rotation8 << "\t" << newpoints._rotation9 << endl;
			pEarth->Geo2Rect(-heilon* CV_PI / 180.0, heilat * CV_PI / 180.0, 100, xvalue, yvalue, zvalue);
			
			//真实经纬度
			double xposvalue, yposvalue, zposvalue;
			///仅计算一次
			if (!groundexist)
			{
				FromXY2LonLatTest(a[2], a[3], 100, heilonnext, heilatnext, 66);
				newpoints._longitude = heilonnext;
				newpoints._latitude = heilatnext;

				pEarth->Geo2Rect(-heilonnext * CV_PI / 180.0, heilatnext * CV_PI / 180.0, 100, xposvalue, yposvalue, zposvalue);
				heiground << newpoints._longitude << "\t" << newpoints._latitude << "\t" << xposvalue << "\t" << yposvalue << "\t" << zposvalue << endl;
			}
			else{
				newpoints._longitude = b[0];
				newpoints._latitude = b[1];
				xposvalue = b[2];
				yposvalue = b[3];
				zposvalue = b[4];
			}

			if (!bwithoutground)
			{
				newpoints._rotation4 = posinfo[0]; newpoints._rotation5 = posinfo[1];
				newpoints._rotation6 = posinfo[2];
			}

			double test1, test2, test3;
			newpoints.OriX(xvalue);
			newpoints.OriY(yvalue);
			newpoints.OriZ(zvalue);
			newpoints.X(xposvalue); newpoints.Y(yposvalue); newpoints.Z(zposvalue);
			//newpoints.Serialize("d:\\testPos\\0410\\");
			double heisqrt = 0;
			heisqrt += pow(newpoints.X() - newpoints.OriX(), 2.0);
			heisqrt += pow(newpoints.Y() - newpoints.OriY(), 2.0);
			heisqrt += pow(newpoints.Z() - newpoints.OriZ(), 2.0);
			heisqrt = sqrt(heisqrt);
			constream << heisqrt << "\t" << heilon << "\t" << heilat << "\t" << newpoints.Lon() << "\t" << newpoints.Lat() << endl;
			if (newpoints.Lat()< 54 || newpoints.Lat() > 58 || newpoints.Lon() < -5 || newpoints.Lon() > -2)
			{
				cout << newpoints.Lat() << "\t" << newpoints.Lon() << endl;
				cout << "控制点坐标计算错误！！" << endl;
				continue;
			}
			allinfo.push_back(ImgPoint(newpoints,false));
		} while (heimatch.peek() != EOF);

		constream.close();
	}

	vector<SwingCablib> _gofuckyourselfSwingError;
	void ReadSwingError()
	{
		_gofuckyourselfSwingError.clear();
		ifstream readsw("c:\\Temp\\match\\restresult.xls");
		if (!readsw.is_open())
		{
			return;
		}

		else
		{
			do
			{
				char reader[1024];
				readsw.getline(reader, 1024);
				double a[4];
				sscanf_s(reader, "%lf %lf %lf %lf", &a[0], &a[1], &a[2], &a[3]);
				SwingCablib calib;
				calib._start = a[0]; calib._end = a[1];
				calib._avg = a[2]; calib._rms = a[3];

				_gofuckyourselfSwingError.push_back(calib);
				endswingvalue = 10494.43;// _gofuckyourselfSwingError[_gofuckyourselfSwingError.size() - 1]._end;
			} while (readsw.peek() != EOF);
		}
	}


public:
	virtual string IGetPosition(double _time);
	int SingHeight()
	{
		return _singleframeHeight;
	}

	int SingWidth()
	{
		return _singleframeWidth;
	}
private:
	double endswingvalue = 0;
	virtual AttitudeAPI GetAttitude(double _utctime);
	vector<TimeCalibStruct>* p_TimeCalib;
	vector<TimeCalibStruct>* p_TimeCalib_X;
	double GetAlongAngle(double yvalue,bool bxarray = false);

	//ImgPoint::ImgPoint(vector<double> swingangleTable, vector<double> lookangleXtable, vector<double> lookangleYtable, double* utctime, double* rotationmatrix)
	ImgPoint SetPoint()
	{
		vector<double> _swing;
		_swing.resize(10786);
		for (int i = 0; i < 10786; i++)
		{
			_swing[i] = -2.0875 + (2.0875 *2 )/ 10786 * i;
		}
		vector<double> _lookangle;
		_lookangle.resize(480);
		for (int i = 0; i < 480; i++)
		{
			_lookangle[i] = 0;
		}

		double time = 0;
		double* rotaMatrix = new double[9];
		memset(rotaMatrix, 0, sizeof(double)* 9);
		ImgPoint imgp(_swing, _lookangle, _lookangle, &time, rotaMatrix);
		delete[] rotaMatrix;
		return imgp;
	}

private:
	bool GetNeightbor(int indexid,vector<Point2d>& _leftmatch,vector<Point2d>& _rightmatch);
	bool PreciseBasedAffine(int indexid, double &x, double &y, double latitude, double longitude, double H, double dx, double dy);
		bool PredictbasedAffine(int indexid, double latitude, double longitude, double H, double &x, double &y);
	CSatOrbit orbithelper;

public:
	double distanceOnePixel(int indexid, const double& x, const double& y, const double& H, const double& dx, const double& dy);
	bool FromlatlonH2xy(double lat, double lon, double H, double &x, double &y, int ccdID);
	void FromXY2LONLATRPC(const double& imgx, const double& imgy, const double& height, double& lon, double&lat,const int& index);
	void FromXY2LonLatTest(const double& imgx, const double& imgy, const double& height, double& lon, double& lat, int imageid);
	void Geo2Rect(const double& lon, const double& lat, const double& height, double& xinfo, double& yinfo, double& zinfo)
	{
		pEarth->Geo2Rect(lon, lat, height, xinfo, yinfo, zinfo);
	}
private:
	void FillPointWithModelValue(vector<ImgPoint> pointarray);

	ImgPoint SetPoint(int swingrange, int ccdrange)
	{
		ccdrange = 480;
		vector<double> _swing;
		_swing.resize(swingrange);
		for (int i = 0; i < swingrange; i++)
		{
			_swing[i] = -2.0875 + (2.0875 * 2) / swingrange * i;
		}
		vector<double> _lookanglex;
		vector<double> _lookangley;
		_lookanglex.resize(ccdrange + 400);
		_lookangley.resize(ccdrange+ 400);
		double focuslength = 2.4571908686868686868686868686869;// *38 / 50.8;
		
		for (int i = 0; i < 880; i++)
		{
			int number_i = i;
			i -= 200;
			double yshift = 0;
			if (true/*p_TimeCalib_X != nullptr*/)
			{
				//yshift = GetAlongAngle(i, true);
				yshift = 19.214285714285714285714285714286 - static_cast<double>(i)* (24.0 / 420.0) - (0 - static_cast<double>(i)* (6.6 / 480.0));
				yshift /= 1;
				//cout << "!!";
			}



			if (i == 240)
			{
				cout << "探元:" << i << "畸变值:" << yshift << endl;
			}
			_lookangley[number_i] = ((232.4 / 2 * 10e-7 - yshift  * 38.0 * 10e-7 /*+ 49.8 * 5*/) / focuslength);
			_lookanglex[number_i] = ((-25.4 - 50.8 * (239) + i * 50.8) * 10e-7 / focuslength);
			i = number_i;
		}

		double time = 0;
		double* rotaMatrix = new double[9];
		memset(rotaMatrix, 0, sizeof(double)* 9);
		rotaMatrix[0] = rotaMatrix[4] = rotaMatrix[8] = 1.0;
		ImgPoint imgp(_swing, _lookanglex, _lookangley, &time, rotaMatrix);
		delete[] rotaMatrix;
		return imgp;
	}

	double perT;
	double Vel(int imageid)
	{
		imageid -= 41;
		///根据行号和景ID确定时间 
		int _lineindex = imageid * _singleframeHeight;
		int _oriline = imageid * _singleframeHeight + _singleframeHeight;
		double time = _simpleTimes[_lineindex].utctime;  //时间的设置 
		return m_pImgModel->VelX(time);
	}

	double CalcCrabAngle(int imageid)
	{
		imageid -= 41;
		//currentImageID = imageid;
		if (imageid == -1) imageid = currentImageID;

		///根据行号和景ID确定时间 
		int _lineindex = imageid * _singleframeHeight;
		int _oriline = imageid * _singleframeHeight + _singleframeHeight;
		double time = _simpleTimes[_lineindex].utctime;  //时间的设置 
		double timenext = time + perT;

		return (m_pImgModel->VelAngle(timenext) - m_pImgModel->VelAngle(time));
	}

private:
	virtual void SetPropertys();
	GeoBase m_base;
	ImageModel FormImageModel();
	vector<SimpleTime> _simpleTimes;
	int currentindex;
public:
	ImageModel* m_pImgModel;

public:
	int Current()
	{
		return currentindex;
	}
	VirtualMachine(string auxxml, bool withimage = false);

	VirtualMachine(string auxxml) :base_point(SetPoint())
	{
		/*
		SetPropertys();
		///开始规格的读取或者配置
		cout << "请给出配置文件信息：" << endl;
		char _configinfo[1024];
		cin >> _configinfo;
		cout << "FRED文件路径" << endl;
		char _fredpathp[1024];
		cin >> _fredpathp;
		EnvelopeAPI api1(_configinfo, _fredpathp);
		if (!api1.Is_Inited())
		{
			EnvelopeAPI api1(_fredpathp, true);
			///现在我们已经获取了一帧信息，
			api1.IFindStart();
			///接下来的工作，是指定几帧构成一景影像信息，而每一帧中辅助信息如何组织
			cout << "请输入一景的长度（范围1~20）:" << endl;
			int num = 0; 
			cin >> num;
			api1.FrameForOneScene(num);
			///现在的问题是，如何让程序识别辅助数据的组织形式
		}
		*/
	}

	void Calibration();

	VirtualMachine(string dirpath, CEnumDomType _jsonobject);


public:
	VirtualMachine(string fredconfig, string fredpath);

public:
	int GetSize(CEnumMachineCompent);
	
private:
	vector<LineTimeAPI> _timeapi;
	vector<OrbitAPI> _orbitApi;
	vector<AttitudeAPI> _attitudeapi;
	vector<AuxFactory> _myFactory;
	ImgPoint base_point;
	int currentImageID;
	EarthModel* pEarth;
	void SetEarthModel();
	string OutPutTime();
	void SetSimpleTime(vector<int> frameheight);
	int _singleframeHeight;
	int _singleframeWidth;
public:
	void FromXY2LonLat(const double& imgx, const double& imgy, const double& height, double& lon, double& lat, int imageid = -1);
	void FromXY2LonLatUnitTest(const double& imgx, const double& imgy, const double& height, double& lon, double& lat, int imageid, string inputer);
public:
	virtual ~VirtualMachine() {}
	
protected:
	void DetailedLineTimeInfo(vector<SimpleTime> _info,int frameAgg);
};
#endif