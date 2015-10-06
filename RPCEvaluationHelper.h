#pragma once

#include <vector>
#include "BasicStructHeader.h"
#include "SatOrbitNew.h"
#include <fstream> 
#include <iostream>
#include <exception>
using namespace std;

enum RPCalculator
{
	Linear_POW1_DIFFER,
	Linear_POW2_DIFFER,
	Linear_POW2_DIFFER_Iteration,
	Linear_POW3_DIFFER,
	Linear_POW3_DIFFER_Iteration,
	Linear_POW1_Same,
	Linear_POW2_Same,
	Linear_POW2_Same_Iteration,
	Linear_POW3_Same,
	Linear_POW3_Same_Iteration,
	Linear_POW1_TrueSame,
	Linear_POW2_TrueSame,
	Linear_POW2_TrueSame_Iteration,
	Linear_POW3_TrueSame,
	Linear_POW3_TrueSame_Iteration,
	Linear_POW3_AllIteration,
};

class RPCEvaluationHelper
{
private:
	double ridgeconfigurationX;
	double ridgeconfigurationY;
	bool _modelinited;
public:
	bool Is_Inited()
	{
		return _modelinited;
	}

	void SetRidgeX(double ridgevalue)
	{
		ridgeconfigurationX = ridgevalue;
	}

	void SetRidgeY(double ridgevalue)
	{
		ridgeconfigurationY = ridgevalue;
	}
	RPCEvaluationHelper(void);
	~RPCEvaluationHelper(void);

public:
	void INewRPC(vector<FPoint> _newpoints, string pathnew);

	void FileName(string _name)
	{
		_checkfile = _name;
	}

	void ISelection(int* coeffmatrix)
	{
		RPCModel_Selected(coeffmatrix);
	}

public:
	ErrorType_1 IFromInnerAndOuterOrientation(string _descriptionJSon);

private:
	//static string _donominatortype[4];
	//static string _processtype[5];
	void RPCModel_Selected(int* coeffmatrix);
	void PreciseBasedAffine(double &latitude, double &longitude, double dlat, double dlon, double x0, double y0, double h0);
	double distanceOnePixel(double latitude, double longitude);
public:
	bool GetGoodBoy(const char* _rpcpath); /// get rpc file  and read it
	bool PushGoodBoy(const char* _rpcpath); ///output rpc paras into file with a standardformat;

	bool GetCommander(vector<FPoint>& _commanderinfo)
	{
		_checkfile = "c:\\Temp\\match\hei_likeit";
		PointList.clear();
		CPointList.clear();
		_CPointList.clear();
		for (int i = 0; i < _commanderinfo.size(); i++)
		{
			PointList.push_back(_commanderinfo[i]);
		}

		_commanderinfo.clear();
		return true;
	}
	bool GetCommander(const char*);   ///Get control point 
	bool GetCPoint(const char*);/// Get checkpoint
	void Compute_npara();//AAAADDDIINNGG
	void RPCModel_AllInter_Simplified(int* mark);
	double GetAvgH()
	{
		return np_array[4];
	}

protected:
	bool _P_F_AccCheck(double xcheck, double s, double ycheck, double t, double lon, double lat)
	{
		if (fabs(xcheck - s) > 0.01 || fabs(ycheck - t) > 0.01)
		{
			exit(EXIT_FAILURE);
			return false;
			char reader[1024];
			sprintf_s(reader, "RPC模型正反算异常,(%lf,%lf) 获取得到(%lf,%lf),经纬度(%lf,%lf)", s, t, xcheck, ycheck, lon, lat);
			cout << reader << endl;
			/*xception e(reader);
			throw e;*/

		}
		return true;
	}

private:
	double _imgwidth;
	double _imgheight;

public:
	void ISetWidth(double _width)
	{
		_imgwidth = _width;
	}

	void ISetHeight(double _height)
	{
		_imgheight = _height;
	}

	bool IRange(double _xvalue, double _yvalue, double _buffer)
	{
		double xstart = -_buffer;
		double xend = _imgwidth + _buffer;
		double ystart = -_buffer;
		double yend = _imgheight + _buffer;

		if (_xvalue < xstart || _xvalue > xend) return false;
		if (_yvalue < ystart || _yvalue > yend) return false;

		return true;
	}

private:
	void TranslateMark(int* mark, int* marks, int * markl, int * marka, int* markb);
	bool CheckSelection(int* mark, int& sizex, int & sizey, int size = 80)
	{
		sizex = sizey = 0;
		for (int i = 0; i < size / 2; i++)
		{
			if (mark[i] != 0)
			{
				sizey++;
			}
			else if (mark[i] == 0)
			{
				cout << "Y方向，第" << i + 1 << "列由于列相关被删除" << endl;
			}

			if (mark[i + size / 2] != 0)
			{
				sizex++;
			}
			else if (mark[i + size / 2] == 0)
			{
				cout << "X方向，第" << i + 1 << "列由于列相关被删除" << endl;
			}
		}

		if (sizey <= 28 || sizex <= 28)
			return false;
		else
			return true;
	}
	string _checkfile;
	bool m_bflag;
	RPCalculator _currentcalculator;
	

public:
	RPCEvaluationHelper(const RPCEvaluationHelper& _input)
	{
		memcpy(px, _input.px, sizeof(double)* 3);
		memcpy(py, _input.py, sizeof(double)* 3);
		memcpy(_s, _input._s, sizeof(double)* 20);
		memcpy(_l, _input._l, sizeof(double)* 20);
		memcpy(_a, _input._a, sizeof(double)* 20);
		memcpy(_b, _input._b, sizeof(double)* 20);

		memcpy(_ds, _input._ds, sizeof(double)* 20);
		memcpy(_dl, _input._dl, sizeof(double)* 20);
		memcpy(_da, _input._da, sizeof(double)* 20);
		memcpy(_db, _input._db, sizeof(double)* 20);

		memcpy(np_array, _input.np_array, sizeof(double)* 20);
		memcpy(np_array_Ind, _input.np_array_Ind, sizeof(double)* 20);

		_bInded = _input._bInded;
		bScaled = _input.bScaled;
		bSwitched = _input.bSwitched;

		 _checkfile = _input._checkfile;
		 m_bflag = _input.m_bflag;
		 _currentcalculator = _input._currentcalculator;

		  ridgeconfigurationX = _input.ridgeconfigurationX;
		 ridgeconfigurationY = _input.ridgeconfigurationY;
		  _modelinited = _input._modelinited;
	}

	RPCEvaluationHelper operator=(const RPCEvaluationHelper& _input)
	{
		memcpy(px, _input.px, sizeof(double)* 3);
		memcpy(py, _input.py, sizeof(double)* 3);
		memcpy(_s, _input._s, sizeof(double)* 20);
		memcpy(_l, _input._l, sizeof(double)* 20);
		memcpy(_a, _input._a, sizeof(double)* 20);
		memcpy(_b, _input._b, sizeof(double)* 20);

		memcpy(_ds, _input._ds, sizeof(double)* 20);
		memcpy(_dl, _input._dl, sizeof(double)* 20);
		memcpy(_da, _input._da, sizeof(double)* 20);
		memcpy(_db, _input._db, sizeof(double)* 20);

		memcpy(np_array, _input.np_array, sizeof(double)* 20);
		memcpy(np_array_Ind, _input.np_array_Ind, sizeof(double)* 20);

		_bInded = _input._bInded;
		bScaled = _input.bScaled;
		bSwitched = _input.bSwitched;

		_checkfile = _input._checkfile;
		m_bflag = _input.m_bflag;
		_currentcalculator = _input._currentcalculator;

		ridgeconfigurationX = _input.ridgeconfigurationX;
		ridgeconfigurationY = _input.ridgeconfigurationY;
		_modelinited = _input._modelinited;

		return *this;
	}

private:
	double px[3], py[3];
	double _s[20], _l[20], _a[20], _b[20];                      ///rpc paras
	double _ds[20], _dl[20], _da[20], _db[20];
	double np_array[10];                     //The normalized parameters
	double np_array_Ind[10];
	vector<FPoint> PointList;  //Controll points, keep them coming
	vector<SPoint> CPointList; //Bunch of points that tells your solution are awesome, or aweful

	//	CSatOrbit  _satHelper;  ///tool1
	CSatOrbit satOrbitHelper;

private:
	void GoodBoy1();       ///A method to solve rpc paras,Whatever it's ,  
	void GoodBoy2();     /// Whatever 2
	void GoodBoy3();    ///

private:
	void GetAccuracy(); ///Calculation of accuracy, easy and necessary

private:
	bool _bInded;
	bool bScaled;
	bool bSwitched;
	vector<FPoint> _CPointList;
private:

	void RPCPLH2a(double P, double L, double H, double *a);
	bool DoRPC1();
	bool DoRPC1Same();

	bool DoRPC1TrueSame();

	bool DoRPC2TrueSame(bool Iteration);

	bool DoRPC3TrueSame(bool Iteration);

	bool DoRPC2(bool bIteration = true);

	bool DoRPC3(bool bIteration = true);

	bool DoRPC2Same(bool bIteration = true);

	bool DoRPC3Same(bool bIteration = true);

	bool DoRPCAllInter();

	bool IterationRpc3();

	void Clear();



public:
	///These parts are what we call interface, 
	bool Switched()
	{
		return bSwitched;
	}

	bool ISWitch();

	bool IConductor(RPCalculator _messageer);

	bool IReport(bool _brpc)
	{
		//Checker();
		UpdateBasicInformation();
		if (_brpc)
			OutPutRpc();//if you dont want rpc,ok,delete this line 
		return true;
	}

	bool IFile(string _imgpath, bool _filein = true, bool _reversed = false);


	bool LightIn(double lon, double lat, double h, double& x, double& y);
	void LightOut(double x, double y, double h, double& lon, double& lat);

private:
	void Checker();
	void OutPutRpc();
	void UpdateBasicInformation()
	{

	}


public:
	bool OutputAccuracyFile(const char*); ///Show accuracy scripts (in some standard formats, xml? txt? it depends)

};

