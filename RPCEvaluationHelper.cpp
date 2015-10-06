#include "stdafx.h"
#include "RPCEvaluationHelper.h"


RPCEvaluationHelper::RPCEvaluationHelper(void)
{
	std::memset(_s, 0, sizeof(double)* 20);
	std::memset(_l, 0, sizeof(double)* 20);
	std::memset(_a, 0, sizeof(double)* 20);
	std::memset(_b, 0, sizeof(double)* 20);
	bScaled = false;
	bSwitched = false;
	m_bflag = true;
	ridgeconfigurationX = 0.000225; ridgeconfigurationY = 0.000225; //岭估计参数默认值
	_modelinited = false;
}

//string RPCEvaluationHelper::_donominatortype[4] = { "分母不相同", "分母相同", "分母为1", "出错" };
//string RPCEvaluationHelper::_processtype[5] = { "一次", "二次不迭代", "二次迭代", "三次不迭代", "三次迭代" };

int counter = 0;
#pragma omp threadprivate(counter)
int increment_counter()
{
	counter++;
	return(counter);
}

RPCEvaluationHelper::~RPCEvaluationHelper(void)
{
	PointList.clear();
	CPointList.clear();
}


bool RPCEvaluationHelper::LightIn(double lon, double lat, double h, double& x, double& y)
{
	double _P = (lat - np_array[2]) / np_array[3];
	double _L = (lon - np_array[0]) / np_array[1];
	double _H = (h - np_array[4]) / np_array[5];

	double para[20] = { 0.0 };
	para[0] = 1; para[1] = _L; para[2] = _P; para[3] = _H;
	para[4] = _L*_P; para[5] = _L * _H; para[6] = _P * _H;
	para[7] = _L * _L; para[8] = _P * _P; para[9] = _H * _H;
	para[10] = _P * _L * _H; para[11] = _L * _L * _L;
	para[12] = _L * _P * _P; para[13] = _L * _H * _H;
	para[14] = _L * _L * _P; para[15] = _P * _P * _P;
	para[16] = _P * _H * _H;
	para[17] = _L * _L * _H;
	para[18] = _P * _P * _H;
	para[19] = _H * _H * _H;

	double LineNum = 0, LineDen = 0, SampNum = 0, SampDen = 0;
	for (int i = 0; i < 20; i++)
	{
		LineNum += _s[i] * para[i];
		LineDen += _l[i] * para[i];

		SampNum += _a[i] * para[i];
		SampDen += _b[i] * para[i];
	}

	x = SampNum / SampDen * np_array[7] + np_array[6];
	y = LineNum / LineDen * np_array[9] + np_array[8];

	//double range = 1000;
	return IRange(x, y, 1000);
}

bool RPCEvaluationHelper::GetCommander(const char* ctrfilepath)
{
	//Open Control Point file
	ifstream reader;
	reader.open(ctrfilepath, ios::in);
	PointList.clear();
	//
	ifstream readproj;
	char readerproj[2048];
	readproj.open("projinfo.txt", ios::in);
	if (!readproj.is_open())
	{
		cout << "缺失投影信息文件，请将projinfo.txt文件放置到执行程序RPCModel.exe所在目录" << endl;
		cout << "若为调试模式，请将projinfo.txt放置到代码所在目录下" << endl;
		cout << "程序即将退出..." << endl;
		abort();
	}
	readproj.getline(readerproj, 2048);
	cout << "Infomation about reference map" << endl;
	cout << readerproj << endl;
	readproj.close();



	bool _isLonlat = false;
	do
	{
		char readerstring[1024];
		reader.getline(readerstring, 1024);

		cout << "开始读取";
		FPoint fp;
		//249              1365.115723         3731.536212       525588.389400      4768326.068400 
		int idd;
		sscanf_s(readerstring, "%d %lf %lf %lf %lf %lf", &idd, &fp.Pixel.x, &fp.Pixel.y, &fp.ground.X, &fp.ground.Y, &fp.ground.Z);
		//{
		//cout << "读取信息不正确" << endl;
		//break;
		//}
		sprintf_s(fp.id, "%d", idd);
		fp.Pixel.y -= 0;
		//cout << "读取正常" << endl;

		//转经纬度
		//DATUM WGS84[] = {6378137.0, 6356752.314, 0.08181919, 0.006694379, 0.0, 0.0, 0.0};
		//satOrbitHelper.rect2geograph(fp.ground.X,fp.ground.Y,fp.ground.Z,WGS84,&fp.ground.lat,&fp.ground.lon,&fp.ground.h);
		///如果看起来像经纬度，就不转了
		if (fabs(fp.ground.X)<180 && fabs(fp.ground.Y)<90)
		{
			_isLonlat = true;
			fp.ground.lon = fp.ground.X;
			fp.ground.lat = fp.ground.Y;
			fp.ground.X = -9999;
			fp.ground.Y = -9999;
		}
		else
		{
			double  lon, lat;
			//imgin.ConvertToWGS84(fp.ground.X, fp.ground.Y, readerproj, fp.ground.lon, fp.ground.lat);
		}

		//satOrbitHelper.MapXY2LatLon(fp.ground.X,fp.ground.Y,fp.ground.lat,fp.ground.lon);
		fp.ground.h = fp.ground.Z;
		//ofs1 << fp.ground.lon << "	" << fp.ground.lat << "	" << fp.ground.h << endl;

		PointList.push_back(fp);
		cout << fp << endl;

	} while (reader.peek() != EOF);

	if (_isLonlat)
	{
		cout << "控制点文件中出现经纬度，程序默认按WGS84经纬度处理，请核实！" << endl;
		cout << "当精度出现异常时，请留意此条信息" << endl;
	}

	else
	{
		cout << "控制点文件地面坐标为平面坐标，请确保exe所在路径下project.txt的投影信息与当前影像配套！" << endl;
		cout << "当程序精度出现异常时，请留意此条信息" << endl;
	}

	//	ofs1.close();
	reader.close(); //Adding,
	cout << "Control Point Got:" << PointList.size() << endl;
	return true;

}

bool RPCEvaluationHelper::ISWitch()
{
	if (!bScaled)
	{
		Compute_npara();
		bScaled = true;
	}

	for (int i = 0; i < PointList.size(); i++)
	{
		PointList[i].ground.lon = (!bSwitched) ? (PointList[i].ground.lon - np_array[0]) / (np_array[1]) : np_array[1] * PointList[i].ground.lon + np_array[0];
		PointList[i].ground.lat = (!bSwitched) ? (PointList[i].ground.lat - np_array[2]) / (np_array[3]) : np_array[3] * PointList[i].ground.lat + np_array[2];//补充
		PointList[i].ground.h = (!bSwitched) ? (PointList[i].ground.h - np_array[4]) / (np_array[5]) : np_array[5] * PointList[i].ground.h + np_array[4];//补充
		PointList[i].Pixel.x = (!bSwitched) ? (PointList[i].Pixel.x - np_array[6]) / (np_array[7]) : np_array[7] * PointList[i].Pixel.x + np_array[6];//补充
		PointList[i].Pixel.y = (!bSwitched) ? (PointList[i].Pixel.y - np_array[8]) / (np_array[9]) : np_array[9] * PointList[i].Pixel.y + np_array[8];//补充
	}

	for (int i = 0; i < CPointList.size(); i++)
	{
		CPointList[i].ground.lon = (!bSwitched) ? (CPointList[i].ground.lon - np_array[0]) / (np_array[1]) : np_array[1] * CPointList[i].ground.lon + np_array[0];
		CPointList[i].ground.lat = (!bSwitched) ? (CPointList[i].ground.lat - np_array[2]) / (np_array[3]) : np_array[3] * CPointList[i].ground.lat + np_array[2];//补充
		CPointList[i].ground.h = (!bSwitched) ? (CPointList[i].ground.h - np_array[4]) / (np_array[5]) : np_array[5] * CPointList[i].ground.h + np_array[4];//补充
		CPointList[i].Pixel.x = (!bSwitched) ? (CPointList[i].Pixel.x - np_array[6]) / (np_array[7]) : np_array[7] * CPointList[i].Pixel.x + np_array[6];//补充
		CPointList[i].Pixel.y = (!bSwitched) ? (CPointList[i].Pixel.y - np_array[8]) / (np_array[9]) : np_array[9] * CPointList[i].Pixel.y + np_array[8];//补充
	}

	for (int i = 0; i < _CPointList.size(); i++)
	{
		_CPointList[i].ground.lon = (!bSwitched) ? (_CPointList[i].ground.lon - np_array[0]) / (np_array[1]) : np_array[1] * _CPointList[i].ground.lon + np_array[0];
		_CPointList[i].ground.lat = (!bSwitched) ? (_CPointList[i].ground.lat - np_array[2]) / (np_array[3]) : np_array[3] * _CPointList[i].ground.lat + np_array[2];//补充
		_CPointList[i].ground.h = (!bSwitched) ? (_CPointList[i].ground.h - np_array[4]) / (np_array[5]) : np_array[5] * _CPointList[i].ground.h + np_array[4];//补充
		_CPointList[i].Pixel.x = (!bSwitched) ? (_CPointList[i].Pixel.x - np_array[6]) / (np_array[7]) : np_array[7] * _CPointList[i].Pixel.x + np_array[6];//补充
		_CPointList[i].Pixel.y = (!bSwitched) ? (_CPointList[i].Pixel.y - np_array[8]) / (np_array[9]) : np_array[9] * _CPointList[i].Pixel.y + np_array[8];//补充
	}

	bSwitched = !bSwitched;

	return bSwitched;  /// So the user of the class can know what's the situation
}

bool RPCEvaluationHelper::DoRPC1()
{
	if (!Switched())
		ISWitch();
	double a[7], aa[49], al[7];
	double bb[49], bl[7];

	std::memset(a, 0, sizeof(double)* 7);
	std::memset(aa, 0, sizeof(double)* 49);
	std::memset(al, 0, sizeof(double)* 7);

	std::memset(bb, 0, sizeof(double)* 49);
	std::memset(bl, 0, sizeof(double)* 7);

	for (int i = 0; i < PointList.size(); i++)
	{
		a[0] = 1.0; a[1] = PointList[i].ground.lon; a[2] = PointList[i].ground.lat; a[3] = PointList[i].ground.h;
		a[4] = -PointList[i].Pixel.y*PointList[i].ground.lon; a[5] = -PointList[i].Pixel.y*PointList[i].ground.lat;
		a[6] = -PointList[i].Pixel.y*PointList[i].ground.h;
		satOrbitHelper.pNormal(a, 7, PointList[i].Pixel.y, aa, al, 1.0);

		a[0] = 1.0; a[1] = PointList[i].ground.lon; a[2] = PointList[i].ground.lat; a[3] = PointList[i].ground.h;
		a[4] = -PointList[i].Pixel.x*PointList[i].ground.lon; a[5] = -PointList[i].Pixel.x*PointList[i].ground.lat;
		a[6] = -PointList[i].Pixel.x*PointList[i].ground.h;
		satOrbitHelper.pNormal(a, 7, PointList[i].Pixel.x, bb, bl, 1.0);
	}
	satOrbitHelper.Gauss(aa, al, 7);
	satOrbitHelper.Gauss(bb, bl, 7);

	std::memset(_s, 0, sizeof(double)* 20);
	std::memset(_l, 0, sizeof(double)* 20);
	std::memset(_a, 0, sizeof(double)* 20);
	for (int i = 0; i < 20; i++)
		_b[i] = 0.0;


	memcpy(_s, al, sizeof(double)* 4);
	memcpy(_l + 1, al + 4, sizeof(double)* 3);
	memcpy(_a, bl, sizeof(double)* 4);
	memcpy(_b + 1, bl + 4, sizeof(double)* 3);

	_l[0] = 1.0;
	_b[0] = 1.0;

	return true;
}

void RPCEvaluationHelper::Clear()
{
	std::memset(_s, 0, sizeof(double));
	std::memset(_l, 0, sizeof(double));
	std::memset(_a, 0, sizeof(double));
	std::memset(_b, 0, sizeof(double));
}

#pragma warning(suppress: 6262)
bool RPCEvaluationHelper::DoRPC3(bool bIteration)
{
	if (!Switched())
		ISWitch();

	Clear();

	if (bIteration)
		DoRPC2(true);

	double a[39], aa[1521], al[39];
	double bb[1521], bl[39];

	std::memset(a, 0, sizeof(double)* 39);
	std::memset(aa, 0, sizeof(double)* 1521);
	std::memset(al, 0, sizeof(double)* 39);
	std::memset(bb, 0, sizeof(double)* 1521);
	std::memset(bl, 0, sizeof(double)* 39);
	double altl = 0, bltl = 0;

	for (int i = 0; i < PointList.size(); i++)
	{
		a[0] = 1.0;
		a[1] = PointList[i].ground.lon;
		a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;
		a[4] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[5] = PointList[i].ground.lon * PointList[i].ground.h;
		a[6] = PointList[i].ground.lat * PointList[i].ground.h;
		a[7] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[8] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[9] = PointList[i].ground.h * PointList[i].ground.h;
		a[10] = PointList[i].ground.lat  * PointList[i].ground.lon * PointList[i].ground.h;
		a[11] = pow(PointList[i].ground.lon, 3);
		a[12] = PointList[i].ground.lon * pow(PointList[i].ground.lat, 2);
		a[13] = PointList[i].ground.lon * pow(PointList[i].ground.h, 2);
		a[14] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.lat;
		a[15] = pow(PointList[i].ground.lat, 3);
		a[16] = PointList[i].ground.lat * pow(PointList[i].ground.h, 2);
		a[17] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.h;
		a[18] = pow(PointList[i].ground.lat, 2) * PointList[i].ground.h;
		a[19] = pow(PointList[i].ground.h, 3);

		for (int j = 20; j < 39; j++)
		{
			a[j] = -a[j - 19] * PointList[i].Pixel.y;
		}

		satOrbitHelper.pNormal(a, 39, PointList[i].Pixel.y, aa, al, 1.0);

		altl += PointList[i].Pixel.y;
		a[0] = 1.0;
		a[1] = PointList[i].ground.lon;
		a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;
		a[4] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[5] = PointList[i].ground.lon * PointList[i].ground.h;
		a[6] = PointList[i].ground.lat * PointList[i].ground.h;
		a[7] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[8] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[9] = PointList[i].ground.h * PointList[i].ground.h;
		a[10] = PointList[i].ground.lat  * PointList[i].ground.lon * PointList[i].ground.h;
		a[11] = pow(PointList[i].ground.lon, 3);
		a[12] = PointList[i].ground.lon * pow(PointList[i].ground.lat, 2);
		a[13] = PointList[i].ground.lon * pow(PointList[i].ground.h, 2);
		a[14] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.lat;
		a[15] = pow(PointList[i].ground.lat, 3);
		a[16] = PointList[i].ground.lat * pow(PointList[i].ground.h, 2);
		a[17] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.h;
		a[18] = pow(PointList[i].ground.lat, 2) * PointList[i].ground.h;
		a[19] = pow(PointList[i].ground.h, 3);

		for (int j = 20; j < 39; j++)
		{
			a[j] = -a[j - 19] * PointList[i].Pixel.x;
		}

		satOrbitHelper.pNormal(a, 39, PointList[i].Pixel.x, bb, bl, 1.0);

		bltl += PointList[i].Pixel.x;
	}

	double x[39];

	memcpy(x, _s, sizeof(double)* 20);
	memcpy(x + 20, _l + 1, sizeof(double)* 19);

	for (int j = 0; j < 39; j++)
	{
		aa[j * 39 + j] += 0;
	}
	//satOrbitHelper.Gauss(aa,al,39);

	satOrbitHelper.GaussExt(aa, al, x, 39);
	//satOrbitHelper.GaussL(aa,al,x,39,altl/PointList.size());
	memcpy(al, x, sizeof(double)* 39);

	memcpy(x, _a, sizeof(double)* 20);
	memcpy(x + 20, _b + 1, sizeof(double)* 19);
	for (int j = 0; j < 39; j++)
		bb[j * 39 + j] += 0;
	//satOrbitHelper.Gauss(bb,bl,39);
	satOrbitHelper.GaussExt(bb, bl, x, 39);
	//satOrbitHelper.GaussL(bb,bl,x,39,bltl/PointList.size());
	memcpy(bl, x, sizeof(double)* 39);

	Clear();

	memcpy(_s, al, sizeof(double)* 20);
	memcpy(_l + 1, al + 20, sizeof(double)* 19);
	memcpy(_a, bl, sizeof(double)* 20);
	memcpy(_b + 1, bl + 20, sizeof(double)* 19);

	_l[0] = 1.0;
	_b[0] = 1.0;

	return true;
}

void RPCEvaluationHelper::INewRPC(vector<FPoint> _newpoints, string pathnew)
{
	_CPointList.clear();
	CPointList.clear();
	PointList.clear();

	_checkfile = pathnew;
	cout << "一共有:" << _newpoints.size() << "个点参与新的RPC解算" << endl;
	for (int i = 0; i < _newpoints.size(); i++)
	{
		PointList.push_back(_newpoints[i]);
	}

	ISWitch();
	IConductor(static_cast<RPCalculator>(15));
	ISWitch();
	OutPutRpc();
}

bool RPCEvaluationHelper::IFile(string _filepath, bool _filein, bool _breversed)
{
	if (_filein)
	{
		_modelinited = true;
		FILE *fp;
		fopen_s(&fp, _filepath.c_str(), "r");
		if (NULL == fp) {
			return false;
		}

		px[0] = 0; px[1] = 1; px[2] = 0;
		py[0] = 0; py[1] = 0; py[2] = 1;

		int i;
		char line[256];
		memset(line, 0, 256);
		//cout << "dfdf" << endl;

		if (_breversed == false)
		{
			fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256);
			np_array[8] = atof(line);
			fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array[6] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array[2] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array[0] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array[4] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array[9] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array[7] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array[3] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array[1] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array[5] = atof(line); fscanf_s(fp, "%s", line, 256);
		}
		else
		{
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array_Ind[2] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array_Ind[0] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array_Ind[6] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array_Ind[8] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array_Ind[4] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array_Ind[3] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array_Ind[1] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array_Ind[7] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array_Ind[9] = atof(line); fscanf_s(fp, "%s", line, 256);
			fscanf_s(fp, "%s", line, 256); fscanf_s(fp, "%s", line, 256); np_array_Ind[5] = atof(line); fscanf_s(fp, "%s", line, 256);
		}

		if (_breversed == false)
		{
			for (i = 0; i < 20; i++)
			{
				fscanf_s(fp, "%s", line, 256);
				fscanf_s(fp, "%s", line, 256);
				_s[i] = atof(line);
			}

			for (i = 0; i < 20; i++)
			{
				fscanf_s(fp, "%s", line, 256);
				fscanf_s(fp, "%s", line, 256);
				_l[i] = atof(line);
			}

			for (i = 0; i < 20; i++)
			{
				fscanf_s(fp, "%s", line, 256);
				fscanf_s(fp, "%s", line, 256);
				_a[i] = atof(line);
			}

			for (i = 0; i < 20; i++)
			{
				fscanf_s(fp, "%s", line, 256);
				fscanf_s(fp, "%s", line, 256);
				_b[i] = atof(line);
			}
		}

		else
		{
			for (i = 0; i < 20; i++)
			{
				fscanf_s(fp, "%s", line, 256);
				fscanf_s(fp, "%s", line, 256);
				_ds[i] = atof(line);
			}

			for (i = 0; i < 20; i++)
			{
				fscanf_s(fp, "%s", line, 256);
				fscanf_s(fp, "%s", line, 256);
				_dl[i] = atof(line);
			}

			for (i = 0; i < 20; i++)
			{
				fscanf_s(fp, "%s", line, 256);
				fscanf_s(fp, "%s", line, 256);
				_da[i] = atof(line);
			}

			for (i = 0; i < 20; i++)
			{
				fscanf_s(fp, "%s", line, 256);
				fscanf_s(fp, "%s", line, 256);
				_db[i] = atof(line);
			}

			_bInded = true;
		}
		fclose(fp);

		///rpc正反算的精度，在研发阶段，这种步骤是必须的
		//ofstream ofs;
		//ofs.open("D:\\RPCheck.txt", ios::out);
		int heightnow = _imgheight; int widthnow = _imgwidth;
#ifdef _DEBUG
		if (_breversed)
		{
			cout << "引入了反算RPC文件，将检查RPC模型正反算精度--要求：最大值小于0.001个像素" << endl;
			for (int t = 0; t < heightnow; t++)
			{
				//if (t == 0) continue;
				bool badline = false;
#pragma omp parallel for
				for (int s = 0; s < widthnow; s++)
				{
					try
					{

						double lontest, latest;
						double xcheck, ycheck;
						LightOut(s, t, 100.0, lontest, latest);
						LightIn(lontest, latest, 100.0, xcheck, ycheck);
						//ofs << s << "	" << t << "	" << xcheck << "	" << ycheck <<"	"<<xcheck-s<<""<<ycheck - t<< endl;
						if (!_P_F_AccCheck(xcheck, s, ycheck, t, lontest, latest))
						{
							badline = true;
							continue;
						}
					}
					catch (exception const&e)
					{
						cout << "an error of type:" << typeid(e).name() << "occured" << endl;
						cout << e.what() << endl;
					}
				}
				if (badline)
					cout << t << endl;
			}
			cout << "正反算精度合格！(不高于0.01个像素)" << endl;

			cout << "继续?" << endl;
			int modedebug = 0;
			cin >> modedebug;
		}
#endif


		/*
		if (NULL == strImgRangeFile)
		{
		DPoint2D ver[4];
		ver[0].x = ver[3].x = SAMP_OFF - SAMP_SCALE;
		ver[1].x = ver[2].x = SAMP_OFF + SAMP_SCALE;
		ver[0].y = ver[1].y = LINE_OFF - LINE_SCALE;
		ver[2].y = ver[3].y = LINE_OFF + LINE_SCALE;
		m_ModelRange = DPolygon(ver, 4);
		}
		else
		{
		FILE* fp;
		fopen_s(&fp, strImgRangeFile, "r");
		unsigned num;
		fscanf_s(fp, "%d", &num);

		DPoint2D* po = new DPoint2D[num];
		int x, y;
		for (unsigned i = 0; i<num; ++i)
		{
		fscanf_s(fp, "%d %d", &x, &y);
		po[i].x = x, po[i].y = y;
		}
		fclose(fp);
		DPolygon DPo(po, static_cast<unsigned> (num));
		delete[]po;
		m_ModelRange = DPo;
		}
		*/
		return true;
	}
}

bool RPCEvaluationHelper::DoRPC2(bool bIteration)
{
	if (!Switched())
		ISWitch();

	Clear();
	if (bIteration)
		DoRPC1();

	double a[19];
	double aa[361], al[19];
	double bb[361], bl[19];
	double altl = 0, bltl = 0;

	std::memset(a, 0, sizeof(double)* 19);
	std::memset(aa, 0, sizeof(double)* 361);
	std::memset(al, 0, sizeof(double)* 19);

	std::memset(bb, 0, sizeof(double)* 361);
	std::memset(bl, 0, sizeof(double)* 19);



	for (int i = 0; i < PointList.size(); i++)
	{
		a[0] = 1.0;
		a[1] = PointList[i].ground.lon;
		a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;
		a[4] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[5] = PointList[i].ground.lon * PointList[i].ground.h;
		a[6] = PointList[i].ground.lat * PointList[i].ground.h;
		a[7] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[8] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[9] = PointList[i].ground.h * PointList[i].ground.h;

		for (int j = 10; j <= 18; j++)
		{
			a[j] = -a[j - 9] * PointList[i].Pixel.y;
		}

		satOrbitHelper.pNormal(a, 19, PointList[i].Pixel.y, aa, al, 1.0);

		altl += PointList[i].Pixel.y;

		a[0] = 1.0;
		a[1] = PointList[i].ground.lon;
		a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;
		a[4] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[5] = PointList[i].ground.lon * PointList[i].ground.h;
		a[6] = PointList[i].ground.lat * PointList[i].ground.h;
		a[7] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[8] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[9] = PointList[i].ground.h * PointList[i].ground.h;

		for (int j = 10; j <= 18; j++)
		{
			a[j] = -a[j - 9] * PointList[i].Pixel.x;
		}

		satOrbitHelper.pNormal(a, 19, PointList[i].Pixel.x, bb, bl, 1.0);

		bltl += PointList[i].Pixel.x;
	}

	double x[39];
	memcpy(x, _s, sizeof(double)* 10);
	memcpy(x + 10, _l + 1, sizeof(double)* 9);
	satOrbitHelper.GaussExt(aa, al, x, 19);
	memcpy(al, x, sizeof(double)* 19);

	memcpy(x, _a, sizeof(double)* 10);
	memcpy(x + 10, _b + 1, sizeof(double)* 9);
	satOrbitHelper.GaussExt(bb, bl, x, 19);
	memcpy(bl, x, sizeof(double)* 19);

	std::memset(_s, 0, sizeof(double)* 20);
	std::memset(_l, 0, sizeof(double)* 20);
	std::memset(_a, 0, sizeof(double)* 20);
	std::memset(_b, 0, sizeof(double)* 20);

	memcpy(_s, al, sizeof(double)* 10);
	memcpy(_l + 1, al + 10, sizeof(double)* 9);
	memcpy(_a, bl, sizeof(double)* 10);
	memcpy(_b + 1, bl + 10, sizeof(double)* 9);

	_l[0] = 1.0;
	_b[0] = 1.0;

	return true;
}

double RPCEvaluationHelper::distanceOnePixel(double latitude, double longitude)
{
	double x1, y1, x2, y2;
	LightIn(latitude, longitude, np_array[4], x1, y1);
	LightIn(latitude + 0.01* np_array[3], longitude, np_array[4], x2, y2);

	return  0.01*np_array[3] / sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

void RPCEvaluationHelper::LightOut(double x, double y, double h, double& lon, double&lat)
{
	if (false)
	{
		double _P = (y - np_array_Ind[2]) / np_array_Ind[3];
		double _L = (x - np_array_Ind[0]) / np_array_Ind[1];
		double _H = (h - np_array_Ind[4]) / np_array_Ind[5];

		double para[20] = { 0.0 };
		para[0] = 1; para[1] = _L; para[2] = _P; para[3] = _H;
		para[4] = _L*_P; para[5] = _L * _H; para[6] = _P * _H;
		para[7] = _L * _L; para[8] = _P * _P; para[9] = _H * _H;
		para[10] = _P * _L * _H; para[11] = _L * _L * _L;
		para[12] = _L * _P * _P; para[13] = _L * _H * _H;
		para[14] = _L * _L * _P; para[15] = _P * _P * _P;
		para[16] = _P * _H * _H;
		para[17] = _L * _L * _H;
		para[18] = _P * _P * _H;
		para[19] = _H * _H * _H;

		double LineNum = 0, LineDen = 0, SampNum = 0, SampDen = 0;
		for (int i = 0; i < 20; i++)
		{
			LineNum += _ds[i] * para[i];
			LineDen += _dl[i] * para[i];

			SampNum += _da[i] * para[i];
			SampDen += _db[i] * para[i];
		}

		lon = SampNum / SampDen * np_array_Ind[9] + np_array_Ind[8];
		lat = LineNum / LineDen * np_array_Ind[7] + np_array_Ind[6];
	}
	else
	{
		lat = np_array[2];
		lon = np_array[0];
		PreciseBasedAffine(lat, lon, np_array[3] * 0.5, np_array[1] * 0.5, x, y, h);

		double xp, yp;
		LightIn(lon, lat, h, xp, yp);
		double e = (xp - x)*(xp - x) + (yp - y)*(yp - y);
		double e0;
		do {
			e0 = e;
			double pixelsize = distanceOnePixel(lat, lon);

			LightIn(lon, lat, h, xp, yp);

			PreciseBasedAffine(lat, lon, fabs(yp - y)*pixelsize, fabs(xp - x)*pixelsize,
				x, y, h);
			LightIn(lon, lat, h, xp, yp);

			e = (xp - x)*(xp - x) + (yp - y)*(yp - y);

			if (e < 0.000001) break;

		} while (e < e0);
	}
}

void RPCEvaluationHelper::PreciseBasedAffine(double &latitude, double &longitude, double dlat, double dlon, double x0, double y0, double h0)
{
	//int i = 0;
	double s[5], l[5], lat[5], lon[5];
	lat[0] = latitude; lon[0] = longitude;
	lat[1] = latitude + dlat; lon[1] = longitude + dlon;
	lat[2] = latitude + dlat; lon[2] = longitude - dlon;
	lat[3] = latitude - dlat; lon[3] = longitude + dlon;
	lat[4] = latitude - dlat; lon[4] = longitude - dlon;

	for (int i = 0; i < 5; i++)
	{
		LightIn(lon[i], lat[i], h0, s[i], l[i]);
	}
	CSatOrbit satobhelper;
	double LAT_OFF, LONG_OFF, LINE_OFF, SAMP_OFF, LAT_SCALE, LONG_SCALE, LINE_SCALE, SAMP_SCALE;
	satobhelper.Compute_avAnddx(lat, 5, LAT_OFF, LAT_SCALE);
	satobhelper.Compute_avAnddx(lon, 5, LONG_OFF, LONG_SCALE);
	satobhelper.Compute_avAnddx(l, 5, LINE_OFF, LINE_SCALE);
	satobhelper.Compute_avAnddx(s, 5, SAMP_OFF, SAMP_SCALE);

	//归一化
	for (int i = 0; i < 5; i++)
	{
		l[i] = (l[i] - LINE_OFF) / LINE_SCALE;
		s[i] = (s[i] - SAMP_OFF) / SAMP_SCALE;

		lat[i] = (lat[i] - LAT_OFF) / LAT_SCALE;
		lon[i] = (lon[i] - LONG_OFF) / LONG_SCALE;
	}

	double ps[3], pl[3];
	//P=ps[0]+ps[1]*s+ps[2]*l;
	//L=pl[0]+pl[1]*s+pl[2]*l;
	double aas[9], aal[9], a[3];
	memset(aas, 0, sizeof(double)* 9);
	memset(aal, 0, sizeof(double)* 9);
	memset(ps, 0, sizeof(double)* 3);
	memset(pl, 0, sizeof(double)* 3);
	for (int i = 0; i < 5; i++)
	{
		a[0] = 1, a[1] = s[i], a[2] = l[i];
		satobhelper.pNormal(a, 3, lat[i], aas, ps, 1);
		satobhelper.pNormal(a, 3, lon[i], aal, pl, 1);
	}
	satobhelper.Gauss(aas, ps, 3);
	satobhelper.Gauss(aal, pl, 3);


	y0 = (y0 - LINE_OFF) / LINE_SCALE;
	x0 = (x0 - SAMP_OFF) / SAMP_SCALE;

	double P = ps[0] + ps[1] * x0 + ps[2] * y0;
	double L = pl[0] + pl[1] * x0 + pl[2] * y0;

	latitude = P*LAT_SCALE + LAT_OFF;
	longitude = L*LONG_SCALE + LONG_OFF;
}

void RPCEvaluationHelper::Checker()
{
	if (Switched())
		ISWitch();

	ofstream ofs1;
	ofs1.open("information.txt");
	ofstream ofs;
	char reader[256];
	string basefile = "";
	if (_currentcalculator == Linear_POW3_AllIteration)
		basefile = _checkfile + "整体迭代" + ".txt";
	else basefile = _checkfile + "_Acc_";// +_donominatortype[(static_cast<int>(_currentcalculator) / 5)] + _processtype[(static_cast<int>(_currentcalculator) % 5)] + ".txt";
	ofs.open(basefile.c_str(), ios::out);

	ofs << "控制点个数：" << PointList.size() << endl;
	ofs << "控制点残差：" << endl;
	ofs << "PointID" << "\t" << "Longitude" << "\t" << "Latitude" << "\t" << "Height" << "\t" << "ProJectionX" << "\t" << "ProjectionY" << "\t" << "ProjectionZ" << "\t" << "PixelX" << "\t" << "PixelY" << "\t" << "XError" << "\t" << "YError" << endl;
	ofs.precision(15);
	double maxerrory; double maxerrorx; maxerrorx = maxerrory = 0;
	double minerrory, minerrorx; minerrorx = minerrory = 99.0;
	char minxerrorindex[256], maxxerrorindex[256], minyerrorindex[256], maxyerrorindex[256];
	//minxerrorindex = maxxerrorindex = minyerrorindex = maxyerrorindex = 0.0;
	double rmsx, rmsy;
	rmsx = rmsy = 0.0;

	if (Switched()) ISWitch();
	for (int i = 0; i < PointList.size(); i++)
	{
		double _x, _y;
		LightIn(PointList[i].ground.lon, PointList[i].ground.lat, PointList[i].ground.h,
			_x, _y);

		ofs << PointList[i].id << "\t" << PointList[i].ground.lon << "\t" << PointList[i].ground.lat << "\t"
			<< PointList[i].ground.h << "\t" << PointList[i].ground.X << "\t" << PointList[i].ground.Y << "\t" << PointList[i].ground.Z
			<< "\t" << PointList[i].Pixel.x << "\t" << PointList[i].Pixel.y << "\t"
			<< _x - PointList[i].Pixel.x << "\t" << _y - PointList[i].Pixel.y << endl;

		rmsx += pow(_x - PointList[i].Pixel.x, 2) / static_cast<double>(PointList.size() - 1);
		rmsy += pow(_y - PointList[i].Pixel.y, 2) / static_cast<double>(PointList.size() - 1);


		//计算最大最小值
		if (fabs(_x - PointList[i].Pixel.x) > fabs(maxerrorx))
		{
			maxerrorx = _x - PointList[i].Pixel.x;
			//maxxerrorindex = PointList[i].id;
			memcpy(maxxerrorindex, PointList[i].id, 256);
		}

		if (fabs(_x - PointList[i].Pixel.x)<fabs(minerrorx))
		{
			minerrorx = _x - PointList[i].Pixel.x;
			//minxerrorindex = PointList[i].id;
			memcpy(minxerrorindex, PointList[i].id, 256);
		}

		if (fabs(_y - PointList[i].Pixel.y) > fabs(maxerrory))
		{
			maxerrory = _y - PointList[i].Pixel.y;
			//maxyerrorindex = PointList[i].id;
			memcpy(maxyerrorindex, PointList[i].id, 256);
		}

		if (fabs(_y - PointList[i].Pixel.y) < fabs(minerrory))
		{
			minerrory = _y - PointList[i].Pixel.y;
			//minyerrorindex = PointList[i].id;
			memcpy(minyerrorindex, PointList[i].id, 256);
		}
	}

	rmsx = sqrt(rmsx);
	rmsy = sqrt(rmsy);

	ofs << "控制点像方残差" << endl;
	ofs << "x方向最大值：" << "\t" << maxerrorx << "\t" << "最小值" << minerrory << endl;
	ofs << "X方向最大点信息(点号):" << "\t" << maxxerrorindex << endl;
	ofs << "X方向最小点信息(点号):" << "\t" << minxerrorindex << endl;
	ofs << "y方向最大值：" << "\t" << maxerrory << "\t" << "最小值" << minerrory << endl;
	ofs << "y方向最大点信息(点号):" << "\t" << maxyerrorindex << endl;
	ofs << "y方向最大点信息(点号):" << "\t" << minyerrorindex << endl;
	ofs << "中误差" << "\t" << "X方向" << "\t" << "Y方向" << "\t" << "平面" << endl;
	ofs << rmsx << "\t" << rmsy << "\t" << "\t" << sqrt(rmsx * rmsx + rmsy * rmsy) << endl;

	cout.precision(15);
	cout << "控制点信息:" << endl;
	cout << "中误差" << "\t" << "X方向" << "\t" << "Y方向" << "\t" << "平面" << endl;
	cout << rmsx << "\t" << rmsy << "\t" << "\t" << sqrt(rmsx * rmsx + rmsy * rmsy) << endl;

	rmsx = rmsy = 0.0;
	memset(maxxerrorindex, 0, 256);
	memset(minxerrorindex, 0, 256);
	memset(maxyerrorindex, 0, 256);
	memset(minyerrorindex, 0, 256);
	maxerrorx = maxerrory = 0.0;
	minerrorx = minerrory = 99.0;

	ofs << "检查点个数:" << CPointList.size() << endl;
	ofs << "检查点残差:" << endl;
	ofs << "PointID" << "\t" << "Longitude" << "\t" << "Latitude" << "\t" << "Height" << "\t" << "ProJectionX" << "\t" << "ProjectionY" << "\t" << "ProjectionZ" << "\t" << "PixelX" << "\t" << "PixelY" << "\t" << "XError" << "\t" << "YError" << endl;
	if (Switched()) ISWitch();
	for (int i = 0; i < CPointList.size(); i++)
	{
		double _x, _y;
		LightIn(CPointList[i].ground.lon, CPointList[i].ground.lat, CPointList[i].ground.h,
			_x, _y);

		ofs << CPointList[i].id << "\t" << CPointList[i].ground.lon << "\t" << CPointList[i].ground.lat << "\t"
			<< CPointList[i].ground.h << "\t" << CPointList[i].ground.X << "\t" << CPointList[i].ground.Y << "\t"
			<< CPointList[i].ground.Z << "\t" << CPointList[i].Pixel.x << "\t" << CPointList[i].Pixel.y << "\t"
			<< _x - CPointList[i].Pixel.x << "\t" << _y - CPointList[i].Pixel.y << endl;

		rmsx += pow(_x - CPointList[i].Pixel.x, 2) / static_cast<double>(CPointList.size() - 1);
		rmsy += pow(_y - CPointList[i].Pixel.y, 2) / static_cast<double>(CPointList.size() - 1);


		//计算最大最小值
		if (fabs(_x - CPointList[i].Pixel.x) > fabs(maxerrorx))
		{
			maxerrorx = _x - CPointList[i].Pixel.x;
			//maxxerrorindex = CPointList[i].id;
			memcpy(maxxerrorindex, CPointList[i].id, 256);
		}

		if (fabs(_x - CPointList[i].Pixel.x)<fabs(minerrorx))
		{
			minerrorx = _x - CPointList[i].Pixel.x;
			//minxerrorindex = CPointList[i].id;
			memcpy(minxerrorindex, CPointList[i].id, 256);
		}

		if (fabs(_y - CPointList[i].Pixel.y) > fabs(maxerrory))
		{
			maxerrory = _y - CPointList[i].Pixel.y;
			//maxyerrorindex = CPointList[i].id;
			memcpy(maxyerrorindex, CPointList[i].id, 256);
		}

		if (fabs(_y - CPointList[i].Pixel.y) < fabs(minerrory))
		{
			minerrory = _y - CPointList[i].Pixel.y;
			//minyerrorindex = CPointList[i].id;
			memcpy(minyerrorindex, CPointList[i].id, 256);
		}
	}

	rmsx = sqrt(rmsx);
	rmsy = sqrt(rmsy);

	ofs << "检查点像方残差" << endl;
	ofs << "x方向最大值：" << "\t" << maxerrorx << "\t" << "最小值" << minerrory << endl;
	ofs << "X方向最大点信息(点号):" << "\t" << maxxerrorindex << endl;
	ofs << "X方向最小点信息(点号):" << "\t" << minxerrorindex << endl;
	ofs << "y方向最大值：" << "\t" << maxerrory << "\t" << "最小值" << minerrory << endl;
	ofs << "y方向最大点信息(点号):" << "\t" << maxyerrorindex << endl;
	ofs << "y方向最大点信息(点号):" << "\t" << minyerrorindex << endl;
	ofs << "中误差" << "\t" << "X方向" << "\t" << "Y方向" << "\t" << "平面" << endl;
	ofs << rmsx << "\t" << rmsy << "\t" << "\t" << sqrt(rmsx * rmsx + rmsy * rmsy) << endl;

	cout.precision(15);
	cout << "检查点信息:" << endl;
	cout << "中误差" << "\t" << "X方向" << "\t" << "Y方向" << "\t" << "平面" << endl;
	cout << rmsx << "\t" << rmsy << "\t" << "\t" << sqrt(rmsx * rmsx + rmsy * rmsy) << endl;

	ofs.close();
}

bool RPCEvaluationHelper::DoRPC1TrueSame()
{
	if (!Switched())
	{
		ISWitch();
		cout << "Warning:Something is wrong with the swithcer.It has been corrected." << endl;
	}
	Clear();
	double a[11], aa[121], al[11];

	std::memset(a, 0, sizeof(double)* 11);
	std::memset(aa, 0, sizeof(double)* 121);
	std::memset(al, 0, sizeof(double)* 11);

	for (int i = 0; i < PointList.size(); i++)
	{
		std::memset(a, 0, sizeof(double)* 11);

		a[0] = 1.0; a[1] = PointList[i].ground.lon; a[2] = PointList[i].ground.lat; a[3] = PointList[i].ground.h;
		a[4] = -PointList[i].Pixel.y*PointList[i].ground.lon; a[5] = -PointList[i].Pixel.y*PointList[i].ground.lat;
		a[6] = -PointList[i].Pixel.y*PointList[i].ground.h;

		satOrbitHelper.pNormal(a, 11, PointList[i].Pixel.y, aa, al, 1.0);

		std::memset(a, 0, sizeof(double)* 11);
		a[7] = 1.0; a[8] = PointList[i].ground.lon; a[9] = PointList[i].ground.lat; a[10] = PointList[i].ground.h;
		a[4] = -PointList[i].Pixel.x*PointList[i].ground.lon; a[5] = -PointList[i].Pixel.x*PointList[i].ground.lat;
		a[6] = -PointList[i].Pixel.x*PointList[i].ground.h;
		satOrbitHelper.pNormal(a, 11, PointList[i].Pixel.x, aa, al, 1.0);
	}

	satOrbitHelper.Gauss(aa, al, 11);
	Clear();

	memcpy(_s, al, sizeof(double)* 4);
	memcpy(_l, al + 4, sizeof(double)* 3);
	memcpy(_a, al + 7, sizeof(double)* 4);
	memcpy(_b, al + 4, sizeof(double)* 3);

	_l[0] = 1.0;
	_b[0] = 1.0;


	return true;
}

bool RPCEvaluationHelper::DoRPC2TrueSame(bool iteration)
{
	if (!Switched())
	{
		ISWitch();
		cout << "Some thing is wrong with the plan watch out!" << endl;
	}

	if (iteration)
		DoRPC1TrueSame();

	double a[29], aa[841], al[29];

	std::memset(a, 0, sizeof(double)* 29);
	std::memset(aa, 0, sizeof(double)* 841);
	std::memset(al, 0, sizeof(double)* 29);

	for (int i = 0; i < PointList.size(); i++)
	{
		std::memset(a, 0, sizeof(double)* 29);

		a[0] = 1.0;
		a[1] = PointList[i].ground.lon;
		a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;
		a[4] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[5] = PointList[i].ground.lon * PointList[i].ground.h;
		a[6] = PointList[i].ground.lat * PointList[i].ground.h;
		a[7] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[8] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[9] = PointList[i].ground.h * PointList[i].ground.h;

		for (int j = 10; j <= 18; j++)
		{
			a[j] = -a[j - 9] * PointList[i].Pixel.y;
		}

		satOrbitHelper.pNormal(a, 29, PointList[i].Pixel.y, aa, al, 1.0);
		std::memset(a, 0, sizeof(double)* 29);

		a[19] = 1.0;
		a[20] = PointList[i].ground.lon;
		a[21] = PointList[i].ground.lat;
		a[22] = PointList[i].ground.h;
		a[23] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[24] = PointList[i].ground.lon * PointList[i].ground.h;
		a[25] = PointList[i].ground.lat * PointList[i].ground.h;
		a[26] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[27] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[28] = PointList[i].ground.h * PointList[i].ground.h;

		for (int j = 10; j <= 18; j++)
		{
			a[j] = -a[j + 10] * PointList[i].Pixel.x;
		}

		satOrbitHelper.pNormal(a, 29, PointList[i].Pixel.x, aa, al, 1.0);
	}

	cout << "pNormal Finished" << endl;

	double x[29];

	memcpy(x, _s, sizeof(double)* 10);
	memcpy(x + 10, _l + 1, sizeof(double)* 9);
	memcpy(x + 19, _a, sizeof(double)* 10);
	memcpy(x + 10, _b + 1, sizeof(double)* 9);

	satOrbitHelper.GaussExt(aa, al, x, 29);
	memcpy(al, x, sizeof(double)* 29);

	memcpy(_s, al, sizeof(double)* 10);
	memcpy(_l + 1, al + 10, sizeof(double)* 9);
	memcpy(_a, al + 19, sizeof(double)* 10);
	memcpy(_b + 1, al + 10, sizeof(double)* 9);

	_l[0] = 1.0;
	_b[0] = 1.0;

	return true;
}

bool RPCEvaluationHelper::DoRPC3TrueSame(bool iteration)
{
	if (!Switched())
	{
		ISWitch();
		cout << "Siutation should be noticed that. Switcher is wrong" << endl;
	}
	if (true)
	{
		DoRPC2TrueSame(true);
	}

	double a[59] = { 0.0 };
	double aa[3481] = { 0.0 };
	double al[59] = { 0.0 };


	for (int i = 0; i < PointList.size(); i++)
	{
		std::memset(a, 0.0, sizeof(double)* 59);
		a[0] = 1.0;
		a[1] = PointList[i].ground.lon;
		a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;
		a[4] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[5] = PointList[i].ground.lon * PointList[i].ground.h;
		a[6] = PointList[i].ground.lat * PointList[i].ground.h;
		a[7] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[8] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[9] = PointList[i].ground.h * PointList[i].ground.h;
		a[10] = PointList[i].ground.lat  * PointList[i].ground.lon * PointList[i].ground.h;
		a[11] = pow(PointList[i].ground.lon, 3);
		a[12] = PointList[i].ground.lon * pow(PointList[i].ground.lat, 2);
		a[13] = PointList[i].ground.lon * pow(PointList[i].ground.h, 2);
		a[14] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.lat;
		a[15] = pow(PointList[i].ground.lat, 3);
		a[16] = PointList[i].ground.lat * pow(PointList[i].ground.h, 2);
		a[17] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.h;
		a[18] = pow(PointList[i].ground.lat, 2) * PointList[i].ground.h;
		a[19] = pow(PointList[i].ground.h, 3);

		for (int j = 20; j <= 38; j++)
		{
			a[j] = -a[j - 19] * PointList[i].Pixel.y;
		}

		satOrbitHelper.pNormal(a, 59, PointList[i].Pixel.y, aa, al, 10.0);

		std::memset(a, 0.0, sizeof(double)* 59);

		a[39] = 1.0;
		a[40] = PointList[i].ground.lon;
		a[41] = PointList[i].ground.lat;
		a[42] = PointList[i].ground.h;
		a[43] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[44] = PointList[i].ground.lon * PointList[i].ground.h;
		a[45] = PointList[i].ground.lat * PointList[i].ground.h;
		a[46] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[47] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[48] = PointList[i].ground.h * PointList[i].ground.h;
		a[49] = PointList[i].ground.lat  * PointList[i].ground.lon * PointList[i].ground.h;
		a[50] = pow(PointList[i].ground.lon, 3);
		a[51] = PointList[i].ground.lon * pow(PointList[i].ground.lat, 2);
		a[52] = PointList[i].ground.lon * pow(PointList[i].ground.h, 2);
		a[53] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.lat;
		a[54] = pow(PointList[i].ground.lat, 3);
		a[55] = PointList[i].ground.lat * pow(PointList[i].ground.h, 2);
		a[56] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.h;
		a[57] = pow(PointList[i].ground.lat, 2) * PointList[i].ground.h;
		a[58] = pow(PointList[i].ground.h, 3);

		for (int j = 20; j <= 38; j++)
		{
			a[j] = -a[j + 20] * PointList[i].Pixel.x;
		}

		satOrbitHelper.pNormal(a, 59, PointList[i].Pixel.x, aa, al, 10.0);
	}

	double x[59];
	memcpy(x, _s, sizeof(double)* 20);
	memcpy(x + 20, _l + 1, sizeof(double)* 19);
	memcpy(x + 39, _a, sizeof(double)* 20);
	memcpy(x + 20, _b + 1, sizeof(double)* 19);

	satOrbitHelper.GaussExt(aa, al, x, 59);
	memcpy(al, x, sizeof(double)* 59);

	memcpy(_s, al, sizeof(double)* 20);
	memcpy(_l + 1, al + 20, sizeof(double)* 19);
	memcpy(_a, al + 39, sizeof(double)* 20);
	memcpy(_b, al + 20, sizeof(double)* 19);

	_l[0] = 1.0;
	_b[0] = 1.0;

	return true;
}

void RPCEvaluationHelper::RPCPLH2a(double P, double L, double H, double *a)
{
	a[0] = 1; a[1] = L; a[2] = P; a[3] = H; a[4] = L*P; a[5] = L*H; a[6] = P*H; a[7] = L*L; a[8] = P*P; a[9] = H*H;
	a[10] = P*L*H; a[11] = L*L*L; a[12] = L*P*P; a[13] = L*H*H; a[14] = L*L*P; a[15] = P*P*P; a[16] = P*H*H;
	a[17] = L*L*H; a[18] = P*P*H; a[19] = H*H*H;
}

void RPCEvaluationHelper::TranslateMark(int* mark, int* marks, int* markl, int* marka, int* markb)
{
	/*
	memcpy(marks, mark, sizeof(double)* 20);
	memcpy(markl, mark + 20,sizeof(double)* 20);
	memcpy(marka, mark + 40,sizeof(double)* 20);
	memcpy(markb, mark + 60,sizeof(double)* 20);
	*/
	for (int i = 0; i < 20; i++)
	{
		marks[i] = mark[i];
		markl[i] = mark[i + 20];
		marka[i] = mark[i + 40];
		markb[i] = mark[i + 60];
	}
}

void RPCEvaluationHelper::RPCModel_Selected(int* coeffmatrix)
{
	double aext[20], NumL, DenL, NumS, DenS, P, L, H;
	//double x[40];
	int marks[20], markl[20], marka[20], markb[20];
	TranslateMark(coeffmatrix, marks, markl, marka, markb);

	int sizeofpoint = PointList.size();
	int coeffY[40];
	int coeffX[40];
	double* coeffmatrix4X = new double[40 * sizeofpoint];//偏导项
	double* coeffmatrix4Y = new double[40 * sizeofpoint];
	memcpy(coeffY, coeffmatrix, sizeof(int)* 40);
	memcpy(coeffX, coeffmatrix + 40, sizeof(int)* 40);
	//std::memset(coeffmatrix4X, 0, sizeof(double)* 40 * sizeofpoint);
	//std::memset(coeffmatrix4Y, 0, sizeof(double)* 40 * sizeofpoint);

	for (int i = 0; i < sizeofpoint; i++)
	{
		P = PointList[i].ground.lat;
		L = PointList[i].ground.lon;
		H = PointList[i].ground.h;

		RPCPLH2a(P, L, H, aext);
		NumL = satOrbitHelper.RPCPLHDynamic(_s, P, L, H, marks);
		DenL = satOrbitHelper.RPCPLHDynamic(_l, P, L, H, markl);
		NumS = satOrbitHelper.RPCPLHDynamic(_a, P, L, H, marka);
		DenS = satOrbitHelper.RPCPLHDynamic(_b, P, L, H, markb);

		for (int j = 0; j < 20; j++)
		{
			coeffmatrix4Y[i * 40 + j] = aext[j] / DenL;
			coeffmatrix4X[i * 40 + j] = aext[j] / DenS;
		}

		coeffmatrix4Y[i * 40 + 20] = 0;
		coeffmatrix4X[i * 40 + 20] = 0;

		for (int j = 21; j < 40; j++)
		{
			coeffmatrix4Y[i * 40 + j] = -aext[j - 20] * NumL / DenL / DenL;
			coeffmatrix4X[i * 40 + j] = -aext[j - 20] * NumS / DenS / DenS;
		}
	}

	double cosMx = 0, cosMy = 0;
	double cosD1x = 0, cosD2x = 0;
	double cosD1y = 0, cosD2y = 0;
	double cos_corx = 0, cos_cory = 0;
	double maxrelationx, maxrelationy;
	maxrelationx = maxrelationy = -99;
	for (int i = 0; i < 40; i++)//for X
	{
		for (int j = 0; j < 40; j++)
		{
			if (i == 20 || j == 20) continue;
			if (i == j)
				continue;  ///相同列之间不做检查 //continue 作用：跳出本层循环，本层循环余下代码在本次循环中不执行
			///已经被评判为相关的不检查
			if (coeffX[i] == 0 || coeffX[j] == 0)
				continue;
			cos_corx = 0; cosMx = 0;
			cosD2x = cosD1x = 0;
			for (int t = 0; t < sizeofpoint; t++)
			{
				cosMx += coeffmatrix4X[i + t * 40] * coeffmatrix4X[j + t * 40];
				cosD1x += (coeffmatrix4X[i + t * 40] * coeffmatrix4X[i + t * 40]);
				cosD2x += (coeffmatrix4X[j + t * 40] * coeffmatrix4X[j + t * 40]);
			}
			cos_corx = cosMx / sqrt(cosD1x) / sqrt(cosD2x);
			if (fabs(cos_corx) > maxrelationx) maxrelationx = fabs(cos_corx);
			if (cos_corx > 0.9)
			{
				////mark[i+1] = 0;
				coeffX[i] = (i > j) ? (0) : 1;
				coeffX[j] = (j > i) ? (0) : 1;
			}
			else
			{
				///Nothing will be done 
			}
		}
	}

	for (int i = 0; i < 40; i++)//for Y
	{
		for (int j = 0; j < 40; j++)
		{
			if (i == 20 || j == 20)
				continue;
			if (i == j)
				continue;
#pragma warning(suppress: 6385)
			if (coeffY[i + 40] == 0 || coeffY[j + 40] == 0)
				continue;
			cosMy = cosD1y = cosD1y = 0;
			for (int t = 0; t < sizeofpoint; t++)
			{
				cosMy += coeffmatrix4Y[i + t * 40] * coeffmatrix4Y[j + t * 40];
				cosD1y += (coeffmatrix4Y[i + t * 40] * coeffmatrix4Y[i + t * 40]);
				cosD2y += (coeffmatrix4Y[j + t * 40] * coeffmatrix4Y[j + t * 40]);
			}
			cos_cory = cosMy / sqrt(cosD1y) / sqrt(cosD2y);
			if (fabs(cos_cory)  > maxrelationy) maxrelationy = fabs(cos_cory);
			if (cos_cory > 0.9)
			{
				////mark[i+1] = 0;
				coeffY[i] = (i > j) ? (0) : 1;
				coeffY[j] = (j > i) ? (0) : 1;
			}
			else
			{
				///Nothing will be done 
			}
		}
	}

	assert(maxrelationx < 1);
	assert(maxrelationy < 1);
	std::memcpy(coeffmatrix, coeffY, sizeof(int)* 40);
	std::memcpy(coeffmatrix + 40, coeffX, sizeof(int)* 40);

	cout << "X参数最大列方向性余弦为" << maxrelationx << endl;
	cout << "Y参数最大列方向性余弦为" << maxrelationy << endl;

	delete[] coeffmatrix4X;
	delete[] coeffmatrix4Y;
	return;
}



void RPCEvaluationHelper::RPCModel_AllInter_Simplified(int* mark)
{
	//First let's clear out
	Clear();
	for (int i = 0; i < 80; i++)
	{
		mark[i] = 1;
	}
	DoRPC3(true);

	////RPCModel_Selected(mark);  no,not here
	do{
		///Translate mark
		int marks[20], markl[20], marka[20], markb[20];
		TranslateMark(mark, marks, markl, marka, markb);
		int _parasize = 0; int _parasize1 = 0;
		for (int i = 0; i < 40; i++)
		{
			if (mark[i] == 1)
				_parasize++;
		}
		for (int i = 40; i < 80; i++)
		{
			if (mark[i] == 1)
				_parasize1++;
		}
		_parasize--;
		_parasize1--;
		double* a = new double[_parasize];
		double* aa = new double[_parasize*_parasize];
		double* al = new double[_parasize];
		double* b = new double[_parasize1]; double* bb = new double[_parasize1 * _parasize1]; double* bl = new double[_parasize1];
		double aext[20], NumL, DenL, NumS, DenS, P, L, H;
		double* x = new double[_parasize];  memset(x, 0, sizeof(double)*_parasize);
		double* x1 = new double[_parasize1]; memset(x1, 0, sizeof(double)*_parasize1);
		double sigema0 = 1e-10;
		double sigema = 1e-10;

		int calculatedtime = 0;
		do
		{
			if (!Switched())
			{
				ISWitch();
			}
			sigema0 = sigema;
			memset(a, 0, sizeof(double)* _parasize); memset(aa, 0, sizeof(double)* _parasize * _parasize); memset(al, 0, sizeof(double)* _parasize);
			memset(b, 0, sizeof(double)* _parasize1); memset(bb, 0, sizeof(double)* _parasize1 * _parasize1); memset(bl, 0, sizeof(double)* _parasize1);
			sigema = 0;

			for (int i = 0; i < PointList.size(); i++)
			{
				P = PointList[i].ground.lat;
				L = PointList[i].ground.lon;
				H = PointList[i].ground.h;
				RPCPLH2a(P, L, H, aext);
				NumL = satOrbitHelper.RPCPLH(_s, P, L, H);
				DenL = satOrbitHelper.RPCPLH(_l, P, L, H);
				NumS = satOrbitHelper.RPCPLH(_a, P, L, H);
				DenS = satOrbitHelper.RPCPLH(_b, P, L, H);

				double atemp[39] = { 0.0 };
				for (int j = 0; j < 20; j++)
				{
					atemp[j] = aext[j] * marks[j] / DenL;
				}

				for (int j = 20; j < 39; j++)
				{
					atemp[j] = -aext[j - 19] * markl[j - 19] * NumL / DenL / DenL;
				}
				int reader = 0;

				for (int j = 0; j < 39; j++)
				{
					if (atemp[j] != 0)
					{
						a[reader] = atemp[j];
						reader++;
					}
				}

				double chekcervalue = (PointList[i].Pixel.y - NumL / DenL);
				satOrbitHelper.pNormal(a, _parasize, chekcervalue, aa, al, 1.0);

				double btemp[39] = { 0.0 };
				for (int j = 0; j < 20; j++)
				{
					btemp[j] = aext[j] * marka[j] / DenS;
				}

				for (int j = 20; j < 39; j++)
				{
					btemp[j] = -aext[j - 19] * markb[j - 19] * NumS / DenS / DenS;
				}

				reader = 0;
				for (int j = 0; j < 39; j++)
				{
					if (btemp[j] != 0)
					{
						b[reader] = btemp[j];
						reader++;
					}
				}

				double valuechecker = PointList[i].Pixel.x - NumS / DenS;
				satOrbitHelper.pNormal(b, _parasize1, PointList[i].Pixel.x - NumS / DenS, bb, bl, 1.0);
				sigema += (PointList[i].Pixel.x - NumS / DenS)*(PointList[i].Pixel.x - NumS / DenS);
			}


			if (sigema < 1e-10) break;
			memset(x, 0, sizeof(double)* _parasize);

			satOrbitHelper.GaussExt(aa, al, x, _parasize);

			int realcount = 0;
			for (int j = 0; j < 20; j++)
			{
				_s[j] = (marks[j] == 0) ? 0 : _s[j] + x[realcount];
				realcount = (marks[j] == 0) ? realcount : realcount + 1;
			}

			//realcount = 0;
			for (int j = 20; j < 39; j++)
			{
				//_l[j - 19] += x[realcount];
				_l[j - 19] = (markl[j - 19] == 0) ? 0 : _l[j - 19] + x[realcount];
				realcount = (markl[j] == 0) ? realcount : realcount + 1;
			}

			memset(x1, 0, sizeof(double)* _parasize1);
			satOrbitHelper.GaussExt(bb, bl, x1, _parasize1);

			realcount = 0;
			for (int j = 0; j < 20; j++)
			{
				_a[j] = (marka[j] == 0) ? 0 : _a[j] + x1[realcount];
				realcount = (marka[j] == 0) ? realcount : realcount + 1;
			}

			//realcount = 0;
			for (int j = 20; j < 39; j++)
			{
				if (markb[j - 19] == 0) _b[j - 19] += markb[j - 19];
				else
				{
					_b[j - 19] = _b[j - 19] + x1[realcount];
					realcount++;
				}

			}

			_l[0] = 1.0;
			_b[0] = 1.0;

			calculatedtime++;

		} while (calculatedtime < 20);

		int xsizebefore, ysizebefore;
		int xsizeafter, ysizeafter;
		xsizebefore = ysizebefore = xsizeafter = ysizeafter = 0;
		if (!CheckSelection(mark, xsizebefore, ysizebefore))
			break;
		RPCModel_Selected(mark);
		if (!CheckSelection(mark, xsizeafter, ysizeafter))
			break;
		if (xsizebefore == xsizeafter && ysizebefore == ysizeafter)
			break;
	} while (true);  ///stop iteration when it is necessary
}

#pragma warning(suppress: 6262)
bool RPCEvaluationHelper::DoRPCAllInter()
{
	Clear();
	DoRPC3(true);

	for (int i = 0; i < 20; i++)
	{
		cout << _s[i] << endl;
		cout << _l[i] << endl;
		cout << _a[i] << endl;
		cout << _b[i] << endl;
	}

	
	double a[39], aa[1521], al[39];
	double bb[1521], bl[39];
	double aext[20], NumL, DenL, NumS, DenS, P, L, H;
	double x[39];

	double sigema0 = 1e-10;
	double sigema = 1e-10;


	int calculatedtime = 0;
	double sigmeold = 1000;
	do
	{
		if (!Switched())
		{
			ISWitch(); //Dont forget that!
		}

		sigema0 = sigema;
		memset(a, 0, sizeof(double)* 39);
		memset(aa, 0, sizeof(double)* 1521);
		memset(al, 0, sizeof(double)* 39);
		memset(bb, 0, sizeof(double)* 1521);
		memset(bl, 0, sizeof(double)* 39);

		sigema = 0;

		for (int i = 0; i < PointList.size(); i++)
		{
			P = PointList[i].ground.lat;
			L = PointList[i].ground.lon;
			H = PointList[i].ground.h;

			RPCPLH2a(P, L, H, aext);
			NumL = satOrbitHelper.RPCPLH(_s, P, L, H);
			DenL = satOrbitHelper.RPCPLH(_l, P, L, H);
			NumS = satOrbitHelper.RPCPLH(_a, P, L, H);
			DenS = satOrbitHelper.RPCPLH(_b, P, L, H);

			for (int j = 0; j < 20; j++)
			{
				a[j] = aext[j] / DenL;
			}

			for (int j = 20; j < 39; j++)
			{
				a[j] = -aext[j - 19] * NumL / DenL / DenL;
			}

			satOrbitHelper.pNormal(a, 39, (PointList[i].Pixel.y - NumL / DenL), aa, al, 1.0);

			sigema += (PointList[i].Pixel.y - NumL / DenL)*(PointList[i].Pixel.y - NumL / DenL);

			for (int j = 0; j < 20; j++)
			{
				a[j] = aext[j] / DenS;
			}

			for (int j = 20; j < 39; j++)
			{
				a[j] = -aext[j - 19] * NumS / DenS / DenS;
			}

			satOrbitHelper.pNormal(a, 39, PointList[i].Pixel.x - NumS / DenS, bb, bl, 1.0);
			sigema += (PointList[i].Pixel.x - NumS / DenS)*(PointList[i].Pixel.x - NumS / DenS);
		}

		sigema /= PointList.size();

		cout << sigema << endl;
		//break;
		//if (sigema < 0.005&&sigema > sigmeold) break;
		
		if (sigema < 1e-10) break;

		sigmeold = sigema;
		//cout<<sigema<<endl;

		memset(x, 0, sizeof(double)* 39);

		satOrbitHelper.GaussExt(aa, al, x, 39);

		for (int j = 0; j < 20; j++)
		{
			_s[j] += x[j];
		}

		for (int j = 20; j < 39; j++)
		{
			_l[j - 19] += x[j];
		}

		memset(x, 0, sizeof(double)* 39);
		satOrbitHelper.GaussExt(bb, bl, x, 39);


		for (int j = 0; j < 20; j++)
		{
			_a[j] += x[j];
		}

		for (int j = 20; j < 39; j++)
		{
			_b[j - 19] += x[j];
		}

		_l[0] = 1.0;
		_b[0] = 1.0;

		calculatedtime++;


		//Checker();
	} while (calculatedtime < 10);

	
	return true;
}

bool RPCEvaluationHelper::DoRPC1Same()
{
	Clear();

	double a[4], aa[16], al[4];
	std::memset(a, 0, sizeof(double)* 4);
	std::memset(aa, 0, sizeof(double)* 16);
	std::memset(al, 0, sizeof(double)* 4);

	double bb[16], bl[4];

	std::memset(bb, 0, sizeof(double)* 16);
	std::memset(bl, 0, sizeof(double)* 4);

	for (int i = 0; i < PointList.size(); i++)
	{
		a[0] = 1.0; a[1] = PointList[i].ground.lon; a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;

		satOrbitHelper.pNormal(a, 4, PointList[i].Pixel.y, aa, al, 1.0);

		a[0] = 1.0; a[1] = PointList[i].ground.lon; a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;

		satOrbitHelper.pNormal(a, 4, PointList[i].Pixel.x, bb, bl, 1.0);
	}

	satOrbitHelper.Gauss(aa, al, 4);
	satOrbitHelper.Gauss(bb, bl, 4);

	Clear();

	memcpy(_s, al, sizeof(double)* 4);
	memcpy(_a, bl, sizeof(double)* 4);

	_l[0] = 1.0;
	_b[0] = 1.0;

	return true;
}

bool RPCEvaluationHelper::DoRPC2Same(bool bIteration)
{
	if (!Switched())
		ISWitch();

	Clear();


	if (bIteration)
		DoRPC1Same();

	double a[10], aa[100], al[10];
	double bb[100], bl[10];

	std::memset(a, 0, sizeof(double)* 10);
	std::memset(aa, 0, sizeof(double)* 100);
	std::memset(al, 0, sizeof(double)* 10);
	std::memset(bb, 0, sizeof(double)* 100);
	std::memset(bl, 0, sizeof(double)* 10);

	for (int i = 0; i < PointList.size(); i++)
	{
		a[0] = 1.0;
		a[1] = PointList[i].ground.lon;
		a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;
		a[4] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[5] = PointList[i].ground.lon * PointList[i].ground.h;
		a[6] = PointList[i].ground.lat * PointList[i].ground.h;
		a[7] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[8] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[9] = PointList[i].ground.h * PointList[i].ground.h;

		satOrbitHelper.pNormal(a, 10, PointList[i].Pixel.y, aa, al, 1.0);

		a[0] = 1.0;
		a[1] = PointList[i].ground.lon;
		a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;
		a[4] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[5] = PointList[i].ground.lon * PointList[i].ground.h;
		a[6] = PointList[i].ground.lat * PointList[i].ground.h;
		a[7] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[8] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[9] = PointList[i].ground.h * PointList[i].ground.h;

		satOrbitHelper.pNormal(a, 10, PointList[i].Pixel.x, bb, bl, 1.0);
	}


	satOrbitHelper.Gauss(aa, al, 10);
	satOrbitHelper.Gauss(bb, bl, 10);

	Clear();

	memcpy(_s, al, sizeof(double)* 10);
	memcpy(_a, bl, sizeof(double)* 10);

	_l[0] = 1.0;
	_b[0] = 1.0;

	return true;
}

bool RPCEvaluationHelper::DoRPC3Same(bool bIteration)
{
	if (!Switched())
		ISWitch();

	Clear();

	if (bIteration)
		DoRPC2Same(true);

	double a[20], aa[400], al[20];
	double bb[400], bl[20];

	std::memset(a, 0, sizeof(double)* 20);
	std::memset(aa, 0, sizeof(double)* 400);
	std::memset(al, 0, sizeof(double)* 20);
	std::memset(bb, 0, sizeof(double)* 400);
	std::memset(bl, 0, sizeof(double)* 20);

	for (int i = 0; i < PointList.size(); i++)
	{
		a[0] = 1.0;
		a[1] = PointList[i].ground.lon;
		a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;
		a[4] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[5] = PointList[i].ground.lon * PointList[i].ground.h;
		a[6] = PointList[i].ground.lat * PointList[i].ground.h;
		a[7] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[8] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[9] = PointList[i].ground.h * PointList[i].ground.h;
		a[10] = PointList[i].ground.lat  * PointList[i].ground.lon * PointList[i].ground.h;
		a[11] = pow(PointList[i].ground.lon, 3);
		a[12] = PointList[i].ground.lon * pow(PointList[i].ground.lat, 2);
		a[13] = PointList[i].ground.lon * pow(PointList[i].ground.h, 2);
		a[14] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.lat;
		a[15] = pow(PointList[i].ground.lat, 3);
		a[16] = PointList[i].ground.lat * pow(PointList[i].ground.h, 2);
		a[17] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.h;
		a[18] = pow(PointList[i].ground.lat, 2) * PointList[i].ground.h;
		a[19] = pow(PointList[i].ground.h, 3);

		satOrbitHelper.pNormal(a, 20, PointList[i].Pixel.y, aa, al, 1.0);

		a[0] = 1.0;
		a[1] = PointList[i].ground.lon;
		a[2] = PointList[i].ground.lat;
		a[3] = PointList[i].ground.h;
		a[4] = PointList[i].ground.lon * PointList[i].ground.lat;
		a[5] = PointList[i].ground.lon * PointList[i].ground.h;
		a[6] = PointList[i].ground.lat * PointList[i].ground.h;
		a[7] = PointList[i].ground.lon * PointList[i].ground.lon;
		a[8] = PointList[i].ground.lat * PointList[i].ground.lat;
		a[9] = PointList[i].ground.h * PointList[i].ground.h;
		a[10] = PointList[i].ground.lat  * PointList[i].ground.lon * PointList[i].ground.h;
		a[11] = pow(PointList[i].ground.lon, 3);
		a[12] = PointList[i].ground.lon * pow(PointList[i].ground.lat, 2);
		a[13] = PointList[i].ground.lon * pow(PointList[i].ground.h, 2);
		a[14] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.lat;
		a[15] = pow(PointList[i].ground.lat, 3);
		a[16] = PointList[i].ground.lat * pow(PointList[i].ground.h, 2);
		a[17] = pow(PointList[i].ground.lon, 2) * PointList[i].ground.h;
		a[18] = pow(PointList[i].ground.lat, 2) * PointList[i].ground.h;
		a[19] = pow(PointList[i].ground.h, 3);

		satOrbitHelper.pNormal(a, 20, PointList[i].Pixel.x, bb, bl, 1.0);
	}

	double x[20];
	memcpy(x, _l, sizeof(double)* 20);
	satOrbitHelper.GaussExt(aa, al, x, 20);

	memcpy(al, x, sizeof(double)* 20);

	std::memset(_l, 0, sizeof(double)* 20);
	memcpy(_s, al, sizeof(double)* 20);

	memcpy(x, _a, sizeof(double)* 20);

	satOrbitHelper.GaussExt(bb, bl, x, 20);
	memcpy(bl, x, sizeof(double)* 20);
	memcpy(_a, bl, sizeof(double)* 20);

	_l[0] = 1.0;
	_b[0] = 1.0;

	return true;
}

void RPCEvaluationHelper::OutPutRpc()
{
	char reader[1024];
	sprintf_s(reader, "%s", _checkfile.c_str());

	FILE *fp;
	fopen_s(&fp, reader, "w");
	if (fp == NULL) return;

	fprintf(fp, "LINE_OFF: %+010.2lf pixels\n", np_array[8]);
	fprintf(fp, "SAMP_OFF: %+010.2lf pixels\n", np_array[6]);
	fprintf(fp, "LAT_OFF: %+011.8lf degrees\n", np_array[2]);
	fprintf(fp, "LONG_OFF: %+012.8lf degrees\n", np_array[0]);
	fprintf(fp, "HEIGHT_OFF: %+08.3lf meters\n", np_array[4]);

	fprintf(fp, "LINE_SCALE: %+010.2lf pixels\n", np_array[9]);
	fprintf(fp, "SAMP_SCALE: %+010.2lf pixels\n", np_array[7]);
	fprintf(fp, "LAT_SCALE: %0+12.8lf degrees\n", np_array[3]);
	fprintf(fp, "LONG_SCALE: %+012.8lf degrees\n", np_array[1]);
	fprintf(fp, "HEIGHT_SCALE: %+08.3lf meters\n", np_array[5]);


	int i;
	for (i = 0; i < 20; i++)
	{
		fprintf(fp, "LINE_NUM_COEFF_%d:  %+50.40le\n", i + 1, _s[i]);
	}

	for (i = 0; i < 20; i++)
	{
		fprintf(fp, "LINE_DEN_COEFF_%d:  %+50.40le\n", i + 1, _l[i]);
	}

	for (i = 0; i < 20; i++)
	{
		fprintf(fp, "SAMP_NUM_COEFF_%d:  %+50.40le\n", i + 1, _a[i]);
	}
	for (i = 0; i < 20; i++)
	{
		fprintf(fp, "SAMP_DEN_COEFF_%d:  %+50.40le\n", i + 1, _b[i]);
	}
	fclose(fp);
	return;
}

bool RPCEvaluationHelper::IConductor(RPCalculator _messager)
{
	_bInded = true;
	_modelinited = true;
	_currentcalculator = _messager;
	switch (_messager)
	{
	case Linear_POW1_DIFFER:
	{

							   cout << "Method: Linear, POW:1. Differ:true." << endl;
							   DoRPC1();
	}
		break;
	case Linear_POW2_DIFFER:
	{
							   cout << "Method:Linear, POW:2. Differ: true.Iteration:false" << endl;
							   DoRPC2(false);
	}
		break;
	case Linear_POW2_DIFFER_Iteration:
	{
										 cout << "Method:Linear, POW:2. Differ: true. Iteration:true" << endl;
										 DoRPC2(true);
	}
		break;
	case Linear_POW3_DIFFER:
	{
							   cout << "Method: Linear, POW:3. Differ:true" << endl;
							   DoRPC3(false);
	}
		break;
	case Linear_POW3_DIFFER_Iteration:
	{
										 cout << "Method:Linear, POW:3. Differr:true." << endl;
										 DoRPC3(true);
	}
		break;
	case Linear_POW1_Same:
	{
							 cout << "Method:Linear,POW:1.Differ:false.Iteration:false" << endl;
							 DoRPC1Same();
	}
		break;
	case Linear_POW2_Same:
	{
							 cout << "Method: Linear, POW:2. Differ:false,Iteration:false" << endl;
							 DoRPC2Same(false);
	}
		break;
	case Linear_POW2_Same_Iteration:
	{
									   cout << "Method: Linear. POW:2. Differ:false.Iteration:true" << endl;
									   DoRPC2Same(true);
	}
		break;
	case Linear_POW3_Same:
	{
							 cout << "Method: Linear. POW:3. Differ:false.Iteration:false" << endl;
							 DoRPC3Same(false);
	}
		break;
	case Linear_POW3_Same_Iteration:
	{
									   cout << "Method:Linear. POW:3. Differ:false. Iteration:true" << endl;
									   DoRPC3Same(true);
	}
		break;
	case Linear_POW1_TrueSame:
	{
								 cout << "Method:Linear,POW:1.Differ:true true.Iteration:false." << endl;
								 DoRPC1TrueSame();
	}
		break;
	case Linear_POW2_TrueSame:
	{
								 cout << "Method:Linear.POW2.Differ:true,true.Iteration:false" << endl;
								 DoRPC2TrueSame(false);
	}
		break;
	case Linear_POW2_TrueSame_Iteration:
	{
										   cout << "Method:Linear.POW2.Differ:true,true.Iteration:true" << endl;
										   DoRPC2TrueSame(true);
	}
		break;
	case Linear_POW3_TrueSame:
	{
								 cout << "Method:Linear.POW3.Differ:true,true.Iteration:false" << endl;
								 DoRPC3TrueSame(false);
	}
		break;
	case Linear_POW3_TrueSame_Iteration:
	{
										   cout << "Method:Linear.POW3.Differ:true,true.Iteration:true" << endl;
										   DoRPC3TrueSame(true);
	}
		break;
	case Linear_POW3_AllIteration:
	{
									 cout << "All interation" << endl;
									 DoRPCAllInter();
	}
		break;
	default:
		break;
	}

	return true;
}

///<Purpose>Read Check point</Purpose>
///<Paras>
/// <Para name="filepath">Path to Check points file</Para>
///<Paras>
bool RPCEvaluationHelper::GetCPoint(const char* chkfilepath)
{
	ifstream readproj;
	char readerproj[2048];
	readproj.open("projinfo.txt", ios::in);
	if (!readproj.is_open())
	{
		cout << "缺失投影信息文件，请将projinfo.txt文件放置到执行程序RPCModel.exe所在目录" << endl;
		cout << "若为调试模式，请将projinfo.txt放置到代码所在目录下" << endl;
		cout << "程序即将退出..." << endl;
		abort();
	}
	readproj.getline(readerproj, 2048);
	cout << "Infomation about reference map" << endl;
	cout << readerproj << endl;
	readproj.close();


	_checkfile = static_cast<string>(chkfilepath);

	//Open Checkpoint file
	ifstream reader;
	reader.open(chkfilepath, ios::in);
	CPointList.clear();

	bool _isLonLat = false;
	do
	{
		char readerstring[1024];
		reader.getline(readerstring, 1024);

		SPoint sp;
		//81               1015.218246         1843.687502       525588.389400      4825926.068400         1329.872029
		//对返回值进行判断
		int idd;
		sscanf_s(readerstring, "%d %lf %lf %lf %lf %lf", &idd
			, &sp.Pixel.x, &sp.Pixel.y, &sp.ground.X, &sp.ground.Y, &sp.ground.Z);

		sprintf_s(sp.id, "%d", idd);

		if (fabs(sp.ground.X) < 180 && fabs(sp.ground.Y) < 90)
		{
			_isLonLat = true;
			sp.ground.lon = sp.ground.X;
			sp.ground.lat = sp.ground.Y;
			sp.ground.X = -9999;
			sp.ground.Y = -9999;
		}
		else
			;//	imgin.ConvertToWGS84(sp.ground.X, sp.ground.Y, readerproj, sp.ground.lon, sp.ground.lat);
		sp.Pixel.y -= 0;

		//satOrbitHelper.MapXY2LatLon(sp.ground.X,sp.ground.Y,sp.ground.lat,sp.ground.lon);

		///Object obtained, print it( Printing is not necessary)
		cout << sp << endl;
		sp.ground.h = sp.ground.Z;

		//Store it into container,faily easy
		CPointList.push_back(sp);

	} while (reader.peek() != EOF); // if we have not reach the end, one more time 

	///Tell us how many points we got
	cout << "Check Points Got:" << CPointList.size() << endl;

	reader.close();

	if (_isLonLat)
	{
		cout << "检查点文件中出现经纬度，程序默认按WGS84经纬度处理，请核实！" << endl;
		cout << "当精度出现异常时，请留意此条信息" << endl;
	}

	else
	{
		cout << "检查点文件地面坐标为平面坐标，请确保exe所在路径下project.txt的投影信息与当前影像配套！" << endl;
		cout << "当程序精度出现异常时，请留意此条信息" << endl;
	}

	_CPointList.clear();
	for (int i = 0; i < PointList.size(); i++)
	{
		_CPointList.push_back(PointList[i]);
	}

	for (int i = 0; i < CPointList.size(); i++)
	{
		FPoint newpoint;
		newpoint.ground.X = CPointList[i].ground.X;
		newpoint.ground.Y = CPointList[i].ground.Y;
		newpoint.ground.Z = CPointList[i].ground.Z;
		newpoint.ground.lat = CPointList[i].ground.lat;
		newpoint.ground.lon = CPointList[i].ground.lon;
		newpoint.ground.h = CPointList[i].ground.h;
		memcpy(newpoint.id, CPointList[i].id, 256);
		newpoint.Pixel.x = CPointList[i].Pixel.x;
		newpoint.Pixel.y = CPointList[i].Pixel.y;
		_CPointList.push_back(newpoint);
	}
	return true;
}

void RPCEvaluationHelper::Compute_npara()
{
	///normalization
	double amax[5], amin[5];
	std::memset(np_array, 0, sizeof(double)* 10);
	amax[0] = PointList[0].ground.lon; amin[0] = PointList[0].ground.lon;
	amax[1] = PointList[0].ground.lat; amin[1] = PointList[0].ground.lat;
	amax[2] = PointList[0].ground.h; amin[2] = PointList[0].ground.h;
	amax[3] = PointList[0].Pixel.x; amin[3] = PointList[0].Pixel.x;
	amax[4] = PointList[0].Pixel.y; amin[4] = PointList[0].Pixel.y;

	for (int i = 0; i < PointList.size(); i++)
	{
		np_array[0] += PointList[i].ground.lon;
		np_array[2] += PointList[i].ground.lat;
		np_array[4] += PointList[i].ground.h;
		np_array[6] += PointList[i].Pixel.x;
		np_array[8] += PointList[i].Pixel.y;
		amax[0] = (PointList[i].ground.lon>amax[0]) ? PointList[i].ground.lon : amax[0];
		amin[0] = (PointList[i].ground.lon <= amin[0]) ? PointList[i].ground.lon : amin[0];
		amax[1] = (PointList[i].ground.lat > amax[1]) ? PointList[i].ground.lat : amax[1];
		amin[1] = (PointList[i].ground.lat <= amin[1]) ? PointList[i].ground.lat : amin[1];
		amax[2] = (PointList[i].ground.h > amax[2]) ? PointList[i].ground.h : amax[2];
		amin[2] = (PointList[i].ground.h <= amin[2]) ? PointList[i].ground.h : amin[2];
		amax[3] = (PointList[i].Pixel.x > amax[3]) ? PointList[i].Pixel.x : amax[3];
		amin[3] = (PointList[i].Pixel.x <= amin[3]) ? PointList[i].Pixel.x : amin[3];
		amax[4] = (PointList[i].Pixel.y > amax[4]) ? PointList[i].Pixel.y : amax[4];
		amin[4] = (PointList[i].Pixel.y <= amin[4]) ? PointList[i].Pixel.y : amin[4];
	}

	for (int i = 0; i < CPointList.size(); i++)
	{
		np_array[0] += CPointList[i].ground.lon;
		np_array[2] += CPointList[i].ground.lat;
		np_array[4] += CPointList[i].ground.h;
		np_array[6] += CPointList[i].Pixel.x;
		np_array[8] += CPointList[i].Pixel.y;
		amax[0] = (CPointList[i].ground.lon>amax[0]) ? PointList[i].ground.lon : amax[0];
		amin[0] = (CPointList[i].ground.lon <= amin[0]) ? CPointList[i].ground.lon : amin[0];
		amax[1] = (CPointList[i].ground.lat > amax[1]) ? CPointList[i].ground.lat : amax[1];
		amin[1] = (CPointList[i].ground.lat <= amin[1]) ? CPointList[i].ground.lat : amin[1];
		amax[2] = (CPointList[i].ground.h > amax[2]) ? CPointList[i].ground.h : amax[2];
		amin[2] = (CPointList[i].ground.h <= amin[2]) ? CPointList[i].ground.h : amin[2];
		amax[3] = (CPointList[i].Pixel.x > amax[3]) ? CPointList[i].Pixel.x : amax[3];
		amin[3] = (CPointList[i].Pixel.x <= amin[3]) ? CPointList[i].Pixel.x : amin[3];
		amax[4] = (CPointList[i].Pixel.y > amax[4]) ? CPointList[i].Pixel.y : amax[4];
		amin[4] = (CPointList[i].Pixel.y <= amin[4]) ? CPointList[i].Pixel.y : amin[4];
	}

	for (int i = 0; i < 5; i++)
	{
		np_array[i * 2] /= static_cast<double>(PointList.size() + CPointList.size()); ///这个是中心值
		np_array[i * 2 + 1] = max(fabs(amax[i] - np_array[i * 2]), fabs(amin[i] - np_array[i * 2])); ///这个是范围
		///对np_array[i * 2 + 1]加个判断，如果它趋近于0，将其设置以为适当的非0正值
		///对结果不会有影响，但是可以避免做除法时除0（导致中心化之后参数变成极大值）
		///add your code below here
		if (fabs(np_array[i * 2 + 1]) < 0.1)
			np_array[i * 2 + 1] = 1.0;
		///modification done

	}

}