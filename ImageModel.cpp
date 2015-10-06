#include "stdafx.h"
#include "ImageModel.h"
#include "EquinoxBasedT2C.h"

ImageModel::ImageModel(OrbitAPI* orbApi, AttitudeAPI* attApi, LineTimeAPI* ltimeAPI) :currentOrbit(OrbitAPI("currentOrbit")), currentAttitude(AttitudeAPI("currentAttitude"))
{
	attType = 0;
	int componentCounter = 0;
	do{
		if (orbApi->IGetTime() <= 0)
		{
			cout << "轨道数据加载完毕" << endl;
			LogAPI logapi;
			logapi.IAddLog("轨道数据加载完毕");
			break;
		}

		_orbitComponents.push_back(*orbApi);
		orbApi++;
	} while (true);

	componentCounter = 0;
	do
	{
		if (attApi->IGetTime() <= 0)
		{
			cout << "姿态数据加载完毕" << endl;
			LogAPI logapi;
			logapi.IAddLog("姿态数据加载完毕");                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
			break;
		}

		_attComponents.push_back(*attApi);
		attApi++;
	} while (true);

	componentCounter = 0;
	do
	{
		if (ltimeAPI->IGetTime() <= 0)
		{
			cout << "行时数据加载完毕" << endl;
			LogAPI logapi;
			logapi.IAddLog("行时数据加载完毕");
			break;
		}

		_timeComponents.push_back(*ltimeAPI);
		ltimeAPI++;
	} while (true);

	LogAPI logapi;
	cout << "辅助数据加载完毕" << endl;
	logapi.IAddLog("辅助数据加载完毕...");
}


ImageModel::~ImageModel()
{
	
}

bool ImageModel::SelfCheck()
{
	bool returnstatus = true;
	if (attType != 0)
	{
		IAddLog("目前只支持本体到轨道的姿态角的处理");
		returnstatus = false;
	}

	//所有辅助数据里目前机制下只有轨道数据会重复 
	//只对重复的情况做优化
	safestring::ClearVectorSame<OrbitAPI>(_orbitComponents);
	safestring::ClearVectorSame<AttitudeAPI>(_attComponents);
	safestring::ClearVectorSame<LineTimeAPI>(_timeComponents);

	returnstatus = safestring::CheckVectorSequence<OrbitAPI>(_orbitComponents,"轨道");
	returnstatus = safestring::CheckVectorSequence<AttitudeAPI>(_attComponents,"姿态");
	returnstatus = safestring::CheckVectorSequence<LineTimeAPI>(_timeComponents,"行时");

	return returnstatus;
}

void ImageModel::IBody2WGS84(string _configfilename, double* body2WGS84, double* positionXYZ)
{
	ifstream ifs;
	
	if (_configfilename.rfind("\\") != _configfilename.size() - 1)
	{
		_configfilename += "\\";
	}
	string attApiDescription = _configfilename + "姿态API.json";
	string orbApiDescription = _configfilename + "轨道API.json";
	//ifstream ifs2;
	ifs.open(orbApiDescription.c_str(), ios::in);
	if (!ifs.is_open())
	{
		cout << "成像测试配置文件不存在，请检查" << endl;
		IAddLog("成像测试配置文件，请检查，指定的配置文件路径为%s",_configfilename.c_str());
		return;
	}

	char reader[1024];
	ifs.getline(reader, 1024);

	OrbitAPI orbitest(static_cast<string>(reader), asJSon);
	double j20002wgs84[9] = { 0 };
	double j20002wgs84dot[9] = { 0 };
	UTCAPI timeapi = orbitest.utcapi;
	Safeauxilary::GetJ2000ToWGS84(timeapi.year, timeapi.month, timeapi.day, timeapi.hour, timeapi.minute, timeapi.second, j20002wgs84, j20002wgs84dot);

	ifstream ifs2;
	ifs2.open(attApiDescription.c_str(), ios::in);
	ifs2.getline(reader, 1024);
	AttitudeAPI attest(static_cast<string>(reader), asJSon);
	double attRot[9] = { 0.0 };
	double q1, q2, q3, q4;
	CSatOrbit orbithelpernow;
	orbithelpernow.quat2matrix(attest.AttitudeQ1, attest.AttitudeQ2, attest.AttitudeQ3, attest.AttitudeQ0, attRot);
	orbithelpernow.invers_matrix(attRot, 3);
	orbithelpernow.matrix2quat(attRot, q1, q2, q3, q4);
	positionXYZ[0] = orbitest.GPSX;
	positionXYZ[1] = orbitest.GPSY;
	positionXYZ[2] = orbitest.GPSZ;

	double position[6] = { orbitest.GPSX, orbitest.GPSY, orbitest.GPSZ,
		orbitest.GPSXV, orbitest.GPSYV, orbitest.GPSZV };

	m_base.ConvertOrbit2WGS84(attRot,j20002wgs84, position, body2WGS84,true);

	//m_base.invers_matrix(body2WGS84, 3);
	ifs.close(); ifs2.close();
}



void ImageModel::IBody2WGS84Test(double time, double* body2WGS84, double* positionXYZ)
{
	//OrbitAPI currentOrbit("hei");
	m_base.LagrangianInterpolation(_orbitComponents, time, _orbitComponents.size(), currentOrbit);

	m_base.QuatInterpolation(_attComponents, time, _attComponents.size(), currentAttitude);
	double attRot[9];

	////本体到J2000
	helper.quat2matrix(currentAttitude.AttitudeQ1, currentAttitude.AttitudeQ2, currentAttitude.AttitudeQ3,
		currentAttitude.AttitudeQ0, attRot);

	UTCAPI timeapi(time - 3600 * 8 ,Enum_TimeStandard_ZY03);
	double j20002wgs84[9] = {0};
	double j20002wgs84dot[9] = {0};
	Safeauxilary::GetJ2000ToWGS84(timeapi.year, timeapi.month, timeapi.day, timeapi.hour, timeapi.minute, timeapi.second, j20002wgs84, j20002wgs84dot);
	//m_base.invers_matrix(attRot,3);
	positionXYZ[0] = currentOrbit.GPSX;
	positionXYZ[1] = currentOrbit.GPSY;
	positionXYZ[2] = currentOrbit.GPSZ;

	//string path = "d:\\" +  safestring::WholeTimeString() + "\\imagemodel";
	//currentAttitude.Serialize(path);
	//currentOrbit.Serialize(path);

	double position[6] = { currentOrbit.GPSX, currentOrbit.GPSY,
		currentOrbit.GPSZ, currentOrbit.GPSXV, currentOrbit.GPSYV,
		currentOrbit.GPSZV };

	//m_base.invers_matrix(attRot, 3);
	m_base.ConvertOrbit2WGS84(attRot, j20002wgs84, position, body2WGS84,true);
}

void ImageModel::IBody2WGS84(double time, double* body2WGS84,double* positionXYZ)
{
	if (attType == 0)
	{
		OrbitAPI currentOrbit("currentPoint");
		m_base.LagrangianInterpolation(_orbitComponents, time, _orbitComponents.size(), currentOrbit);
		AttitudeAPI currentAttitude("currentAttitude");
		m_base.QuatInterpolation(_attComponents, time, _attComponents.size(), currentAttitude);
		double attRot[9];
		UTCAPI timeapi(time, Enum_TimeStandard_ZY03);
		double j20002wgs84[9] = { 0 };
		double j20002wgs84dot[9] = { 0 };
		Safeauxilary::GetJ2000ToWGS84(timeapi.year, timeapi.month, timeapi.day, timeapi.hour, timeapi.minute, timeapi.second, j20002wgs84, j20002wgs84dot);
		m_base.Quat2Matrix(currentAttitude.AttitudeQ1, currentAttitude.AttitudeQ2, currentAttitude.AttitudeQ3,
			currentAttitude.AttitudeQ0, attRot);
		positionXYZ[0] = currentOrbit.GPSX;
		positionXYZ[1] = currentOrbit.GPSY;
		positionXYZ[2] = currentOrbit.GPSZ;
		double position[6] = { currentOrbit.GPSX, currentOrbit.GPSY,
			currentOrbit.GPSZ, currentOrbit.GPSXV, currentOrbit.GPSYV,
			currentOrbit.GPSZV};
		m_base.ConvertOrbit2WGS84(attRot,j20002wgs84, position, body2WGS84,true);
	}
}
