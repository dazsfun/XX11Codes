#include "stdafx.h"
#include "MetaClass.h"

namespace MetaClass
{
	string Register[150] = { "Undistributed" };
}
bool CheckInput(string inputer, string& thevalue)
{
	cout << "�����롪��" << inputer.c_str() << endl;

	char reader[1024];
	cin >> reader;
	if (safestring::compare(reader, "NEXTMETAINFO", 12))
	{
		thevalue = "StopNode";
		return false;
	}

	thevalue = reader;
	return true;
}

vector<string> RegisterAttitude(string AttitudeRegister[])
{
	vector<string> resulting;
	int ncounter = 0;
	do
	{
		string heivalue;
		if (CheckInput(AttitudeRegister[ncounter], heivalue))
		{
			resulting.push_back(heivalue);
			ncounter++;
		}
		else
		{
			for (unsigned int i = 0; i < (15 - ncounter); i++)
			{
				resulting.push_back("NULLCONTENT");
			}
		}
	} while (ncounter < 15);

	return resulting;
}

///�������ļ�������Ϊ�������
//��ȡ˳��:
/// 
bool MetaClass::MetaClass(string _configFile)
{
	string _content = safestring::ReadFile(_configFile);

	Json::Reader reader;
	Json::Value root;

	Json::ValueType typehei = root.type();
	cout << _content.c_str()<<endl;
	bool registernotexisted = true;
	if ((registernotexisted = reader.parse(_content, root)))
	{
		int sizeofroot = root.size();
		int counter = 0;
		for (Json::ValueIterator itr = root.begin(); itr != root.end(); itr++) {
			// Print depth.
			string hei = itr.key().asString();
			string heinew = root.get(hei.c_str(), 0).asString();
			Register[counter] = heinew.c_str();
			counter++;

		}
	}

	cout << registernotexisted << endl;
	registernotexisted = !registernotexisted;
	bool RegisterNew = false;
	std::cout << "��Ҫע���µ�Ԫ��?(1\0)" << endl;
	cin >> RegisterNew;
	if (registernotexisted || RegisterNew)
	{
		string LineTimeRegister[15] = { "HTIME1"
			, "HTIME2AGG", "HTIME3SEG",
			"HTIME4LINE", "HTIME5ActLine", "HTIME6UTCTIME"
			, "HTIME7DATETIME", "HTIME8LINESEG", "HTIMELTIMEReservedPara1",
			"HTIMELTIMEReservedPara2", "HTIMELTIMEGPSReservedPara3", "HTIMELTIMEReservedPara4", "HTIMELTIMEReservedPara5",
			"HTIMELTIMEReservedPara6", "HTIMELTIMEReservedPara7" };
		string CalibRegister[15] = { "IAScan2Cameara0Roll"
			, "IAScan2Cameara1Pitch", "IAScan2Cameara2YAW",
			"IBCamera2Install0Roll", "IBCamera2Install1Pitch", "IBCamera2Install2YAW"
			, "ICInstall2Body0Roll", "ICInstall2Body1PITCH", "ICInstall2Body2YAW",
			"IDCalibratedPara0Roll", "IDCalibratedPara1PITCH", "IDCalibratedPara2YAW", "IEDateTimeString",
			"IEControlPointErrorX", "IEControlPointErrorY" };
		string AttitudeRegister[15] = { "Attitude01UTC"
			, "Attitude02DateTime", "Attitude03Type", "Attitude04ROLL",
			"Attitude05PITCH", "Attitude06YAW", "Attitude07Q0"
			, "Attitude08Q1", "Attitude09Q2", "Attitude10Q3",
			"Attitude11RollVec", "Attitude12PitchVec", "Attitude13YAWVec",
			"AttitudeReservedPara1", "AttitudeReservedPara2" };
		string gpsRegister[15] = { "GPS1UTC"
			, "GPS2DateTime", "GPS3X",
			"GPS4Y", "GPS5Z", "GPS6XV"
			, "GPS7YV", "GPS8ZV", "GPSReservedPara1",
			"GPSReservedPara2", "GPSReservedPara3", "GPSReservedPara4", "GPSReservedPara5",
			"GPSReservedPara6", "GPSReservedPara7" };
		vector<string> LineTimeDefi, attDef, gpsDef, CalibDefi;
		///ԭע���ļ������ڻ����Ѿ��𻵣���������ע��
		bool addnew = false;
		string tempread = Register[0];
		int hei = tempread.find("UnDistributed");
		if (tempread.find("UnDistributed") == string::npos)
		{
			std::cout << "��̬Ԫ���ݹ��ע����Ϣ�Ѵ��ڣ��Ƿ�����ע��?(1\0)" << endl;
			cin >> addnew;
		}
		else addnew = true;
		if (addnew) ///����ע����̬
		{
			cout << "�������ʾ������̬Ԫ�����ļ�ע�ᣬ" << endl;
			cout << "ע�������ʶ:NEXTMETAINFO" << endl;
			//ע����̬

			attDef = RegisterAttitude(AttitudeRegister);

			for (int i = 0; i < 15; i++)
			{
				Register[i + 0] = attDef[i].c_str();
			}
		}

		bool addOrbitNew = false;
		tempread = Register[15];
		if (tempread.find("Undistributed") == string::npos)
		{
			cout << "����Ѿ�ע�ᣬ��Ҫ����ע�᣿��1\0��" << endl;
			cin >> addOrbitNew;
		}
		else addOrbitNew = true;
		if (addOrbitNew)
		{
			//ע����
			//ע����̬

			gpsDef = RegisterAttitude(gpsRegister);

			for (int i = 0; i < 15; i++)
			{
				Register[i + 15] = gpsDef[i].c_str();
			}
		}

		bool addTimeNew = false;
		tempread = Register[30];
		if (tempread.find("UnDistributed") == string::npos)
		{
			cout << "��ʱ��Ϣ��ע�ᣬ��Ҫ����ע��(1\0)" << endl;
			cin >> addTimeNew;
		}
		else addTimeNew = true;
		if (addTimeNew)
		{
			//ע����ʱ
			//ע����̬


			LineTimeDefi = RegisterAttitude(LineTimeRegister);

			for (int i = 0; i < 15; i++)
			{
				Register[i + 30] = LineTimeDefi[i].c_str();
			}
		}

		//ע����������
		bool addCalibNew = false;
		tempread = Register[45];
		if (tempread.find("UnDistributed") == string::npos)
		{
			cout << "����������ע�ᣬ����ע�᣿(1\0)" << endl;
			cin >> addCalibNew;
		}
		else addCalibNew = true;
		if (addCalibNew)
		{


			CalibDefi = RegisterAttitude(CalibRegister);
			for (int i = 0; i < 15; i++)
			{
				Register[i + 45] = CalibDefi[i];
			}
		}

		for (int i = 60; i < 150; i++)
		{
			Register[i] = "UnDistributed";
		}
		Json::Value new_item;
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 15; j++)
			{
				if (i == 1)
				{
					if (attDef.size() < 15)
					{
						new_item[AttitudeRegister[j]] = Register[j];
					}
					else
						new_item[AttitudeRegister[j]] = attDef[j];
				}
				else if (i == 2)
				{
					if (gpsDef.size() < 15)
					{
						new_item[gpsRegister[j]] = Register[15 + j];
					}
					else
						new_item[gpsRegister[j]] = gpsDef[j];
				}
				else if (i == 0)
				{
					if (LineTimeDefi.size() < 15)
					{
						new_item[LineTimeRegister[j]] = Register[30 + j];
					}
					else
						new_item[LineTimeRegister[j]] = LineTimeDefi[j];
				}
				else if (i == 3)
				{
					if (CalibDefi.size() < 15)
					{
						new_item[CalibRegister[j]] = Register[45 + j];
					}
					else
						new_item[CalibRegister[j]] = CalibDefi[j];
				}
				else
				{
					char tempcharsprintf[1024];
					sprintf_s(tempcharsprintf, "%d", i * 15 + j);
					string tempallocate = "UnDistributed" + string(tempcharsprintf);
					new_item[tempallocate] = Register[i * 15 + j];
				}
			}
		}

		Json::FastWriter writer;
		string out = writer.write(new_item);

		safestring::WriteFile(_configFile, out);

		cout << out.c_str() << endl;

		return false;

	}
}
