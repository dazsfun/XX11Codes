#include "stdafx.h"
#include "AuxFactory.h"


AuxFactory::AuxFactory(string _configfilepath) :currentTime("NotAllocated"),
timeBackup(currentTime),
curOribt("NotAllocated"),
orbackup(curOribt),
curAtt("NotAllocated"),
attBackup(curAtt)
{
	SetPropertys(); curcolnumber = 0;
	string filecontent = safestring::ReadFile(_configfilepath);

	ifstream ifs(_configfilepath.c_str());
	if (filecontent == "Nothing" && !ifs.is_open())
	{
		///�����ļ������ڣ���ʼ�����������
		string definitionstrings[11]
			= { "singlecount", "mainIDPosition", "mainIDLength",
			"satTimePosition", "satTimeLength", "satAttPosition", "satAttLength",
			"satGPSPosition", "satGPSLength", "satTimeCheckPosition",
			"satTimeCheckLength" };
		string translation[11] = { "�ֽ���", "֡����λ��", "֡�����ֽ���", "ʱ��������ʼλ��",
			"ʱ������ռ���ֽ���", "��̬������ʼλ��", "��̬����ռ���ֽ���",
			"���������ʼλ��", "�������ռ���ֽ���", "��ʱ������ʼλ��", "��ʱ����ռ���ֽ���" };
		for (int i = 0; i < 11; i++)
		{
			cout << definitionstrings[i].c_str() << "( " << translation[i].c_str() << " )" << endl;
			int valuereader = 0;
			cin >> valuereader;
			void* pAddr = m_listPropertyAddr[i];
			(*(int*)pAddr) = valuereader;
		}

		//Serialize(_configfilepath);
		ofstream ofs(_configfilepath.c_str());
		ofs << Serialize().c_str();
		ofs.close();
		return;
	}

	///�������ļ��й���
	DeSerialize(filecontent.c_str());


}

AuxFactory::AuxFactory(const AuxFactory& _newfactory) :curAtt(_newfactory.curAtt),
attBackup(_newfactory.attBackup), curOribt(_newfactory.curOribt),
orbackup(_newfactory.orbackup),
currentTime(_newfactory.currentTime),
timeBackup(_newfactory.timeBackup)
{
	singlecount = _newfactory.singlecount;
	mainIDPosition = _newfactory.mainIDPosition;
	mainIDLength = _newfactory.mainIDLength;
	satTimeLength = _newfactory.satTimeLength;
	satTimePosition = _newfactory.satTimePosition;
	satAttLength = _newfactory.satAttLength;
	satAttPosition = _newfactory.satAttPosition;
	satGPSLength = _newfactory.satGPSLength;
	satGPSPosition = _newfactory.satAttPosition;
	satTimeCheckLength = _newfactory.satTimeCheckLength;
	satTimeCheckPosition = _newfactory.satTimeCheckPosition;
	curcolnumber = _newfactory.curcolnumber;
	framenum = _newfactory.framenum;

}

int AuxFactory::framecount;

AuxFactory::~AuxFactory()
{
}

void AuxFactory::SetPropertys()
{
	m_listName.clear();
	m_listPropertyAddr.clear();
	m_listType.clear();
	SetProperty("singlecount", asInt, &singlecount);
	SetProperty("mainIDPosition", asInt, &mainIDPosition);
	SetProperty("mainIDLength", asInt, &mainIDLength);
	SetProperty("satTimePosition", asInt, &satTimePosition);
	SetProperty("satTimeLength", asInt, &satTimeLength);
	SetProperty("satAttPosition", asInt, &satAttPosition);
	SetProperty("satAttLength", asInt, &satAttLength);
	SetProperty("satGPSPosition", asInt, &satGPSPosition);
	SetProperty("satGPSLength", asInt, &satGPSLength);
	SetProperty("satTimeCheckPosition", asInt, &satTimeCheckPosition);
	SetProperty("satTimeCheckLength", asInt, &satTimeCheckLength);
}

///����ӿ��Ѿ���������עʱ�䣺2015��03��08 ��ע�ڣ�3����
bool AuxFactory::IGiveMeFive(char* _auxcontent)
{
	///����MainID���
	double rate = (static_cast<double>(singlecount) / 8.0);
	int mainIDbytecounts = (static_cast<double>(singlecount) / 8.0) * mainIDLength;
	unsigned char* temppointer = new unsigned char[4];
	memset(temppointer, 0, sizeof(char)* 4);
	memcpy(temppointer + (4 - mainIDbytecounts), _auxcontent + mainIDPosition, sizeof(char)* mainIDLength);
	int heitemp = 0;
	memcpy(&heitemp, temppointer, sizeof(char)* 4);

	cout << heitemp << endl;

	vector<unsigned char> vtemp;
	safestring::TrimBytes(_auxcontent + static_cast<int>(static_cast<double>(singlecount) / 8.0 * satTimePosition),
		vtemp, static_cast<int>(static_cast<double>(singlecount) / 8.0 * satTimeLength), 12);

	char* timeinitial = new char[static_cast<int>(static_cast<double>(singlecount) / 8.0 * satTimeLength)];
	//LineTimeAPI ltimeapi(timeinitial);


	return true;
}

bool AuxFactory::IEquipment(vector<unsigned char> _information)
{
	///��ȡMainID
	unsigned char* temparray = new unsigned char[_information.size()];
	memset(temparray, 0, _information.size());
	for (int i = 0; i < _information.size(); i++)
	{
		temparray[i] = _information[i];
	}

	///��ȡ�к�
	curcolnumber = safestring::FormatInt(temparray + 98, 2);
	cout << "��ǰ�к�" << endl;
	cout << curcolnumber << endl;
	IAddLog("��ǰ�к�:%d", curcolnumber);

	int mainIdIS = safestring::FormatInt(temparray + mainIDPosition, mainIDLength);
	cout << mainIdIS << endl;

	LineTimeAPI linetimeAPI(temparray + satTimePosition);
	currentTime = linetimeAPI;
	currentTime.ACTLine(curcolnumber);
	OrbitAPI orbitAPI(temparray + satGPSPosition);
	curOribt = orbitAPI;
	AttitudeAPI attitudeAPI(temparray + satAttPosition);
	curAtt = attitudeAPI;
	delete[] temparray;
	return true;
}

int AuxFactory::GetCol()
{
	return curcolnumber;
}

LineTimeAPI& AuxFactory::LineTime()
{
	return currentTime;
}

AttitudeAPI& AuxFactory::Attitude()
{
	return curAtt;
}

OrbitAPI& AuxFactory::Orbit()
{
	return curOribt;
}

