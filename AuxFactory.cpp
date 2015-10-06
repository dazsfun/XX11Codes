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
		///配置文件不存在，开始自助定义过程
		string definitionstrings[11]
			= { "singlecount", "mainIDPosition", "mainIDLength",
			"satTimePosition", "satTimeLength", "satAttPosition", "satAttLength",
			"satGPSPosition", "satGPSLength", "satTimeCheckPosition",
			"satTimeCheckLength" };
		string translation[11] = { "字节数", "帧计数位置", "帧计数字节数", "时间数据起始位置",
			"时间数据占用字节数", "姿态数据起始位置", "姿态数据占用字节数",
			"轨道数据起始位置", "轨道数据占用字节数", "行时数据起始位置", "行时数据占用字节数" };
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

	///从配置文件中构建
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

///这个接口已经废弃，标注时间：2015、03、08 标注期：3个月
bool AuxFactory::IGiveMeFive(char* _auxcontent)
{
	///解析MainID编号
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
	///获取MainID
	unsigned char* temparray = new unsigned char[_information.size()];
	memset(temparray, 0, _information.size());
	for (int i = 0; i < _information.size(); i++)
	{
		temparray[i] = _information[i];
	}

	///获取列号
	curcolnumber = safestring::FormatInt(temparray + 98, 2);
	cout << "当前列号" << endl;
	cout << curcolnumber << endl;
	IAddLog("当前列号:%d", curcolnumber);

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

