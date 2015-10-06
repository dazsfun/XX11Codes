#pragma once
#include <vector>
#include <string>
using namespace std;
#include "jsonSerializer.h"
#include "log.h"
#include "AttitudeAPI.h"
#include "OrbitAPI.h"
#include "LineTimeAPI.h"
#include "safestring.h"

///����ඨ���˸������ݵĽ���ģ��
///��Ҫ�Ĳ���˵����ע��
///ע�������ݲ���Զ������ȷ�������������쳣����
///    ֻ��֤���ϸ����սӿ����õ�����»����ȷ��������
///ԭ�� 1.����ǻ����࣬����ͼ����һ�����ڿ�������֮ǰ��Ȼ�ñ�֤����������Ч��
///      2. safestring�������У��ṩ��������������Ч�Լ�鹦�ܺ��ַ���ת������
class AuxFactory : public JSONOBJECTSerializer
{
public:
	AuxFactory(string _configpath);
	AuxFactory(int _singlecount, int _mainIDPosition,
		int _mainIDLength, int _satTimePosition, int _satTimeLength,
		int _satAttPosition, int _satAttLength,
		int _satGPSPosition, int _satGPSLength,
		int _satTimeCheckPosition, int _satTimeCheckLength)
		: singlecount(_singlecount), mainIDPosition(_mainIDPosition), mainIDLength(_mainIDLength),
		satTimePosition(_satTimePosition), satTimeLength(_satTimeLength), satAttPosition(_satAttPosition),
		satAttLength(_satAttLength), satGPSPosition(_satAttPosition), satGPSLength(_satGPSLength),
		satTimeCheckPosition(_satTimeCheckPosition), satTimeCheckLength(_satTimeCheckLength),
		currentTime("{\"TIMEAGG\":0,\"TIMESEG\":0,\"LINE\":0,\"ActLine\":0,\"UTCTIME\":0,\"DATETIME\":\"NOT ALLOCATED\",\"LINESEG\":0,\"ReservedPara1\":\"ReservedPara\",\"ReservedPara2\":\"ReservedPara\",\"ReservedPara3\":\"ReservedPara\",\"ReservedPara4\":\"ReservedPara\",\"ReservedPara5\":\"ReservedPara\",\"ReservedPara6\":\"ReservedPara\",\"ReservedPara7\":\"ReservedPara\",\"ReservedPara0\":\"ReservedPara\"}", asJSon),
		timeBackup(currentTime),
		curOribt("{\"GPSUTC\":0,\"GPSDateTime\":\"NotAllocated\",\"GPSX\":0,\"GPSY\":0,\"GPSZ\":0,\"GPSXV\":0,\"GPSYV\":0,\"GPSZV\":0,\"ReservedPara1\":\"ReservedPara\",\"ReservedPara2\":\"ReservedPara\",\"ReservedPara3\":\"ReservedPara\",\"ReservedPara4\":\"ReservedPara\",\"ReservedPara5\":\"ReservedPara\",\"ReservedPara6\":\"ReservedPara\",\"ReservedPara7\":\"ReservedPara\"}", asJSon),
		orbackup(curOribt),
		curAtt("{\"AttitudeUTC\":0,\"AttitudeDateTime\":\"NotAllocated\",\"AttitudeType\":0,\"AttitudeROLL\":0,\"AttitudePITCH\":0,\"AttitudeYAW\":0,\"AttitudeQ0\":0,\"AttitudeQ1\":0,\"AttitudeQ2\":0,\"AttitudeQ3\":0,\"AttitudeROLLVec\":0,\"AttitudePITCHVec\":0,\"AttitudeYAWVect\":0,\"ReservedPara1\":\"ReservedPara\",\"ReservedPara2\":\"ReservedPara\"}", asJSon),
		attBackup(curAtt)
	{
		framecount = 0;
		SetPropertys();
		curcolnumber = 0;
	}

	AuxFactory(const AuxFactory& _newfactory);
	~AuxFactory();

private:
	string timeeApiDesriptor;
	string gpsAPIDescriptor;
	string attAPIDescriptor;
	int mainFrameNumber;
	

private:
	LineTimeAPI currentTime, timeBackup;
	OrbitAPI curOribt, orbackup;
	AttitudeAPI curAtt, attBackup;
	int curcolnumber;
	int framenum;
public:
	int singlecount;   ///һ����׼�ֳ���Ӧ��λ���������XX11IRS��12bit�����ģ�����12
	int mainIDPosition;
	int mainIDLength;
	int satTimePosition;
	int satTimeLength;
	int satAttPosition;
	int satAttLength;
	int satGPSPosition;
	int satGPSLength;
	int satTimeCheckPosition;
	int satTimeCheckLength;
	int GetCol();
	static int framecount;
	void FrameNum(const int& _framenum)
	{
		framenum = _framenum;
	}

	static void AddFrame()
	{
		framecount++;
	}

	static void ClearFrame()
	{
		framecount = 0;
	}

	///���Ǹ������õĽӿ�,�����ǳ����Ͻ���ȡ��Ҳ���Կ�����
	virtual bool IGiveMeFive(char* _auxcontent);
	OrbitAPI& Orbit();
	AttitudeAPI& Attitude();
	LineTimeAPI& LineTime();

	virtual bool IEquipment(vector<unsigned char> _information);
	virtual void SetUTCAPIByTimestring() {}
protected:
	void SetPropertys();

};

