#pragma once
#include <vector>
#include <string>
#include "json\json.h"
#include "safestring.h"
#include <iostream>
using namespace std;

namespace MetaClass
{
	////�漰���߳�ʱ��ֻʹ�������C������Ϊȫ�ֱ���������������STL������string�� ��
	///Register�������ڳ�ʼ��ʱΪÿ��Ԫ���ݣ�����������̬��Ԥ��15�ֻ������ԣ�����ʱ
	///ע���ļ����Զ�����ΪNULLCONTENT
	extern string Register[150];

	bool MetaClass(string _configfile);
	//~MetaClass();
};

