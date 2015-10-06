#pragma once
#include <vector>
#include <string>
#include "json\json.h"
#include "safestring.h"
#include <iostream>
using namespace std;

namespace MetaClass
{
	////涉及多线程时，只使用最基本C类型作为全局变量（包括不允许STL容器，string等 ）
	///Register变量会在初始化时为每种元数据（例如轨道，姿态）预留15种基本属性，不足时
	///注册文件会自动补空为NULLCONTENT
	extern string Register[150];

	bool MetaClass(string _configfile);
	//~MetaClass();
};

