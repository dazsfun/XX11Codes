#include "stdafx.h"
#include <direct.h>
#include "jsonSerializer.h"


string JSONOBJECTSerializer::Serialize(string dirpath, bool bfile)
{
	string hei = dirpath + "\\" + objectname + ".json";
	ofstream ofs(hei.c_str(), ios::app|ios::ate);
	if (!ofs.is_open())
	{
		string dir = hei.substr(0, hei.rfind("\\"));
		if (!safestring::createDir(dir)) return "InvalidDir";
		ofs.open(hei.c_str(), ios::app | ios::ate);
		if (!ofs.is_open())
		{
			cout << "创建" << objectname.c_str() << "序列化结果失败" << endl;
			cout << "请检查" << endl;
			IAddLog("创建%s序列化结果失败，请检查", objectname.c_str());
			return "InValid";
		}
	}
	ofs << Serialize().c_str();
	ofs.close();
	return hei;
}

bool JSONOBJECTSerializer::operator<(const JSONOBJECTSerializer & objectinfo) const
{
	return ((this->utcapi) < (objectinfo.utcapi));
}

bool JSONOBJECTSerializer::operator>(const JSONOBJECTSerializer & objectinfo) const
{
	return ((this->utcapi) > (objectinfo.utcapi));
}

bool JSONOBJECTSerializer::operator<=(const JSONOBJECTSerializer& objectinfo) const
{
	return ((this->utcapi) <= (objectinfo.utcapi));
}

bool JSONOBJECTSerializer::operator>=(const JSONOBJECTSerializer& objectinfo) const
{
	return ((this->utcapi) >= (objectinfo.utcapi));
}

bool JSONOBJECTSerializer::operator==(const JSONOBJECTSerializer& objetinfo) const
{
	return ((this->utcapi) == (objetinfo.utcapi));
}


///序列化过程
string JSONOBJECTSerializer::Serialize() const
{
	Json::Value new_item;
	int nSize = m_listName.size();
	for (int i = 0; i < nSize; i++)
	{
		void* pAddr = m_listPropertyAddr[i];
		switch (m_listType[i])
		{

		case asInt:
			new_item[m_listName[i]] = (*(int*)pAddr);
			break;
		case asUInt:
			new_item[m_listName[i]] = (*(unsigned int*)pAddr);
			break;
		case asString:
			new_item[m_listName[i]] = (*(string*)pAddr);
			break;
		case asDouble:
			new_item[m_listName[i]] = (*(double*)pAddr);
			break;
		default:
			break;
		}
	}

	Json::FastWriter writer;
	string out2 = writer.write(new_item);
	return out2;
}

void JSONOBJECTSerializer::SetProperty(int index, CEumJsonTypeMap type, void* addr)
{
	m_listName.push_back(MetaClass::Register[index]);
	m_listPropertyAddr.push_back(addr);
	m_listType.push_back(type);
}

void JSONOBJECTSerializer::SetProperty(string name, CEumJsonTypeMap type, void *addr)
{
	m_listName.push_back(name);
	m_listPropertyAddr.push_back(addr);
	m_listType.push_back(type);
}

ErrorType JSONOBJECTSerializer::XmlDeSerialization(const char* str)
{
	int sizeofit = m_listName.size();
	//cout << _xmlRegex.c_str() << endl;
	int* _regexgroup = new int[sizeofit];
	for (int i = 0; i < sizeofit; i++)
	{
		_regexgroup[i] = (2 * i + 1);
	}

	vector<string> _results = safestring::GetRegexResult(str, _xmlRegex, _regexgroup, sizeofit);

	Function2(_results);
	delete[] _regexgroup;
	return Process_Success;
}


///暂时没有时间了，但是用autoptr或者模版是不可以支持更多类型
///答案应该是肯定的，至少可以支持自定义类型，以后再来重构
///？？ 什么时候可以支持STL容器的反序列化？
void JSONOBJECTSerializer::Function2(vector<string> _contents)
{
	for (int i = 0; i < _contents.size(); i++)
	{
		CEumJsonTypeMap temptype = m_listType[i];
		if (temptype == asInt)
		{
			void* pAddr = m_listPropertyAddr[i];
			(*(int*)pAddr) = atoi(_contents[i].c_str());
		}
		else if (temptype == asDouble)
		{
			void* pAddr = m_listPropertyAddr[i];
			(*(double*)pAddr) = atof(_contents[i].c_str());
		}
		else if (temptype == asInt64)
		{
			void* pAddr = m_listPropertyAddr[i];
			(*(long*)pAddr) = atol(_contents[i].c_str());
		}
		else if (temptype == asUInt)
		{
			void* pAddr = m_listPropertyAddr[i];
			(*(unsigned int*)pAddr) = static_cast<unsigned int>(atoi(_contents[i].c_str()));
		}
		else if (temptype == asUInt64)
		{
			void* pAddr = m_listPropertyAddr[i];
			(*(unsigned long*)pAddr) = static_cast<unsigned long>(atol(_contents[i].c_str()));
		}
		else if (temptype == asString)
		{
			void* pAddr = m_listPropertyAddr[i];
			(*(string*)pAddr) = _contents[i];
		}
	}
}

void JSONOBJECTSerializer::Function1(int _index1)
{
	CEumJsonTypeMap tempType = m_listType[_index1];

	if (tempType == asString)
	{
		char tempreader[1024];

		sprintf_s(tempreader, 1024, "<%s>(([\\s\\S]\*?))<\/%s>\\s*\\S*\\s*", m_listName[_index1].c_str(), m_listName[_index1].c_str());

		//_xmlRegex = "<ImageMode>";
		_xmlRegex += tempreader;
	}
	else
	{
		char tempreader[1024];
		sprintf_s(tempreader, 1024, "<%s>(-?\\d+(\\.\\d+)?)<\/%s>\\s*\\S*\\s*", m_listName[_index1].c_str(), m_listName[_index1].c_str());
		_xmlRegex += tempreader;
	}

}

CJson_Ojbect_Base_Serializer::CJson_Ojbect_Base_Serializer()
{

}
//<GpsTime>([\s\S]*?)</GpsTime>[\s\S]?<PosX>(-?\d+(\.\d+)?)</PosX>[\s\S]?<PosY>(-?\d+(\.\d+)?)</PosY>[\s\S]?<PosZ>(-?\d+(\.\d+)?)</PosZ>[\s\S]?<VelX>(-?\d+(\.\d+)?)</VelX>[\s\S]?<VelY>(-?\d+(\.\d+)?)</VelY>[\s\S]?<VelZ>(-?\d+(\.\d+)?)</VelZ>[\s\S]?
ErrorType JSONOBJECTSerializer::FormatXml()
{
	_xmlRegex = "[\\s\\S]?[\s\S]?";
	int sizeofit = m_listName.size();

	for (int i = 0; i < sizeofit; i++)
	{
		Function1(i);
	}
	return Process_Success;
}

bool JSONOBJECTSerializer::DeSerialize(const char* str)
{
	if (_DOMType == asXml)
	{
		XmlDeSerialization(str);
		SetUTCAPIByTimestring();
		return true;
	}
	Json::Reader reader;
	Json::Value root;

	if (reader.parse(str, root))
	{
		int nSize = m_listName.size();
		for (int i = 0; i < nSize; i++)
		{
			void* pAddr = m_listPropertyAddr[i];

			switch (m_listType[i])
			{
			case asInt:
			{
						 // string temp = root.get(m_listName[i], "").asString();
						 //int heitemp;  sscanf_s(temp.c_str(), "%d", &heitemp);
						  //int heitemp = root.get(m_listName[i], "").asInt();
						  (*(int*)pAddr) = root.get(m_listName[i], 0).asInt();
			}
				break;
			case asUInt:
			{
						   //string temp = root.get(m_listName[i], "").asString();
						   //unsigned int heitemp;  sscanf_s(temp.c_str(), "%d", &heitemp);
						   (*(unsigned int*)pAddr) = root.get(m_listName[i], 0).asUInt();
			}
			case asString:
				(*(string*)pAddr) = root.get(m_listName[i], "").asString();
				break;
			case asDouble:
			{
							 //string temp = root.get(m_listName[i], "").asString();
							 //double dheitemp;  sscanf_s(temp.c_str(), "%lf", &dheitemp);
							 (*(double*)pAddr) = root.get(m_listName[i], "").asDouble();
			}
			default:
				break;
			}
		}
		//SetUTCAPIByTimestring();
		return true;
	}

	return false;
}

void JSONOBJECTSerializer::ClearPropertys()
{
	m_listName.clear();
	m_listPropertyAddr.clear();
	m_listType.clear();
}