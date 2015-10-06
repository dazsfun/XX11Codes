#pragma once
#include <string>
#include <vector>
#include "log.h"
#include "safestring.h"
#include "UTCAPI.h"
#include "json\json.h"
using namespace std;

#include "MetaClass.h"
#pragma comment(lib,"json_mtd.lib")


enum CEnumDomType
{
	asXml = 1,
	asJSon,
};
///本类目前支持，对有string，int ，long double 等基本类型数据的
///Xml，Json标准的序列化与反序列化
///目前还不支持STL容器，暂时看也不需要
typedef class CJson_Ojbect_Base_Serializer : public LogAPI
{
public:
	CJson_Ojbect_Base_Serializer();

	 CJson_Ojbect_Base_Serializer(unsigned char* _content)
	{
		 IGiveMeFive(_content);
	}
	///目前支持的类型
	enum CEumJsonTypeMap
	{
		asInt = 1,
		asUInt,
		asString,
		asInt64,
		asUInt64,
		asDouble,
	};

	string Serialize() const;
	string Serialize(string dirpath , bool bfile = true);
	bool DeSerialize(const char* str);

public:
	UTCAPI utcapi;
	bool operator<=(const CJson_Ojbect_Base_Serializer& orbitino) const;
	bool operator>=(const CJson_Ojbect_Base_Serializer& orbitino) const;
	bool operator<(const CJson_Ojbect_Base_Serializer& orbitino) const;
	bool operator>(const CJson_Ojbect_Base_Serializer& orbitino) const;
	bool operator==(const CJson_Ojbect_Base_Serializer& orbitino) const;

	///For test
public:
	bool GetListSize(int& sizeofit)
	{
		int size1 = m_listName.size();
		int size2 = m_listType.size();
		int size3 = m_listPropertyAddr.size();
		if (size1 != size2)
			return false;
		if (size2 != size3)
			return false;

		sizeofit = size1;

		return true;
	}


public:
	void SetType(CEnumDomType _dom)
	{
		_XmlFormated = false;
		_DOMType = _dom;
		SetPropertys();
		if (_dom == asXml)
			FormatXml();
	}

	virtual void ClearPropertys();

public:
	virtual ~CJson_Ojbect_Base_Serializer(void){
		int sizeofit = m_listPropertyAddr.size();
		for (int i = 0; i < sizeofit; i++)
		{
			if (m_listPropertyAddr[i] != nullptr)
				m_listPropertyAddr[i] = nullptr;
		}
	}

private:
	CEnumDomType _DOMType;
	ErrorType XmlDeSerialization(const char*);
	ErrorType FormatXml();

public:
	string objectname;
	void SetProperty(int index, CEumJsonTypeMap, void *addr);

	void SetProperty(string name, CEumJsonTypeMap type, void *addr);

	virtual void SetPropertys() = 0;
	virtual bool IGiveMeFive(unsigned char* _auxcontent) { return false; }
	virtual void SetUTCAPIByTimestring() = 0;
	vector<string> m_listName;
	vector<void*> m_listPropertyAddr;
	vector<CEumJsonTypeMap>  m_listType;
	string _xmlRegex;

protected:
	bool _XmlFormated;
	void Function1(int _index);
	void Function2(vector<string> _contents);
}JSONOBJECTSerializer;