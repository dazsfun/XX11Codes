#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "Log.h"
#include "jsonSerializer.h"
#include "Safestring.h"
#include "AuxFactory.h"
using namespace std;


typedef class FredEnvelop : public LogAPI
{
public:
	int singlebyte;
	int _headerlength;
	char* _header;
	int _length;
	char* _infobyte;
	string _headerdescprtion;

private:
	int realdatawidth;

public:
	int IGetWidth()
	{
		return realdatawidth;
	}

	string formstring()
	{
		int realength = static_cast<double>(singlebyte / 8.0) * _length;
		string returnvalue = "";
		char reader[1024];
		for (int i = 0; i < realength; i++)
		{
			unsigned char hei = _infobyte[i];
			sprintf_s(reader, "%x", hei);
			returnvalue += reader;
		}

		return returnvalue;
	}
public:
	///避免使用默认构造函数，贪便宜会有小bug
	///让编译器帮助提醒你，变量没有完整初始化
	/*
	FredEnvelop()
	{
		_header = nullptr;

	}
	*/
	FredEnvelop(const FredEnvelop& newinfo)
	{
		*this = newinfo;		
	}

	void operator=(const FredEnvelop& newinfo)
	{
		singlebyte = newinfo.singlebyte;
		_headerlength = newinfo._headerlength;
		_length = newinfo._length;
		_headerdescprtion = newinfo._headerdescprtion;

		bool needtobeadded = false;
		if (((singlebyte / 4) * (_headerlength)) % 2 != 0)
			needtobeadded = true;

		double rate = static_cast<double>(singlebyte) / 8.0;
		int realheaderlength, reallength;
		if (!needtobeadded)
		{
			realheaderlength = static_cast<int>(rate * _headerlength);
			reallength = static_cast<int>(rate * _length);
		}
		else
		{
			realheaderlength = static_cast<int>(rate * (_headerlength + 1));
			reallength = static_cast<int>(rate * (_length-1));
		}

		if (newinfo._header != nullptr)
		{
			_header = new char[realheaderlength];
			memcpy(_header, newinfo._header, sizeof(char)* realheaderlength);
		}

		else
			_header = nullptr;

		if (newinfo._infobyte != nullptr)
		{
			_infobyte = new char[reallength];
			memcpy(_infobyte, newinfo._infobyte, sizeof(char)* (reallength));
		}

		else
			_infobyte = nullptr;

		realdatawidth = reallength / rate;
	}

	FredEnvelop(int hlength, char* cheader, int bodylength, int singlebytelength) :_headerlength(hlength), _header(cheader), _length(bodylength), singlebyte(singlebytelength)
	{
		if (((singlebyte / 4)*hlength) % 2 != 0)
		{
			realdatawidth = (bodylength - 1);
		}
		_infobyte = nullptr;
	}

	ErrorType ISet(int headerlength, int length, char* header)
	{
		_length = length;
		_headerlength = headerlength;
		_header = header;
		///拖延法则，直到非要取用这个Envelop里的数据不可了，或者非要向其中写入有效数据了，才分配空间。
		///避免无效FRED数据占用不必要的内存空间
		///_infobye = new char[length];
	}

public:
	///if you want this class to be classic,you gotta work much harder

	///组成影像内容
	bool IFormImage(unsigned int* imgcontent);
	bool IFormAux(AuxFactory& auxfactory);
	~FredEnvelop()
	{
		if (_header != nullptr)
		{
			delete[] _header;
			_header = nullptr;
		}
		if (_infobyte != nullptr)
		{
			delete[] _infobyte;
			_infobyte = nullptr;
		}
		
	}
}Envelop;

struct EnvelopeDescription
{

	EnvelopeDescription(Envelop _data, int dlength) :_dataformat(_data), Length(dlength)
	{

	}

	
	EnvelopeDescription(const EnvelopeDescription& newinfo) :_dataformat(newinfo._dataformat)
	{
		Length = newinfo.Length;
	}
	

	void operator= (const EnvelopeDescription& newinfo)
	{
		_dataformat = newinfo._dataformat;
		Length = newinfo.Length;
	}

	Envelop _dataformat;
	int Length;
};



class EnvelopeAPI : public JSONOBJECTSerializer
{
public:
	//不定义没有参数的默认构造函数
	//根据配置文件来判断基本读取规则
	EnvelopeAPI(string configInfo, string fredinfo);
	//根据用户的输入来定义改对象的基本读取规则
	EnvelopeAPI(string fredinfo, string enveloptp,bool GetClientinput);
	virtual ~EnvelopeAPI() {}
	bool Is_Inited();

public:
	void ICheckContent()
	{
		for (int i = 0; i < _content.size(); i++)
		{
			if (_content[i].Length == 1)
			{
				string hei = _content[i]._dataformat.formstring();
				//cout << i << endl;
				//cout << hei.c_str() << endl;
				IAddLog("%d", i);
				IAddLog(hei.c_str());
			}
		}
	}

	virtual void SetUTCAPIByTimestring(){};
private:
	string _enveloptype;
	string _fredpath;
	vector<EnvelopeDescription> _description;
	vector<EnvelopeDescription> _content;
	FredEnvelop* satalliteMetaInfo;
	long _startPos;
	long _endPos;
	int length;
	bool _inited;
	bool TransInput(const char* tempchar, char* outputer, int ncount, int singlecount);
	bool GetInput(int bufferlength,int);
	bool SafeConfig(string types_,string id_);

	int auxprintnum;
	vector<int> _auxindex;
	int imgoutputcount;
public:
	int IAuxNum()
	{
		int resultnum = 0;
		for (int i = 0; i < _content.size(); i++)
		{
			if (_content[i].Length == 1)
			{
				_auxindex.push_back(i);
				resultnum++;
			}
		}

		for (int i = 0; i < _auxindex.size(); i++)
		{
			IAddLog("辅助数据所在的帧编号为:%d", _auxindex[i]);
		}

		///比较辅助数据是否有重复
		for (int i = 0; i < _auxindex.size() - 1; i++)
		{
			int indexofaux = _auxindex[i];
			int indexofnexaux = _auxindex[i + 1];
			bool contentcompare = safestring::compare_s(_content[indexofaux]._dataformat._infobyte + 3, _content[indexofnexaux]._dataformat._infobyte + 3, 144);
			bool contentcompare1 = safestring::compare_s(_content[indexofaux]._dataformat._infobyte + 3, _content[indexofaux]._dataformat._infobyte + 3, 144);
			if (contentcompare1)
			{
				cout << "测试通过" << endl;
			}
			if (contentcompare)
			{
				cout << "辅助数据相同" << endl;
				cout << "相同辅助数据对应的索引号分别为" << indexofaux << "\t" << indexofnexaux<<endl;
				IAddLog("辅助数据出现相同，相同辅助数据对应的索引号为,%d\t%d", indexofaux, indexofnexaux);
			}
		}

		return resultnum;
	}

	int IFrameNum()
	{
		int resultnum = 0;
		for (int i = 0; i < _content.size(); i++)
		{
			if (_content[i].Length != 1)
			{
				resultnum++;
			}
		}

		return resultnum / _imgframeheight;
	}

	bool IFormImage(unsigned int* content)
	{
		if(content == nullptr)
		{
			///使用这种方式调用时，将会跳过无效数据 
			///哪些数据无效由调用者决定，本类直接执行任务
			imgoutputcount++;
			return true;
		}
		for (int i = 0; i < _imgframeheight; )
		{
			if (imgoutputcount >= _content.size())
			{
				cout << "FRED中所有已读取到的影像数据已经输出完毕！" << endl;
				IAddLog("FRED中所有已读取到的影像数据已经输出完毕！");
				return true;
			}
			if (_content[imgoutputcount].Length == 1)
			{
				imgoutputcount++;
				continue;
			}

			_content[imgoutputcount]._dataformat.IFormImage(content + i * _imgframewidth);
			imgoutputcount++;
			i++;
		}
		return true;
	}

	AuxFactory IFormAux();

	long GetStartPos();
	long GetEndPos();
	int FrameForOneScene()
	{
		return _FrameNum_One_Scene;
	}

	void FrameForOneScene(int num)
	{
		_FrameNum_One_Scene = num;
	}
public:
	ErrorType IJson(EnvelopeDescription* _descriptionObject, string _chardescription);
	ErrorType IFindStart();

	ErrorType IScan();
	ErrorType IRead();
	virtual ErrorType ICheck(){ return Process_Success; };

	int IWidth();
	int IHeight();
	int IByteCount();

protected:
	ErrorType IScan(ifstream& ifs, int id);
	ErrorType FrameScan(ifstream& ifs);

	//"{\"bytecount\":\"12\",\"standardstring\":\"AAA AAA AAA AAA\",\"headercount\":\"4\",\"sub_framelength\",\"484\",\"times\":\"5000\"}";
	int bytecount;
	int headercount;
	int times;
	int subframelength;

	int _imgframewidth;
	int _imgframeheight;
	string standardstring;
	virtual void SetPropertys();

private:
	int _FrameNum_One_Scene;

};



class FREDStructure
{
public:
	FREDStructure(string _input) :objectname(_input)
	{

	}

	~FREDStructure()
	{

	}

	///根据初始化的名称寻找Fred的配置文件，如果没有找到，则进行交互式创建
	ErrorType IGeneration();
	ErrorType ILocate(string _filepath);

private:
	string objectname;
	//对于Satellite Fred数据，接口文件中定义的一个“字长”所包含的字节数
	int _bytelenth;
	vector<FredEnvelop> _envelop;
	int* _sequence;
	long _startpos;
	long _endpos;

};

class CScanner
{
public:
	 map<string, string> basicinfomap;

public:
	string objectname;

public:
	///根据objectname获取对应的配置文件，从而对对应型号的原始数据进行读取
	virtual ErrorType IStreamIn(string inputfilepath);

	//从读取得到的原始数据缓存中，得到对应星的标准数据包集合,纯虚函数，子类必须分别实现这个模块
	virtual ErrorType IConduct() = 0;

	//从数据包集合中，输出该段数据的姿态，轨道，行时等必要的信息
	virtual ErrorType IMetaAux();

	///根据配置文件完成逻辑分景，这个模块同样是交给各个子类去实现的
	virtual ErrorType IZone() const =0;

	///完成物理切景  并分配每景对应的辅助数据
	virtual ErrorType IDivision();
};