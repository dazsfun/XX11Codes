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
	///����ʹ��Ĭ�Ϲ��캯����̰���˻���Сbug
	///�ñ��������������㣬����û��������ʼ��
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
		///���ӷ���ֱ����Ҫȡ�����Envelop������ݲ����ˣ����߷�Ҫ������д����Ч�����ˣ��ŷ���ռ䡣
		///������ЧFRED����ռ�ò���Ҫ���ڴ�ռ�
		///_infobye = new char[length];
	}

public:
	///if you want this class to be classic,you gotta work much harder

	///���Ӱ������
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
	//������û�в�����Ĭ�Ϲ��캯��
	//���������ļ����жϻ�����ȡ����
	EnvelopeAPI(string configInfo, string fredinfo);
	//�����û�������������Ķ���Ļ�����ȡ����
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
			IAddLog("�����������ڵ�֡���Ϊ:%d", _auxindex[i]);
		}

		///�Ƚϸ��������Ƿ����ظ�
		for (int i = 0; i < _auxindex.size() - 1; i++)
		{
			int indexofaux = _auxindex[i];
			int indexofnexaux = _auxindex[i + 1];
			bool contentcompare = safestring::compare_s(_content[indexofaux]._dataformat._infobyte + 3, _content[indexofnexaux]._dataformat._infobyte + 3, 144);
			bool contentcompare1 = safestring::compare_s(_content[indexofaux]._dataformat._infobyte + 3, _content[indexofaux]._dataformat._infobyte + 3, 144);
			if (contentcompare1)
			{
				cout << "����ͨ��" << endl;
			}
			if (contentcompare)
			{
				cout << "����������ͬ" << endl;
				cout << "��ͬ�������ݶ�Ӧ�������ŷֱ�Ϊ" << indexofaux << "\t" << indexofnexaux<<endl;
				IAddLog("�������ݳ�����ͬ����ͬ�������ݶ�Ӧ��������Ϊ,%d\t%d", indexofaux, indexofnexaux);
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
			///ʹ�����ַ�ʽ����ʱ������������Ч���� 
			///��Щ������Ч�ɵ����߾���������ֱ��ִ������
			imgoutputcount++;
			return true;
		}
		for (int i = 0; i < _imgframeheight; )
		{
			if (imgoutputcount >= _content.size())
			{
				cout << "FRED�������Ѷ�ȡ����Ӱ�������Ѿ������ϣ�" << endl;
				IAddLog("FRED�������Ѷ�ȡ����Ӱ�������Ѿ������ϣ�");
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

	///���ݳ�ʼ��������Ѱ��Fred�������ļ������û���ҵ�������н���ʽ����
	ErrorType IGeneration();
	ErrorType ILocate(string _filepath);

private:
	string objectname;
	//����Satellite Fred���ݣ��ӿ��ļ��ж����һ�����ֳ������������ֽ���
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
	///����objectname��ȡ��Ӧ�������ļ����Ӷ��Զ�Ӧ�ͺŵ�ԭʼ���ݽ��ж�ȡ
	virtual ErrorType IStreamIn(string inputfilepath);

	//�Ӷ�ȡ�õ���ԭʼ���ݻ����У��õ���Ӧ�ǵı�׼���ݰ�����,���麯�����������ֱ�ʵ�����ģ��
	virtual ErrorType IConduct() = 0;

	//�����ݰ������У�����ö����ݵ���̬���������ʱ�ȱ�Ҫ����Ϣ
	virtual ErrorType IMetaAux();

	///���������ļ�����߼��־������ģ��ͬ���ǽ�����������ȥʵ�ֵ�
	virtual ErrorType IZone() const =0;

	///��������о�  ������ÿ����Ӧ�ĸ�������
	virtual ErrorType IDivision();
};