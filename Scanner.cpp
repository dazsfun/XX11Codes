#include "stdafx.h"
#include "Scanner.h"

bool EnvelopeAPI::Is_Inited()
{
	return _inited;
}

//"{\"bytecount\":\"12\",\"standardstring\":\"AAA AAA AAA AAA\",\"headercount\":\"4\",\"sub_framelength\",\"484\",\"times\":\"5000\"}";
void EnvelopeAPI::SetPropertys()
{
	SetProperty("bytecount", asInt, &bytecount);
	SetProperty("standardstring", asString, &standardstring);
	SetProperty("headercount", asInt, &headercount);
	SetProperty("sub_framelength", asInt, &subframelength);
	SetProperty("times", asInt, &times);
}

bool EnvelopeAPI::TransInput(const char* tempchar, char* outputer, int ncount, int singlecount)
{
	IAddLog("��ʼת����ʼ��־����");
	if (safestring::FredBytes2char(tempchar, ncount, outputer, singlecount))
	{
		IAddLog("��ֹ��־����%sת���ɹ�", tempchar);
		_inited = true;
		return true;
	}

	else
	{
		IAddLog("����%sת��ʧ�ܣ�������־�еĴ�����Ϣ", tempchar);
		printf("����%sת��ʧ�ܣ�������־�еĴ�����Ϣ", tempchar);
		return false;
	}

}

bool EnvelopeAPI::GetInput(int bufferlength, int singlecount)
{
	int ncounter = 0;
	char endmarker[1024] = "end_register";
	do
	{
		char tempreader[1024];
		cin.getline(tempreader, 1024);

		if (safestring::compare(const_cast<char*>(tempreader), const_cast<char*>(endmarker), 12))
		{
			cout << "����������" << endl << "FRED��������..." << endl << "����Ϊ�����ļ����" << endl;
			return true;
		}

		//string heihei = "{\"bytecount\":\"12\",\"standardstring\":\"AAA AAA AAA AAA\",\"headercount\":\"4\",\"sub_framelength\":\"484\",\"times\":\"5000\"}";

		if (!DeSerialize(tempreader))
		{
			IAddLog("��񲻷���Ҫ�󣬴��û�����������");
			cout << "��񲻷���Ҫ������������" << endl << "�������Ϊ:" << endl;
			cout << "{\"bytecount\":\"12\",\"standardstring\":\"AAA AAA AAA AAA\",\"headercount\":\"4\",\"sub_framelength\",\"484\",\"times\":\"5000\"}" << endl;
		}
		else if (subframelength != bufferlength)
		{
			IAddLog("��񳤶Ȳ�����һ����Ҫ��Ҫ����Ϊ%d,�û��ṩ���Ϊ%d", bufferlength, subframelength);
			cout << "��֡������ͷ��֡���뱣��һ�£�ͷ��֡����Ϊ" << bufferlength << "��ǰ���Ϊ:" << subframelength << endl;
			cout << "����������" << endl;
			subframelength = bufferlength;
		}
		else
		{
			int heisize = ((singlecount * headercount % 2) == 0) ? singlecount * headercount / 2 : singlecount * (headercount+1) / 2;
			char* _tempcharoutput = new char[heisize];
			memset(_tempcharoutput, 0, sizeof(char)*heisize);
			if (!TransInput(standardstring.c_str(), _tempcharoutput, headercount, singlecount))
			{
				IAddLog("ͷ����%sת�����󣬵ȴ��û�����������", standardstring.c_str());
				cout << "ͷ����" << standardstring.c_str() << "ת��ʧ��,����ݳ�����Ϣ����¼����" << endl;
			}
			else
			{
				ncounter++;
				char tempchar1[256];
				sprintf_s(tempchar1, "%d", ncounter);
				SafeConfig(_enveloptype, tempchar1);
				Envelop tempev(headercount, _tempcharoutput, subframelength - headercount, singlecount * 4);
				EnvelopeDescription tempds(tempev, times);
				_description.push_back(tempds);
			}
		}
	} while (true);
}

bool EnvelopeAPI::SafeConfig(string types_, string id)
{
	string hei = Serialize();
	ConfigAPI capi(types_, id, hei);
	return true;
}

EnvelopeAPI::EnvelopeAPI(string configinfo, string fredinfo)
{
	SetPropertys();
	_inited = false;
	imgoutputcount = 0;
	auxprintnum = 0;
	///��ȡConfig�ļ�
	ifstream ifs;
	_fredpath = fredinfo;
	ifs.open(configinfo.c_str(), ios::in);

	if (!ifs.is_open())
	{
		cout << "ָ���������ļ������ڣ�������ָ�������ļ�·�������߽���ע��" << endl;
		IAddLog("ָ���������ļ������ڣ�������ָ�������ļ�·�������߽���ע��");
		return;
	}

	do
	{
		char reader[1024];
		ifs.getline(reader, 1024);
		ifs.getline(reader, 1024);
		if (!DeSerialize(reader))
		{
			cout << "�����ļ����𻵣�����������½���ע��" << endl;
			return;
		}

		int singlecount = bytecount / 4;
		if (bytecount % 4 != 0)
		{
			cout << "FRED�ļ������ֳ���λ��������4�ı����������ļ���Ϣ��������������½���ע��" << endl;
		}

		int heisize = ((singlecount * headercount % 2) == 0) ? singlecount * headercount / 2 : singlecount * headercount / 2 + 1;
		char* _tempcharoutput = new char[heisize];
		if (!TransInput(standardstring.c_str(), _tempcharoutput, headercount, singlecount))
		{
			IAddLog("ͷ����%sת�����󣬵ȴ��û�����������", standardstring.c_str());
			cout << "ͷ����" << standardstring.c_str() << "ת��ʧ��,����ݳ�����Ϣ����¼����" << endl;
		}
		else
		{
			Envelop tempev(headercount, _tempcharoutput, subframelength - headercount, bytecount);
			EnvelopeDescription tempds(tempev, times);
			if (times > 1)
			{
				_imgframeheight = times;
				_imgframewidth = tempev.IGetWidth();
			}
			_description.push_back(tempds);
		}
		ifs.getline(reader, 1024);
	} while (ifs.peek() != EOF);

	ifs.close();

	_inited = true;
}

int EnvelopeAPI::IWidth()
{
	return _imgframewidth;
}

int EnvelopeAPI::IHeight()
{
	return _imgframeheight;
}

int EnvelopeAPI::IByteCount()
{
	return bytecount;
}

EnvelopeAPI::EnvelopeAPI(string fredinfo, string enveloptypes, bool GetClientinput)
{
	///�趨�ӿ��ļ�����
	_enveloptype = enveloptypes;
	SetPropertys();
	_inited = false;
	_fredpath = fredinfo;
	char tempchar[1024];
	cout << "�������ǰҪ�����FRED���ı�ʶ�룬��ʶ��Ϊ��д��ĸ�����֣��»�����ɵ��ַ���" << endl;
	cin.getline(tempchar, 1024);
	_enveloptype = tempchar;
	cout << "��������ʼ��־���" << endl;

	cin.getline(tempchar, 1024);
	cout << "��������ʼ��־�ֳ�" << endl;
	int ncount = 0;
	cin >> ncount;
	cout << "������һ�ֳ���Ӧ��λ��������Ϊ4�ı�����" << endl;
	int nbyte = 0;
	cin >> nbyte;

	int singlecount = nbyte / 4;
	int counter = singlecount * ncount / 2;
	char* outputer = new char[counter];
	TransInput(tempchar, outputer, ncount, singlecount);
	cout << "������һ����֡�ĳ��ȣ�����ͷ�ڵ�ĳ��ȣ�" << endl;
	int bufferlength = 0;
	cin >> bufferlength;

	Envelop startEv(counter, outputer, bufferlength - counter, nbyte);
	EnvelopeDescription starthei(startEv, 1);
	///char tempchar1[1024];
	///sprintf_s(tempchar, "%d", 0);
	///SafeConfig(_enveloptype, tempchar1);
	_description.push_back(starthei);

	cout << "�밴�����µĸ�ʽ�������µ�ÿ����Ҫ��֡�ĸ�ʽ���Լ����ǽ������ظ��ĳ���" << endl;
	cout << "���磺 һ�ֳ�Ϊ12λ�� ��־ͷ�ַ���Ϊ\"AAA AAA AAA AAA\",��־ͷ�ַ����������ֳ���λ4��,һ��֡�ܳ���Ϊ��484�����ָ�ʽ����֡���ظ�5000��" << endl;
	cout << "���Ӧ��cmdָ��Ϊ:" << endl;
	cout << "{\"bytecount\":\"12\",\"standardstring\":\"AAA AAA AAA AAA\",\"headercount\":\"4\",\"sub_framelength\",\"484\",\"times\":\"5000\"}";
	cout << "�ظ� ���밴�����ϸ�ʽָ����ָ�ֱ���ļ�����������" << endl;
	cout << "�������ʱ������end_register" << endl;
	cout << "ע1��������Լ��1��������֡�ģ������Ѿ����������֡����֡������Ӧ����ͬ���������㣬FRED����彫ʧ��" << endl;
	cout << "ע2����������ΪFRED����ע�ᣬע��ɹ��󣬻���ָ��Ŀ¼�������ض��������ļ�����ͬ���ǵ�FRED��������ظ�ע�ᣬ�ɲ��ö�ȡ�����ļ��Զ�����FRED���" << endl;

	if (GetInput(bufferlength, singlecount))
	{
		IAddLog("��������...��");
		cout << "�������ϣ�" << endl;
	}
}

ErrorType EnvelopeAPI::IJson(EnvelopeDescription* _descriptionObject, string _chardescription)
{
	return Process_Success;
}

ErrorType EnvelopeAPI::IFindStart()
{
	ifstream ifs;
	ifs.open(_fredpath.c_str(), ios::in | ios::_Nocreate | ios::binary);

	if (!ifs.is_open())
	{
		IAddLog("%s�ļ������ڣ�FredAPI��ʼ�� ʧ��", _fredpath.c_str());
		return File_Invalid_Existence_Fred;
	}

	char reader;
	while (ifs.peek() != EOF)
	{
		ifs >> reader;  //һ��ʼֻ��һ��һ��char���ƶ�

		///����ҵ��˻��з� ����ʾ����
		if (reader == '\x0d')    ///ANSI��ʽ�»��з��������ݲ���������FRED�ļ�������õĸ�ʽ
		{
			ifs >> reader;
			if (reader == '\x0a')
			{
				//����
				IAddLog("Fred�ļ��������˻��з����ļ��������ƻ���λ��Ϊ%d", ifs.cur - 2);
				IAddLog("�����˳�");
				return Fred_Unexpected_End;
			}
			else
				ifs.seekg(-1, ifs.cur);
		}

		//����ҵ�������֡�Ŀ�ͷ�����н�һ���ж�
		if (reader == _description[0]._dataformat._header[0])
		{
			ifs.seekg(-1, ifs.cur);
			int reallength = static_cast<int>(_description[0]._dataformat.singlebyte / 8.0 * _description[0]._dataformat._headerlength);
			char* tempreader = new char[reallength];
			ifs.read(tempreader, reallength);
			if (safestring::compare(tempreader, _description[0]._dataformat._header, reallength))
			{
				_startPos = static_cast<long>(ifs.cur);
				ifs.close();
				delete[] tempreader;
				return Process_Success;
			}
			delete[] tempreader;
		}
	}

	IAddLog("��û����Fred�ļ�:%s���ҵ��涨���ļ�ͷ��Ϣ���ļ�ͷ��ϢΪ��");
	return Fred_Cannot_Find_Header;
}

ErrorType FREDStructure::ILocate(string _filepath)
{
	ifstream ifs;
	ifs.open(_filepath.c_str(), ios::in | ios::_Nocreate | ios::binary);

	if (!ifs.is_open())
	{
		cout << "FRED�ļ������ڣ�����м��" << endl;
		throw(exception("FRED�ļ������ڣ�����м��"));
		return File_Invalid_Existence_Fred;
	}

	return Process_Success;
}

long EnvelopeAPI::GetStartPos()
{
	return _startPos;
}

ErrorType EnvelopeAPI::IScan()
{
	ifstream ifs;
	ifs.open(_fredpath.c_str(), ios::binary);
	if (!ifs.is_open())
	{
		IAddLog("Fred�ļ�%s�����ڻ�����", _fredpath.c_str());
		printf("Fred�ļ�%s�����ڻ�����", _fredpath.c_str());
		return File_Invalid_Existence_Fred;
	}

	_endPos = _startPos - 1;
	//ifs.seekg(_endPos, ios::beg);

	int sizeofit = _description.size();
	for (int i = 0; i < sizeofit; i++)
	{
		cout << "Scan����:" << i + 1 << "/" << sizeofit << "\r";
		if (!(IScan(ifs, i) == Process_Success))
		{
			return Fred_Cannot_Find_Header;
		}

		if (_description[i].Length == 1)
		{
			//ICheckContent();
		}
		_endPos = ifs.cur;
	}

	_endPos = ifs.cur;
	ifs.close();
	return Process_Success;
}

bool FredEnvelop::IFormAux(AuxFactory& _factory)
{
	///�и�Fred��Ϣ
	vector<unsigned char> heiinfo;
	safestring::TrimBytes(_infobyte, heiinfo, _length, 12);
	_factory.IEquipment(heiinfo);
	return true;
}

///Make sure imgcontent's memory has been located
bool FredEnvelop::IFormImage(unsigned int* imgcontent)
{
	///����single count �ж���ζ�ȡ����
	double _readstep = static_cast<double>(singlebyte) / 8.0;
	int realstep = static_cast<int>(_readstep + 0.5);

	if (imgcontent == nullptr)
	{
		cout << "����ԭʼ���ݵ����������" << endl;
		IAddLog("����ԭʼ���ݵ����������");
	}

	int totalstep = static_cast<int>(static_cast<double>(realdatawidth) *  _readstep);
	memset(imgcontent, 0, realdatawidth * sizeof(unsigned int));
	for (int i = 0; i < realdatawidth; i++)
	{
		double _startnum = i * _readstep;
		unsigned int _startindex = static_cast<unsigned int>(_startnum);
		int tempint;
		if (_readstep != static_cast<unsigned int>(_readstep)&& i %2 == 0)
			tempint = -12;
		else
			tempint = -8;
		///Ȼ���ж��Ǵӵ�λ���ǴӸ�λ��ʼ��ȡ���Լ��Ⱥ�˳��
		for (unsigned int j = _startindex; j < (_startindex + realstep); j++)
		{
			int moveleftsecond = tempint + 8 * (_startindex + realstep -j);// ///�ڶ������Ƶ�λ��
			moveleftsecond = (moveleftsecond <= 0) ? 0 : moveleftsecond;
			int moverightfirst; ///��һ�����Ƶ�λ��
			int moveleftfirst;// ��һ�����Ƶ�λ��


			if (_startnum != _startindex&& j == _startindex)
			{
				moveleftfirst = 4;
				moverightfirst = 4;
			}
			///����ж������д���
			//else if (j == (_startindex + realstep -1) && (_readstep != static_cast<int>(_readstep)))
			else if (j == (_startindex + realstep - 1) && (i % 2 == 0) && (_readstep != static_cast<int>(_readstep)))
			{
				moveleftfirst = 0;
				moverightfirst = 4;
			}
			else
			{
				moveleftfirst = moverightfirst = 0;
			}

			unsigned char heiuchar = _infobyte[j];
			unsigned char heiucharnew = (heiuchar) << moveleftfirst;
			//checked to be ok
			///A bug detected here
		    ///20150220
			imgcontent[i] += static_cast<unsigned int>(heiucharnew >> moverightfirst) << moveleftsecond;
			//imgcontent[i] = 
		}
		if (imgcontent[i] > pow(2, singlebyte))
		{
			cout << "�������ݴ��ڽӿڶ��壬Error Code: Fatal0E301" << endl;
			IAddLog("�������ݴ��ڽӿڶ��壬Error Code: Fatal0E301");
			return false;
		}
	}

	return true;
}

ErrorType EnvelopeAPI::IScan(ifstream& ifs, int id)
{
	EnvelopeDescription hei(_description[id]);

	///�ж��Ƿ���Ҫ��1
	bool needtobeadded = false;
	double readheaderlength;
	if ((hei._dataformat.singlebyte / 4 * hei._dataformat._headerlength) % 2 != 0)
	{
		needtobeadded = true;
	}

	int times = hei.Length;
	int headerlength = hei._dataformat._headerlength;
	int bodylength = hei._dataformat._length;

	double rate = static_cast<double>(hei._dataformat.singlebyte) / 8;


	if (needtobeadded)
	{
		readheaderlength = headerlength * rate;
		headerlength = static_cast<int>(rate * (headerlength + 1));
		bodylength = static_cast<int>(rate * (bodylength - 1));
	}

	else
	{
		headerlength = static_cast<int>(rate * headerlength);
		bodylength = static_cast<int>(rate * bodylength);
	}

	char* tempreader = new char[headerlength + bodylength];
	char recorder[1024];



	//if (needtobeadded)
		//times *= 2;
	for (int i = 0; i < times; i++)
	{
		ifs.read(tempreader, headerlength + bodylength);
		//ifs.seekg((headerlength + bodylength), ios::cur);  we dont need this line if we already define ios::app
		if (!needtobeadded)
		{
			
			if (!safestring::compare(tempreader, hei._dataformat._header, headerlength))
			{
				ifs.seekg(-(headerlength + bodylength), ios::cur);
				IAddLog("������FRED�ļ���ָ��λ��%dû�з����ļ�ͷ%s", ifs.cur, hei._dataformat._headerdescprtion.c_str());
				printf("������FRED�ļ���ָ��λ��%dû�з����ļ�ͷ%s", ifs.cur, hei._dataformat._headerdescprtion.c_str());
				delete[] tempreader;
				return Fred_Cannot_Find_Header;
			}
		}
		else
		{
			if (!safestring::compare(tempreader, hei._dataformat._header, readheaderlength,true))
			{
				ifs.seekg(-(headerlength + bodylength), ios::cur);
				IAddLog("������FRED�ļ���ָ��λ��%dû�з����ļ�ͷ%s", ifs.cur, hei._dataformat._headerdescprtion.c_str());
				printf("������FRED�ļ���ָ��λ��%dû�з����ļ�ͷ%s", ifs.cur, hei._dataformat._headerdescprtion.c_str());
				delete[] tempreader;
				return Fred_Cannot_Find_Header;
			}
		}

		EnvelopeDescription des = _description[id];

		if (des._dataformat._infobyte != nullptr)
		{
			delete[] des._dataformat._infobyte;
			des._dataformat._infobyte = nullptr;
		}

		des._dataformat._infobyte = new char[bodylength];
		memcpy_s(des._dataformat._infobyte, bodylength, tempreader + headerlength, bodylength);
		_content.push_back(des);

	}

	delete[] tempreader;
	return Process_Success;
}

AuxFactory EnvelopeAPI::IFormAux()
{
	///ͬ��������ط��������ݶ���������ļ�·�����̶��� �պ���Ҫ�޸�
	AuxFactory auxfct("auxBinaryconfig");

	if (auxprintnum >= _auxindex.size())
	{
		cout << "���������Ѿ�ȫ��������" << endl;
		IAddLog("�������ݷ������");
	}

	cout << _auxindex[auxprintnum] << endl;
	//cout << (int)(_content[_auxindex[auxprintnum]]._dataformat._infobyte[100]) << endl;
	Envelop env = _content[_auxindex[auxprintnum]]._dataformat;
	env.IFormAux(auxfct);
	auxprintnum++;
	return auxfct;
}