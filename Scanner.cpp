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
	IAddLog("开始转换起始标志代码");
	if (safestring::FredBytes2char(tempchar, ncount, outputer, singlecount))
	{
		IAddLog("起止标志代码%s转换成功", tempchar);
		_inited = true;
		return true;
	}

	else
	{
		IAddLog("代码%s转换失败，请检查日志中的错误信息", tempchar);
		printf("代码%s转换失败，请检查日志中的错误信息", tempchar);
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
			cout << "规格输入完毕" << endl << "FRED规格建立完毕..." << endl << "保存为配置文件完毕" << endl;
			return true;
		}

		//string heihei = "{\"bytecount\":\"12\",\"standardstring\":\"AAA AAA AAA AAA\",\"headercount\":\"4\",\"sub_framelength\":\"484\",\"times\":\"5000\"}";

		if (!DeSerialize(tempreader))
		{
			IAddLog("规格不符合要求，待用户重新输入中");
			cout << "规格不符合要求，请重新输入" << endl << "范例规格为:" << endl;
			cout << "{\"bytecount\":\"12\",\"standardstring\":\"AAA AAA AAA AAA\",\"headercount\":\"4\",\"sub_framelength\",\"484\",\"times\":\"5000\"}" << endl;
		}
		else if (subframelength != bufferlength)
		{
			IAddLog("规格长度不符合一致性要求，要求规格为%d,用户提供规格为%d", bufferlength, subframelength);
			cout << "子帧长度与头子帧必须保持一致，头子帧长度为" << bufferlength << "当前规格为:" << subframelength << endl;
			cout << "请重新输入" << endl;
			subframelength = bufferlength;
		}
		else
		{
			int heisize = ((singlecount * headercount % 2) == 0) ? singlecount * headercount / 2 : singlecount * (headercount+1) / 2;
			char* _tempcharoutput = new char[heisize];
			memset(_tempcharoutput, 0, sizeof(char)*heisize);
			if (!TransInput(standardstring.c_str(), _tempcharoutput, headercount, singlecount))
			{
				IAddLog("头代码%s转换错误，等待用户重新输入中", standardstring.c_str());
				cout << "头代码" << standardstring.c_str() << "转换失败,请根据出错信息重新录入规格" << endl;
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
	///读取Config文件
	ifstream ifs;
	_fredpath = fredinfo;
	ifs.open(configinfo.c_str(), ios::in);

	if (!ifs.is_open())
	{
		cout << "指定的配置文件不存在，请重新指定配置文件路径，或者进行注册" << endl;
		IAddLog("指定的配置文件不存在，请重新指定配置文件路径，或者进行注册");
		return;
	}

	do
	{
		char reader[1024];
		ifs.getline(reader, 1024);
		ifs.getline(reader, 1024);
		if (!DeSerialize(reader))
		{
			cout << "配置文件被损坏，请检查或者重新进行注册" << endl;
			return;
		}

		int singlecount = bytecount / 4;
		if (bytecount % 4 != 0)
		{
			cout << "FRED文件基本字长的位数必须是4的倍数，配置文件信息错误，请检查或者重新进行注册" << endl;
		}

		int heisize = ((singlecount * headercount % 2) == 0) ? singlecount * headercount / 2 : singlecount * headercount / 2 + 1;
		char* _tempcharoutput = new char[heisize];
		if (!TransInput(standardstring.c_str(), _tempcharoutput, headercount, singlecount))
		{
			IAddLog("头代码%s转换错误，等待用户重新输入中", standardstring.c_str());
			cout << "头代码" << standardstring.c_str() << "转换失败,请根据出错信息重新录入规格" << endl;
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
	///设定接口文件类型
	_enveloptype = enveloptypes;
	SetPropertys();
	_inited = false;
	_fredpath = fredinfo;
	char tempchar[1024];
	cout << "请给出当前要输入的FRED规格的标识码，标识码为大写字母与数字，下划线组成的字符串" << endl;
	cin.getline(tempchar, 1024);
	_enveloptype = tempchar;
	cout << "请输入起始标志命令：" << endl;

	cin.getline(tempchar, 1024);
	cout << "请输入起始标志字长" << endl;
	int ncount = 0;
	cin >> ncount;
	cout << "请输入一字长对应的位数（必须为4的倍数）" << endl;
	int nbyte = 0;
	cin >> nbyte;

	int singlecount = nbyte / 4;
	int counter = singlecount * ncount / 2;
	char* outputer = new char[counter];
	TransInput(tempchar, outputer, ncount, singlecount);
	cout << "请输入一个子帧的长度（包括头节点的长度）" << endl;
	int bufferlength = 0;
	cin >> bufferlength;

	Envelop startEv(counter, outputer, bufferlength - counter, nbyte);
	EnvelopeDescription starthei(startEv, 1);
	///char tempchar1[1024];
	///sprintf_s(tempchar, "%d", 0);
	///SafeConfig(_enveloptype, tempchar1);
	_description.push_back(starthei);

	cout << "请按照如下的格式输入余下的每个需要子帧的格式，以及它们将连续重复的长度" << endl;
	cout << "例如： 一字长为12位， 标志头字符串为\"AAA AAA AAA AAA\",标志头字符串包含的字长数位4，,一子帧总长度为：484，这种格式的子帧将重复5000次" << endl;
	cout << "请对应的cmd指令为:" << endl;
	cout << "{\"bytecount\":\"12\",\"standardstring\":\"AAA AAA AAA AAA\",\"headercount\":\"4\",\"sub_framelength\",\"484\",\"times\":\"5000\"}";
	cout << "重复 输入按照以上格式指定的指令，直到文件规格被叙述完毕" << endl;
	cout << "结束规格时请输入end_register" << endl;
	cout << "注1：完整性约束1，所有子帧的（包括已经输入的首子帧的子帧）长度应当相同，若不满足，FRED规格定义将失败" << endl;
	cout << "注2：以上流程为FRED规格的注册，注册成功后，会在指定目录下生成特定的配置文件，相同卫星的FRED规格无需重复注册，可采用读取配置文件自动建立FRED规格" << endl;

	if (GetInput(bufferlength, singlecount))
	{
		IAddLog("规格建立完毕...！");
		cout << "规格建立完毕！" << endl;
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
		IAddLog("%s文件不存在，FredAPI初始化 失败", _fredpath.c_str());
		return File_Invalid_Existence_Fred;
	}

	char reader;
	while (ifs.peek() != EOF)
	{
		ifs >> reader;  //一开始只能一个一个char的移动

		///如果找到了换行符 ，提示错误
		if (reader == '\x0d')    ///ANSI格式下换行符，本库暂不考虑其他FRED文件不会采用的格式
		{
			ifs >> reader;
			if (reader == '\x0a')
			{
				//报错
				IAddLog("Fred文件中遇到了换行符，文件可能已破坏，位置为%d", ifs.cur - 2);
				IAddLog("程序将退出");
				return Fred_Unexpected_End;
			}
			else
				ifs.seekg(-1, ifs.cur);
		}

		//如果找到了整个帧的开头，进行进一步判断
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

	IAddLog("并没有在Fred文件:%s中找到规定的文件头信息，文件头信息为：");
	return Fred_Cannot_Find_Header;
}

ErrorType FREDStructure::ILocate(string _filepath)
{
	ifstream ifs;
	ifs.open(_filepath.c_str(), ios::in | ios::_Nocreate | ios::binary);

	if (!ifs.is_open())
	{
		cout << "FRED文件不存在，请进行检查" << endl;
		throw(exception("FRED文件不存在，请进行检查"));
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
		IAddLog("Fred文件%s不存在或已损坏", _fredpath.c_str());
		printf("Fred文件%s不存在或已损坏", _fredpath.c_str());
		return File_Invalid_Existence_Fred;
	}

	_endPos = _startPos - 1;
	//ifs.seekg(_endPos, ios::beg);

	int sizeofit = _description.size();
	for (int i = 0; i < sizeofit; i++)
	{
		cout << "Scan进度:" << i + 1 << "/" << sizeofit << "\r";
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
	///切割Fred信息
	vector<unsigned char> heiinfo;
	safestring::TrimBytes(_infobyte, heiinfo, _length, 12);
	_factory.IEquipment(heiinfo);
	return true;
}

///Make sure imgcontent's memory has been located
bool FredEnvelop::IFormImage(unsigned int* imgcontent)
{
	///根据single count 判断如何读取数据
	double _readstep = static_cast<double>(singlebyte) / 8.0;
	int realstep = static_cast<int>(_readstep + 0.5);

	if (imgcontent == nullptr)
	{
		cout << "卫星原始数据导入错误，请检查" << endl;
		IAddLog("卫星原始数据导入错误，请检查");
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
		///然后判断是从低位还是从高位开始读取，以及先后顺序
		for (unsigned int j = _startindex; j < (_startindex + realstep); j++)
		{
			int moveleftsecond = tempint + 8 * (_startindex + realstep -j);// ///第二次左移的位数
			moveleftsecond = (moveleftsecond <= 0) ? 0 : moveleftsecond;
			int moverightfirst; ///第一次右移的位数
			int moveleftfirst;// 第一次左移的位数


			if (_startnum != _startindex&& j == _startindex)
			{
				moveleftfirst = 4;
				moverightfirst = 4;
			}
			///这个判定条件有错误
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
			cout << "数据内容大于接口定义，Error Code: Fatal0E301" << endl;
			IAddLog("数据内容大于接口定义，Error Code: Fatal0E301");
			return false;
		}
	}

	return true;
}

ErrorType EnvelopeAPI::IScan(ifstream& ifs, int id)
{
	EnvelopeDescription hei(_description[id]);

	///判断是否需要补1
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
				IAddLog("错误！在FRED文件的指定位置%d没有发现文件头%s", ifs.cur, hei._dataformat._headerdescprtion.c_str());
				printf("错误！在FRED文件的指定位置%d没有发现文件头%s", ifs.cur, hei._dataformat._headerdescprtion.c_str());
				delete[] tempreader;
				return Fred_Cannot_Find_Header;
			}
		}
		else
		{
			if (!safestring::compare(tempreader, hei._dataformat._header, readheaderlength,true))
			{
				ifs.seekg(-(headerlength + bodylength), ios::cur);
				IAddLog("错误！在FRED文件的指定位置%d没有发现文件头%s", ifs.cur, hei._dataformat._headerdescprtion.c_str());
				printf("错误！在FRED文件的指定位置%d没有发现文件头%s", ifs.cur, hei._dataformat._headerdescprtion.c_str());
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
	///同样，这个地方辅助数据定义的配置文件路径被固定了 日后需要修改
	AuxFactory auxfct("auxBinaryconfig");

	if (auxprintnum >= _auxindex.size())
	{
		cout << "辅助数据已经全部输出完毕" << endl;
		IAddLog("辅助数据分离完毕");
	}

	cout << _auxindex[auxprintnum] << endl;
	//cout << (int)(_content[_auxindex[auxprintnum]]._dataformat._infobyte[100]) << endl;
	Envelop env = _content[_auxindex[auxprintnum]]._dataformat;
	env.IFormAux(auxfct);
	auxprintnum++;
	return auxfct;
}