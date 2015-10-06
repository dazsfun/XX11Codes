#include "stdafx.h"
#include "safestring.h"

#include <iostream>  



///��ȥ�⻹�ÿ��Ǵ��С��ô
///������������ȫ�ֱ����ķ�ʽ�����̴߳����������ϰ�����Ҫʱ�޸�
union {
	int i;
	unsigned char c[4];
} u;

union
{
	float ivalue;
	unsigned char cvalue[4];
}ud;

union
{
	short ivalue;
	unsigned char cvalue[2];

}signedshortvalue;


#include <io.h>
#include <cstring>
#include <iostream>

using namespace std;

const int MAXLEN = 1024;
//�������Ŀ¼����
unsigned long FILECOUNT = 0;
//��¼�ļ�����


void ListDir(const char* pchData,vector<string>& _container,string _aflix)
{
	_finddata_t   fdata;
	//�����ļ����ҽṹ����
	long   done;
	char tempdir[MAXLEN] = { 0 };
	//����һ����ʱ�ַ����飬�洢Ŀ¼
	strcat_s(tempdir, pchData);
	//�����ַ���
	char formator[1024];
	sprintf_s(formator, "%s", _aflix.c_str());
	strcat_s(tempdir, formator);        //ֻ����ָ����ʽ��
	//�����ַ���
	done = _findfirst(tempdir, &fdata);
	//��ʼ�����ļ�
	if (done != -1)
		//�Ƿ���ҳɹ�
	{
		int ret = 0;
		while (ret != -1)
			//����һ��ѭ��
		{
			if (fdata.attrib != _A_SUBDIR)
				//�ж��ļ�����
			{
				if (strcmp(fdata.name, "...") != 0 && strcmp(fdata.name, "..") != 0 && strcmp(fdata.name, ".") != 0)            //����
				{
					string dir = pchData;
					if (dir.rfind("\\") != dir.size() - 1)
					{
						dir += "\\";
					}

					string filefound = dir + fdata.name;
					_container.push_back(filefound);
					FILECOUNT++;
					//�ۼ��ļ�
				}
			}
			ret = _findnext(done, &fdata);
			//������һ���ļ�
			if (fdata.attrib == _A_SUBDIR && ret != -1)      //�ж��ļ����ԣ������Ŀ¼����ݹ����
			{
				if (strcmp(fdata.name, "...") != 0 && strcmp(fdata.name, "..") != 0 && strcmp(fdata.name, ".") != 0)            //����
				{
					char pdir[MAXLEN] = { 0 };            //�����ַ�����

					strcat_s(pdir, pchData);               //�����ַ���
					strcat_s(pdir, "");
					//�����ַ���

					strcat_s(pdir, fdata.name);            //�����ַ���
					ListDir(pdir,_container,_aflix);
					//�ݹ����
				}
			}
		}
	}
}


/*
int main(void)
{
	while (true)         //���һ��ѭ��
	{
		FILECOUNT = 0;
		char  szFileDir[128] = { 0 }; //����һ���ַ����飬�洢Ŀ¼
		cin >> szFileDir;
		if (strcmp(szFileDir, "e") == 0)  //�˳�ϵͳ
		{
			break;
		}
		ListDir(szFileDir);  //����ListDir��������Ŀ¼
		cout << "����" << FILECOUNT << "���ļ�" << endl;
		//ͳ���ļ�����
	}
	return 0;
}
*/
bool safestring::compare(char* left, char* right, int _length)
{
	char* _left = left;
	char* _right = right;
	for (int i = 0; i < _length; i++)
	{
		if (*_left != *_right || *_left == '\0' || *_right == '\0')
		{
			return false;
		}
		_left++;
		_right++;
	}
	return true;
}

void safestring::showMemoryInfo(void)
{
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
	cout << "�ڴ�ʹ��:" << pmc.WorkingSetSize / 1000 << "K/" << pmc.PeakWorkingSetSize / 1000 << "K + " << pmc.PagefileUsage / 1000 << "K/" << pmc.PeakPagefileUsage / 1000 << "K" << endl;
	cout << pmc.WorkingSetSize << endl;
	int pause_debug;
	cin >> pause_debug;
}

bool safestring::compare_s(char* left, char* right, int _length)
{
	char* _left = left;
	char* _right = right;
	for (int i = 0; i < _length; i++)
	{
		if (*_left != *_right)
		{
			return false;
		}
		_left++;
		_right++;
	}
	return true;
}


bool safestring::compare(char* left, char* right, double length, bool hei)
{
	char* _left = left;
	char* _right = right;
	int _length = length;
	for (int i = 0; i < _length; i++)
	{
		if (*_left != *_right || *_left == '\0' || *_right == '\0')
		{
			return false;
		}
		_left++;
		_right++;
	}

	///���һ����ֻ�Ƚϸ�λ
	unsigned char lastleft = *_left;
	unsigned char lastright = *_right;

	if ((lastleft << 4) != (lastright << 4))
	{
		return false;
	}



	return true;
}

string safestring::DateTimeString()
{
	time_t tNowTime;
	time(&tNowTime);
	tm tLocalTime;
	localtime_s(&tLocalTime, &tNowTime);
	//�õ����ڵ��ַ���
	string strDateTime = ValueToStr(tLocalTime.tm_year + 1900) + "-" +
		ValueToStr(tLocalTime.tm_mon + 1) + "-" +
		ValueToStr(tLocalTime.tm_mday);
	return strDateTime;
}

string safestring::WholeTimeString()
{
	time_t tNowTime;
	time(&tNowTime);
	tm tLocalTime;
	localtime_s(&tLocalTime, &tNowTime);
	//�õ����ڵ��ַ���
	string strDateTime = ValueToStr(tLocalTime.tm_year + 1900) + "-" +
		ValueToStr(tLocalTime.tm_mon + 1) + "-" +
		ValueToStr(tLocalTime.tm_mday) + "-" +
		ValueToStr(tLocalTime.tm_hour) + "-" +
		ValueToStr(tLocalTime.tm_min);
	return strDateTime;
}

void safestring::WriteFile(string filepath, string _content)
{
	ofstream ofs;
	ofs.open(filepath.c_str());
	ofs.write(_content.c_str(), _content.size());
	ofs.close();
}

string safestring::ReadFile(string _inputfile)
{
	ifstream ifs;
	ifs.open(_inputfile.c_str(), ios::in);
	if (!ifs.is_open())
	{
		ifs.close();
		return "Nothing";
	}
	std::stringstream buffer;
	buffer << ifs.rdbuf();
	std::string contents(buffer.str());
	ifs.close();
	return contents;
}

vector<string> safestring::GetRegexResult(const char* str, string _regexstring, const int* _group, int size)
{
	vector<string> _returnpara;
	const::std::tr1::regex pattern(_regexstring.c_str());//,std::basic_regex);
	std::smatch result;

	string _inputstring = str;
	//_inputstring = "<ImageMode>IMGMODE<\/ImageMode>\n<MainFrameNum>555<\/MainFrameNum>\n";
	string temptest = SingleLine(_inputstring);
	bool match = std::regex_match(temptest, result, pattern);
	if (!match)
	{
		LogAPI api;
		api.IAddLog("����ƥ�����Ϊ\n%s\n%s", _regexstring.c_str(), _inputstring.c_str());
		cout << "����ƥ��ʧ��" << endl;
		return _returnpara;
	}


	for (int i = 0; i < size; i++)
	{
		_returnpara.push_back(result[_group[i]].str());
	}
	return _returnpara;
}


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}


std::vector<std::string> safestring::_splitstring(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

bool safestring::FredBytes2char(const char* source, int count, char* output, int _singleLength)
{

	LogAPI api;

	bool addAtailplease = false;

	if (_singleLength * count % 2 != 0)
	{

		api.IAddLog("����ı�־ͷָ��û�в����ֽڣ��϶�Ϊ��֡�ڵ㣬ָ��Ϊ:%s,�����Զ�����֡���β���һ���ֳ���\n����ȷ��ң��������ǰһ���ֳ��Ķκ�", source);

		addAtailplease = true;
		//return false;
	}
	string input = source;

	///����Ƿ�Ϊ16��������
	if (!hexCheck(input))
	{
#ifdef _DEBUG
		api.IAddLog("�����ָ�����16���Ʊ�ʾ��Ϊ��Чָ�ָ��Ϊ:%s,�밴��16����Լ������ָ��", source);
#endif
		return false;
	}


	vector<string> _stringarray = _splitstring(source, ' ');
	if (_stringarray.size() != count)
	{
#ifdef _DEBUG
		api.IAddLog("�����Fred��ʶ�ַ�����Ч����������ָ����������������ָ���ַ���Ϊ��%s,ָ����������Ϊ��%d", source, count);
		return false;
#endif
	}

	string cmdstring = "";
	for (int i = 0; i < count; i++)
	{
		if (_stringarray[i].size() != _singleLength)
		{
#ifdef _DEBUG
			api.IAddLog("�Ӳ���%s���Ȳ�����ָ����һ�ֳ����ȣ�%d(��λ��4λ)��", _stringarray[i].c_str(), _singleLength);
#endif
			return false;
		}
		cmdstring += _stringarray[i];
	}

	////����һ���ֳ���0
	if (addAtailplease)
	{
		for (int i = 0; i < _singleLength; i++)
			cmdstring += "0";
	}

	int iii;
	///����־ͷָ��ת���ɶ����ƽ��
	int charlength = count * _singleLength / 2;
	if (addAtailplease)
		charlength++;
	for (int i = 0; i < charlength; i++)
	{
		unsigned char heiyo = 0;
		unsigned char* hei = &heiyo;
		hexstring2byte(cmdstring.substr(i * 2, 2), hei);
		output[i] = *hei;
	}

	return true;
}

bool safestring::hexCheck(string& _input)
{
	Upper(_input);
	int _sizeofit = _input.size();
	for (int i = 0; i < _sizeofit; i++)
	{
		if ((_input[i] - 'A') <= ('F' - 'A') || (_input[i] - '0') < ('9' - '0'))
			continue;
		return false;
	}
	return true;
}

bool safestring::Upper(string& _input)
{
	int sizeofit = _input.size();
	for (int i = 0; i < sizeofit; i++)
	{
		if ((_input[i] - 'a') < 25 && (_input[i] - 'a') >= 0)
		{
			_input[i] += 'A' - 'a';
		}
	}
	return true;
}

bool safestring::char2FredBytes(char* _source, int _length, char* output, int _singleLength)
{
#ifdef _DEBUG
	LogAPI api;
#endif
	if (_length % _singleLength != 0)
	{
#ifdef _DEBUG
		api.IAddLog("Fred�ļ���ʶ����ʧ�ܣ��涨Ϊ%d���ֽڴ洢������ʶͷ����%d��Ϊ%d�ı���",
			_singleLength * 4, _length, _singleLength);
#endif
		return false;
	}

	unsigned char* _tempuc = new unsigned char[_length];
	char2byte(_source, _length, _tempuc);

	string hei = "";

	for (int i = 0; i < _length; i++)
	{
		char temp[16];
		sprintf_s(temp, sizeof(char)* 16, "%02X", _tempuc[i]);
		hei += temp;
	}

	int actlength = (_singleLength + 1) * (_length / _singleLength);
	for (int i = 0; i < _length / _singleLength; i++)
	{
		for (int j = 0; j < _singleLength; j++)
		{
			output[i * (_singleLength + 1) + j] = hei[i * _singleLength + j];
			;
		}
		output[i * (_singleLength + 1) + _singleLength] = ' ';
	}

	output[actlength - 1] = '\0';
	delete[] _tempuc;
	return true;
}

vector<string> safestring::GetRegexResult(string _input, string _regexstring)
{
	vector<string> _returnpara;

	const std::tr1::regex pattern(_regexstring.c_str());
	string contentinput = _input;
	std::smatch result;
	std::string::const_iterator start = contentinput.begin();
	std::string::const_iterator end = contentinput.end();
	while (std::regex_search(start, end, result, pattern))
	{
		_returnpara.push_back(result[1].str());
		start = result[0].second;
	}
	return _returnpara;
}

bool safestring::hexstring2byte(string _input, unsigned char* hei)
{
	if (!hexCheck(_input))
	{
		return false;
	}

	if (_input.size() != 2)
	{
		//hei = nullptr;
		return false;
	}

	*hei = 0;
	for (int i = 0; i < 2; i++)
		*hei += ((_input[i] - '0') <= ('9' - '0')) ? static_cast<int>((_input[i] - '0') * pow(16.0, 1 - i)) : static_cast<int>((_input[i] - 'A' + 10) * pow(16.0, 1 - i));

	return true;
}

bool safestring::char2byte(char* _source, int _length, unsigned char* _destination)
{
	///������ʹ�þֲ�ָ���������_source���ݣ���ջ�����ݼ����ܱ����ǻ�ʧ
	///char source[] = "ab";
	///����ǿ������ת��ʱ���������  

	for (int i = 0; i < _length; i++)
	{
		_destination[i] = _source[i];
	}
	return true;
}

ErrorType safestring::Is_TimeFormated(string _input, int& year, int& month, int& day, int& hour, int& minute, double& second)
{
	//(\\d{2}|\\d{4})(?:\\-)?([0]{1}\\d{1}|[1]{1}[0-2]{1})(?:\\-)?([0-2]{1}\\d{1}|[3]{1}[0-1]{1})(?:T)?([0-1]{1}\\d{1}|[2]{1}[0-3]{1})(?::)?([0-5]{1}\\d{1})(?::)?([0-5]{1}\\d{1}.?(\\d*)?)
	const std::tr1::regex pattern("(\\d{2}|\\d{4})(?:\\-)?([0]{1}\\d{1}|[1]{1}[0-2]{1})(?:\\-)?([0-2]{1}\\d{1}|[3]{1}[0-1]{1})(?:T)?([0-1]{1}\\d{1}|[2]{1}[0-3]{1})(?::)?([0-5]{1}\\d{1})(?::)?([0-5]{1}\\d{1}.?(\\d*)?)");
	string timeinput = _input;
	std::smatch result;
	bool match = std::regex_match(timeinput, result, pattern);
	if (!match)
	{
		return Format_Invalid_TimeString;
	}

	year = atoi(static_cast<string>(result[1]).c_str());
	month = atoi(static_cast<string>(result[2]).c_str());
	day = atoi(static_cast<string>(result[3]).c_str());
	hour = atoi(static_cast<string>(result[4]).c_str());
	minute = atoi(static_cast<string>(result[5]).c_str());
	second = atof(static_cast<string>(result[6]).c_str());

	return Process_Success;
}

string safestring::SingleLine(string _input)
{
	int sizeofit = _input.length();
	for (int i = 0; i < sizeofit; i++)
	{
		if (_input[i] == 0x0A)
			_input[i] = 'a';
	}
	return _input;
}

///�и�����
bool safestring::TrimBytes(char* inputinfo, vector<unsigned char>& output,int length, int singlebyte)
{
	double rate = singlebyte / 8.0;
	int realstep = static_cast<int>(rate + 0.5);
	int acsinglebyte = static_cast<int>(rate);
	int allcount = rate * length;
	output.reserve(length); output.resize(length);

	for (int i = 0; i < length; i++)
	{
		double _doubleindex = i * rate;
		int _index = static_cast<int>(_doubleindex);

		///�Ȱ����ݶ�����
		int tempint;
		if (rate != static_cast<unsigned int>(rate) && i % 2 == 0)
			tempint = -12;
		else
			tempint = -8;

		unsigned int content = 0;
		for (unsigned int j = _index; j < (_index + realstep); j++)
		{
			int moveleftsecond = tempint + 8 * (_index + realstep - j);
			moveleftsecond = (moveleftsecond <= 0) ? 0 : moveleftsecond;
			int moverightfirst;
			int moveleftfirst;

			if (_doubleindex != _index && j == _index)
			{
				moveleftfirst = moverightfirst = 4;
			}

			else if (j == (_index + realstep - 1) && (i % 2 == 0) && (rate != static_cast<int>(rate)))
			{
				moveleftfirst = 0;
				moverightfirst = 4;
			}
			else
			{
				moveleftfirst = moverightfirst = 0;
			}

			unsigned char heiuchar = inputinfo[j];
			unsigned char heiucharnew = (heiuchar) << moveleftfirst;
			content += static_cast<unsigned int>(heiucharnew >> moverightfirst) << moveleftsecond;
		}

		if (content > pow(2, singlebyte))
		{
			cout << "�������ݴ��ڽӿڶ��壬Error Code: Fatal0E301" << endl;
			return false;
		}

		if (singlebyte % 8 != 0)
		{
			content = content >> 4; //ÿ���ֳ��ĸ�8λ��Ч
		}
		unsigned char heitemp = content; output[i] = heitemp;
	}

	return true;
}

short safestring::FormatShort(unsigned char* _infor, int _length)
{
	memset(signedshortvalue.cvalue, 0, 2);
	short returnNum = 0;
	unsigned short wrongReturn = 0;
	returnNum = -801; wrongReturn = returnNum; returnNum = wrongReturn;
	if (_length != 2)
	{
		cout << "�쳣�� ������short���͵��ֽ���ֻ��Ϊ2" << endl;
		LogAPI api;
		api.IAddLog("�쳣�� ������int���͵��ֽ���ֻ��Ϊ2");
		return returnNum;
	}

	signedshortvalue.cvalue[0] = _infor[1];
	signedshortvalue.cvalue[1] = _infor[0];

	returnNum = 0;

	return signedshortvalue.ivalue;
}

int safestring::FormatInt(unsigned char* _infor,int _length)
{
	memset(u.c, 0, 4);
	int returnNum = 0;
	unsigned int wrongReturn = 0;
	returnNum = -1; wrongReturn = returnNum; returnNum = wrongReturn;

	if (_length > 4)
	{
		cout << "�쳣�� ������int���͵��ֽ������Ϊ4" << endl;
		LogAPI api;
		api.IAddLog("�쳣�� ������int���͵��ֽ������Ϊ4");
		return returnNum;
	}

	returnNum = 0;
	unsigned char formation[4] = {0};// = new unsigned char[4];

	memcpy(formation + 4 - _length, _infor, _length);

	//cout << sizeof(int) << endl;
	for (int i = 0; i < 4; i++)
	{
		u.c[i] = formation[3 - i];
	}

	//delete[] formation;
	return u.i;
}

#include <direct.h>
#include <io.h>
bool safestring::createDir(string safepath)
{
	///Check validation
	if (safepath.rfind("\\") != safepath.size() - 1)
	{
		safepath += "\\";
	}
	int marchingplace = 0;
	do
	{
		marchingplace = safepath.find("\\",marchingplace) + 1;
		if (marchingplace <= 0)
		{
			cout << safepath.c_str() << "Creation Succeed" << endl;
			break;
		}
		string currentDir = safepath.substr(0, marchingplace);
		if (currentDir.find(":") == marchingplace -2)			continue;

		//������·���Ѿ�������,������
		if (_access_s(currentDir.c_str(), 0)==0) continue;
		if (_mkdir(currentDir.c_str()) != 0)
		{
			cout << "·�����а�������ϵͳ������ַ�������" << endl;
			LogAPI api;
			api.IAddLog("Check the path: %s, Creation of this dir has been failed", currentDir.c_str());
			return false;
		}
	} while (true);

	return true;
}


float safestring::FormatDouble(unsigned char* _infor, int _length)
{
	memset(ud.cvalue, 0, 4);
	int returnNum = 0;
	unsigned int wrongReturn = 0;
	returnNum = -1; wrongReturn = returnNum; returnNum = wrongReturn;

	if (_length > 4)
	{
		cout << "�쳣�� ������double���͵��ֽ������Ϊ4" << endl;
		LogAPI api;
		api.IAddLog("�쳣�� ������double���͵��ֽ������Ϊ4");
		return returnNum;
	}

	double returnValid1 = 0;
	unsigned char formation[4] = { 0 };// = new unsigned char[4];



	memcpy(formation + 4 - _length, _infor, _length);

	for (int i = 0; i < 4; i++)
	{
		ud.cvalue[i] = formation[3-i];
	}

	float x = *(float *)&formation;

	//memcpy(&returnValid1, formation, 4);

	//delete[] formation;
	return ud.ivalue;
}

void safestring::FindFileSInDir(string _dirpath, vector<string>& _results,string aflix)
{
	ListDir(_dirpath.c_str(), _results, aflix);
}

void safestring::Array_MinMax(double* datainfi, const int& sizedata, double& maxvalue, double& minvalue)
{
	maxvalue = -INFINITY;
	minvalue = INFINITY;
	for (int i = 0; i < sizedata; i++)
	{
		maxvalue = (datainfi[i] > maxvalue) ? datainfi[i] : maxvalue;
		minvalue = (datainfi[i] < minvalue) ? datainfi[i] : minvalue;
	}
}



int safestring::Findneighbour(double number)
{
	return static_cast<int>(number + 0.5);
}