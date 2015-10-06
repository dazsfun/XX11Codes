#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <regex>
#include "log.h"
using namespace std;
#include <afx.h>
#include <psapi.h>  
#pragma comment(lib,"psapi.lib")  
using namespace std;

namespace safestring 
{
	void showMemoryInfo(void);
	bool compare(char*, char*, int);
	bool compare(char*, char*, double, bool);
	bool char2byte(char*, int _length, unsigned char*);
	bool char2FredBytes(char* _source, int _length, char* _output, int _singleLength = 3);
	vector<string> _splitstring(const string&, char);
	bool FredBytes2char(const char* source, int count, char* output, int _singleLength = 3);
	bool Upper(string& _input);

	///检查字符串中是否表示16进制字符以外的数据，若有返回false，正常则返回true
	///无论正常与否，_input中的字母都将被转换成UpperCase
	bool hexCheck(string& _input);
	
	ErrorType Is_TimeFormated(string _input, int& year, int& month, int& day, int& hour, int& minute, double& second);

	vector<string> GetRegexResult(const char* str, string _regexstring, const int* _group,int _size);
	vector<string> GetRegexResult(string _input, string _rgexstring);

	bool hexstring2byte(string _input, unsigned char* hei);
	string ReadFile(string _filepathDSDSS);
	void WriteFile(string _filepath, string content);

	string SingleLine(string _input);

	///将一个被补齐的数据转换成unsigned char数据
	bool TrimBytes(char* inputinfo, vector<unsigned char>& output,int length, int singlebyte);

	///基本工具
	int FormatInt(unsigned char* hei, int _length);
	float FormatDouble(unsigned char* hei, int _length);
	short FormatShort(unsigned char* hei, int _length);
	///
	template<typename T> string ValueToStr(T value)
	{
		ostringstream ost;
		ost << value;
		return ost.str();
	}

	string DateTimeString();

	string WholeTimeString();

	bool compare_s(char* left, char* right, int _length);

	template<typename T> void ClearVectorSame(vector<T>& _content)
	{
		vector<T>::iterator itr = _content.begin();
		do
		{
			if (*itr == *(itr + 1))
			{
				itr = _content.erase(itr);
			}
			else
				itr++;//这里迭代器会更新出错
		} while ((itr+1) != _content.end());
	}

	void Array_MinMax(double* datainfi, const int& sizedata, double& maxvalue, double& minvalue);

	bool createDir(string safepath);

	template<typename T> bool CheckVectorSequence(const vector<T>& _content,string hei )
	{
		bool returnpara = true;;
		LogAPI logapi;
		logapi.IAddLog("进行%s类型数组的顺序检查中", hei.c_str());
		 vector<T>::const_iterator itr = _content.begin();
		int ncounter = 0;
		do
		{

			if (*itr >= *(itr + 1))
			{
				logapi.IAddLog("第%d，%d个元素实际顺序出现异常", (*itr).IGetTime(), (*(itr + 1)).IGetTime());
				returnpara = false;
				itr++;
			}
			else
				itr++;//这里迭代器会更新出错
		} while ((itr+1) != _content.end());

		return returnpara;
	}

	template<typename T> void CopyVector(vector<T>& _destination, const vector<T>& _source)
	{
		_destination.insert(_destination.end(), _source.begin(), _source.end());
	}

	int Findneighbour(double number);

	void FindFileSInDir(string _dirpath,vector<string>& _results,string aflix = "*.*");
};