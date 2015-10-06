#include "stdafx.h"

#include "AGLLIBHELPER.h"





AGLLIBHELPER::~AGLLIBHELPER(void)
{
}

void AGLLIBHELPER::DebugSingle(double x,double y)
{
	_P_SinglePoint_X = x;
	_P_SinglePoint_Y = y;
	_P_SingleCheckActivated = true;
}

void AGLLIBHELPER::CloseDebugSingle()
{
	_P_SingleCheckActivated = false;
}

bool AGLLIBHELPER::_f_GetLineType()
{
	//_theCurve =(_theCurve.blinereport.yfromxrme<((_theCurve.yfromxbarfloaterhormannrpt.rmserror<_theCurve.yfromxsplinerpt.rmserror)?_theCurve.yfromxbarfloaterhormannrpt.rmserror:_theCurve.yfromxsplinerpt.rmserror))?
	//_theCurve.blinereport.yfromxrme: ((_theCurve.yfromxbarfloaterhormannrpt.rmserror<_theCurve.yfromxsplinerpt.rmserror)?_theCurve.yfromxbarfloaterhormannrpt.rmserror:_theCurve.yfromxsplinerpt.rmserror);

	double rmxl1,rmxl2,rmxl3;
	rmxl1 = fabs(_theCurve.blinereport.xfromyrme); rmxl2 = fabs(_theCurve.xfromybarfloaterhormannrpt.rmserror);
	rmxl3 = fabs(_theCurve.xfromysplinerpt.rmserror);
	double rmxr1,rmxr2,rmxr3;
	rmxr1 = fabs(_theCurve.blinereport.yfromxrme); rmxr2 = fabs(_theCurve.yfromxbarfloaterhormannrpt.rmserror);
	rmxr3 = fabs(_theCurve.yfromxsplinerpt.rmserror);
	if(rmxl1<=rmxl2&&rmxl1<=rmxl3)
	{
		_theCurve.typeforlinexfromy = LINE_POLYLINE;
		_theCurve.generalreport.xfromymaxe = _theCurve.blinereport.xfromymaxe;
		_theCurve.generalreport.xfromymine = _theCurve.blinereport.xfromymine;
		_theCurve.generalreport.xfromyrme = _theCurve.blinereport.xfromyrme;
	}

	else if(rmxl2<=rmxl1&&rmxl2<=rmxl3)
	{
		_theCurve.typeforlinexfromy = LINE_RATIONALCURVE;
		_theCurve.generalreport.xfromymaxe = _theCurve.xfromybarfloaterhormannrpt.maxerror;
		_theCurve.generalreport.xfromymaxe = _theCurve.xfromybarfloaterhormannrpt.avgerror;
		_theCurve.generalreport.xfromyrme = _theCurve.xfromybarfloaterhormannrpt.rmserror;
	}

	else if(rmxl3<=rmxl2&&rmxl3<=rmxl1)
	{
		_theCurve.typeforlinexfromy = LINE_SPLINE;
		_theCurve.generalreport.xfromymaxe = _theCurve.xfromysplinerpt.maxerror;
		_theCurve.generalreport.xfromymine = _theCurve.xfromysplinerpt.avgerror;
		_theCurve.generalreport.xfromyrme = _theCurve.xfromysplinerpt.rmserror;
	}

	if(rmxr1<=rmxr2&&rmxr1<=rmxr3)
	{
		_theCurve.typeforlineyfromx = LINE_POLYLINE;
		_theCurve.generalreport.yfromxmaxe = _theCurve.blinereport.yfromxmaxe;
		_theCurve.generalreport.yfromxmine = _theCurve.blinereport.yfromxmine;
		_theCurve.generalreport.yfromxrme = _theCurve.blinereport.yfromxrme;
	}

	else if(rmxr2<=rmxr1&&rmxr2<=rmxr3)
	{
		_theCurve.typeforlineyfromx = LINE_RATIONALCURVE;
		_theCurve.generalreport.yfromxmaxe = _theCurve.yfromxbarfloaterhormannrpt.maxerror;
		_theCurve.generalreport.yfromxmine = _theCurve.yfromxbarfloaterhormannrpt.avgerror;
		_theCurve.generalreport.yfromxrme = _theCurve.yfromxbarfloaterhormannrpt.rmserror;
	}

	else if(rmxr3<=rmxr2&&rmxr3<=rmxr1)
	{
		_theCurve.typeforlineyfromx = LINE_SPLINE;
		_theCurve.generalreport.yfromxmaxe= _theCurve.yfromxsplinerpt.maxerror;
		_theCurve.generalreport.yfromxmine = _theCurve.yfromxsplinerpt.avgerror;
		_theCurve.generalreport.yfromxrme = _theCurve.yfromxsplinerpt.rmserror;
	}

	return true;
}


bool AGLLIBHELPER::_f_GetPOW()
{
	///��ȡ
	try
	{
		polynomialbar2pow(_theCurve.xfromybarycentric,_theCurve.xfromybarycentricPow);
		polynomialbar2pow(_theCurve.yfromxbarycentric,_theCurve.yfromxbarycentricPow);

		polynomialbar2pow(_theCurve.xfromybarfloaterhormann,_theCurve.xfromybarycentricHormannPow);
		polynomialbar2pow(_theCurve.yfromxbarfloaterhormann,_theCurve.yfromxbarycentricHormannPow);


		///��չ��ʽ����ȷ�Խ��в���,����Ҫʹ��������ʽ


	}
	catch(string e)
	{
	}
	return true;
}

bool AGLLIBHELPER::_P_F_GetNormalLine(barycentricinterpolant p,double x,DirectLine& dline)
{
	real_1d_array curvearray;
	polynomialbar2pow(p,curvearray);
	int nums = curvearray.length();

	///��÷���
	double a,b,c; double tempvalue;
	tempvalue = 0;
	for(int i=0;i<nums;i++)
	{
		tempvalue += i * curvearray[i] * pow(x,i-1);
	}

	double y = barycentriccalc(p,x);
	dline = DirectLine(1,tempvalue,-x - y*tempvalue);				
	///���߻�ȡ���

	return false;
}

///�ռ����߷��ߵ�һ��ʽ
///�纯������ʾ
int AGLLIBHELPER::GetNormalVectorFromX(double x,DirectLine& normalline)
{
	try
	{
		switch (_theCurve.typeforlineyfromx)
		{
			//����ʽģʽ
		case LINE_POLYLINE:
			{
				_P_F_GetNormalLine(_theCurve.yfromxbarycentric,x,normalline);
			}
			break; 
		case LINE_RATIONALCURVE:
			{
				_P_F_GetNormalLine(_theCurve.yfromxbarfloaterhormann,x,normalline);
			}
			break;

		case LINE_SPLINE:
			{
				if(fabs(_theCurve.yfromxbarfloaterhormannrpt.rmserror)<fabs(_theCurve.blinereport.yfromxrme))
					_P_F_GetNormalLine(_theCurve.yfromxbarfloaterhormann,x,normalline);
				else
					_P_F_GetNormalLine(_theCurve.yfromxbarycentric,x,normalline);
			}
			break;
		default:
			break;
		}
	}
	catch(string e)
	{
		cout<<e.c_str()<<endl;
		return 0;
	}
	return 1;
}

///GetNORMALVECTORFROMY��GetNormalVectorFromX���õ���ֱ�߶�Ӧ������ϵ���෴��
int AGLLIBHELPER::GetNormalVectorFromY(double y,DirectLine& normalline)
{
	try
	{
		switch (_theCurve.typeforlineyfromx)
		{
			//����ʽģʽ
		case LINE_POLYLINE:
			{
				_P_F_GetNormalLine(_theCurve.xfromybarycentric,y,normalline);
			}
			break; 
		case LINE_RATIONALCURVE:
			{
				_P_F_GetNormalLine(_theCurve.xfromybarfloaterhormann,y,normalline);
			}
			break;

		case LINE_SPLINE:
			{
				if(fabs(_theCurve.xfromybarfloaterhormannrpt.rmserror)<fabs(_theCurve.blinereport.xfromyrme))
					_P_F_GetNormalLine(_theCurve.xfromybarfloaterhormann,y,normalline);
				else
					_P_F_GetNormalLine(_theCurve.xfromybarycentric,y,normalline);
			}
			break;
		default:
			break;
		}
	}
	catch(string e)
	{
		cout<<e.c_str()<<endl;
		return 0;
	}
}



bool AGLLIBHELPER::_P_F_CheckSpread(barycentricinterpolant p ,real_1d_array a,double interpolationvalue)
{
	try{
		///��ȡϵ��
		int paranum = a.length();
		double powresult =0;
		for(int i=0;i<paranum;i++)
		{
			powresult += pow(interpolationvalue,i) * a[i];
		}

		double baryresult  = 0;
		baryresult = barycentriccalc(p,interpolationvalue);
		if(fabs(baryresult- powresult)>0.001*fabs(interpolationvalue))
		{
			string ss = "���󣡺������ʽ��չ��ʽ������������ȴ���ǧ��֮һ������";
			throw ss;
		}
	}

	catch(string e)
	{
		cout<<e.c_str()<<endl;
	}

	return true;
}

double AGLLIBHELPER::GetYFromXValue(double x)
{
	double y;
	///ʹ��try��catch�������쳣
	try
	{
		switch (_theCurve.typeforlineyfromx)
		{
		case LINE_POLYLINE:
			{
				y = barycentriccalc(_theCurve.yfromxbarycentric,x);
			}
			break;
		case LINE_RATIONALCURVE:
			{
				y = barycentriccalc(_theCurve.yfromxbarfloaterhormann,x);
			}
			break;
		case LINE_SPLINE:
			{
				y = spline1dcalc(_theCurve.yfromxspline,x); 
			}
		default:
			{
				y  = -9999;
				///���쳣�ر�����֤Ч�ʵ�ͬʱ�����ɴ�ɱ��ǧ�����Ź�һ���쳣
				string ss = "���񵽴���ֱ��δ����ϵ�����±��ڲ�ȡֵ";
				throw ss;
			}
			break;
		}	
	}
	catch(string e)
	{
		cout<<e.c_str()<<endl;
	}
}


double AGLLIBHELPER::GetXFromYValue(double y)
{
	double x;
	///ʹ��try��catch�������쳣
	try
	{
		switch (_theCurve.typeforlinexfromy)
		{
		case LINE_POLYLINE:
			{
				x = barycentriccalc(_theCurve.xfromybarycentric,y);
			}
			break;
		case LINE_RATIONALCURVE:
			{
				x = barycentriccalc(_theCurve.xfromybarfloaterhormann,y);
			}
			break;
		case LINE_SPLINE:
			{
				x = spline1dcalc(_theCurve.xfromyspline,y); 
			}
		default:
			{
				x = -9999;
				string  ss = "���񵽴���ֱ��δ����ϵ�����±��ڲ�ȡֵ";
				throw ss;
			}
			break;
		}	
	}
	catch(string e)
	{
		cout<<e.c_str()<<endl;
	}

	return x;
}

bool AGLLIBHELPER::InterSectWithDLine(const DirectLine& lineofnv,double startx,bool bdecide ,double* intersectx,double* intersecty)
{
	if(_theCurve.generalreport.xfromyrme<=_theCurve.generalreport.yfromxrme)
	{

				///���ֱ���ϵ�����
				ImgS p1,p2,p3,p4;
				p1.x = 0; p1.x = lineofnv.GetY(0);
				p2.x = _imgWidth; p2.y = lineofnv.GetY(_imgWidth);

				//��ȡ��������ϵ�����
				p3.y = 0; p3.x = GetXFromYValue(0);
				p4.y = _imgHeight; p4.x = GetXFromYValue(_imgHeight);
				///��ȡ��������ϵ�����

				///Ϊ�˱��ַ�װ�ԣ�ʹ�ñ����Լ���װ���ཻ���
				ImgS imgintersection;
				_P_F_GetIntersectionPoint(p1,p2,p3,p4,imgintersection);
				if(bdecide)
				{
					_minxvalue = imgintersection.x;
					_minyvalue = imgintersection.y;
				}
				else
				{
					_maxyvalue = imgintersection.y;
					_maxxvalue = imgintersection.x;
				}

				if(intersectx!=NULL)
				{
					*intersectx = imgintersection.x;
				}
				
				if(intersecty!=NULL)
				{
					*intersecty = imgintersection.y;
				}
	}

	return false;
}


///��ʼ��ģ�飬ֻ��Ϊ�˼���������Ƿ��������
///��ʽģ����Ҫ����1. ˫�����ϣ�X-��Y����Ϻ�Y->X����ϡ�
///2. ���õķ�װ�����������ߣ�����ʽ��ϣ��������߱���ֵ���ַ�����ѡ�񾫶���ߵķ���
///����Ͻ���洢�����Ա����THECURVE��ȥ����Ҫע����ǣ������Ա����Ӧ����˽�е�
///��֤����ֻ���Ⱪ¶��������GETX and GETY�� ����Դ�����е�Curve�����ƣ�
///����ֻ��Ϊ�����Լ����˽⾫�����
///����ά���ԣ���
void AGLLIBHELPER::InitFromImgPoints(vector<ImgS> imgpoints1,int sequence,bool lefts)
{
	string xinitial = "[";
	string yinitial = "[";
	for(int i=0;i<imgpoints1.size();i++)
	{
		char inputvalue[512];

		if(i==(imgpoints1.size() -1))
		{
			sprintf_s(inputvalue,512,"%lf",imgpoints1[i].x);
			xinitial+=string(inputvalue);
			sprintf_s(inputvalue,512,"%lf",imgpoints1[i].y);
			yinitial += string(inputvalue);
		}
		else
		{
			sprintf_s(inputvalue,512,"%lf,",imgpoints1[i].x);
			xinitial+=string(inputvalue);
			sprintf_s(inputvalue,512,"%lf,",imgpoints1[i].y);
			yinitial += string(inputvalue);
		}
	}
	xinitial += "]";
	yinitial += "]";

	real_1d_array x = xinitial.c_str();
	real_1d_array y = yinitial.c_str();

	///�������ǻ����ҵ�y�������С�����ֵ���Է����ҵ�����Ӱ�񲿷ֵı߽�
	int p,q;
	_startyvalue = 999999; _endyvalue = -99999;
	for(p =0;p<y.length();p++)
	{
		if(y[p]>_endyvalue) /*_maxvalue = y[p];*/ {_endyvalue = y[p]; _endxvalue = x[p];}
		if(y[p]<_startyvalue) /*_minvalue = y[p];*/ {_startyvalue = y[p]; _startxvalue = x[p];}
	}
	_valuetype = 0;
	vector<CurveInfo> curves;
	/////ʹ���������ʽ���
	CurveInfo polycurveinfo;
	//barycentricinterpolant p;
	polynomialfitreport polyreport;
	polynomialfitreport polyreportright;
	int infos;
	//polynomialbuild(x, y, _theCurve.yfromxbarycentric);
	//polynomialbuild(y, x, _theCurve.xfromybarycentric);
	polynomialfit(x,y,5,infos,_theCurve.yfromxbarycentric,polyreport);
	polynomialfit(y,x,5,infos,_theCurve.xfromybarycentric,polyreportright);
	///��ö���ʽ��ϵ������������ߵ�һ����ʽ��
	real_1d_array apara;
	polynomialbar2pow(_theCurve.yfromxbarycentric,apara);
	int nnum = apara.length();

	polycurveinfo.Thecurve.polypara = new double[nnum];
	polycurveinfo.Thecurve.type = nnum;
	polycurveinfo.Thecurve.LINETYPE = LINE_POLYLINE;

	for(int i=0;i<nnum;i++)
	{
		polycurveinfo.Thecurve.polypara[i] = apara[i];	//
	}

	/////��ö���ʽ��Ϻ�Ĳв�
	/////�����̫��Ҫ����ʱ�䣬��Ϊ����ʽ���û���ṩһ���������С������������Ž��м���
	//double maxerror = 0;
	//double minerror = 9999;
	//double rmerror = 0;
	//int pointnum = imgpoints1.size();
	//for(int i=2;i<pointnum-5;i++)
	//{
	//	cout<<"From x to y,calculated number:"<<i<<" of total number:"<<pointnum<<"\r";
	//	double v = barycentriccalc(_theCurve.yfromxbarycentric, imgpoints1[i].x);
	//	double errorvalue = v - imgpoints1[i].y;
	//	if(fabs(errorvalue)>fabs(maxerror)) maxerror = errorvalue;
	//	if(fabs(errorvalue)<fabs(minerror)) minerror = errorvalue;
	//	rmerror += (errorvalue*errorvalue)/(pointnum-1);
	//	CurveError er;
	//	er.errorvalue = errorvalue;
	//	er.xpoint = imgpoints1[i].x;
	//	polycurveinfo.TheError.push_back(er);
	//}

	//rmerror = sqrt(rmerror);
	//_theCurve.blinereport.yfromxmaxe =  maxerror;
	//_theCurve.blinereport.yfromxmine = minerror;
	//_theCurve.blinereport.yfromxrme = rmerror ;

	//maxerror = 0;
	//minerror = 99999;
	//rmerror = 0;
	//for(int i=2;i<pointnum-5;i++)
	//{
	//	cout<<"From y to x,calculated number:"<<i<<" of total number:"<<pointnum<<"\r";
	//	double v = barycentriccalc(_theCurve.xfromybarycentric, imgpoints1[i].y);
	//	double errorvalue = v - imgpoints1[i].x;
	//	if(fabs(errorvalue)>fabs(maxerror)) maxerror = errorvalue;
	//	if(fabs(errorvalue)<fabs(minerror)) minerror = errorvalue;
	//	rmerror += (errorvalue*errorvalue)/(pointnum-1);
	//}

	/////��Ԥ��������������������ں��ߺܽӽ���ֱ����ƽ������ɨ����ʱ�����������Ч�����ܺã�ͨ��������ȡ��1
	//rmerror = sqrt(rmerror);

	///����ͳ�ƹ��̱�ȡ����̫ռ��ʱ��

	_theCurve.blinereport.yfromxmine = polyreport.avgerror;
	_theCurve.blinereport.yfromxmaxe = polyreport.maxerror;
	_theCurve.blinereport.yfromxrme = polyreport.rmserror;

	_theCurve.blinereport.xfromymaxe = polyreportright.maxerror;
	_theCurve.blinereport.xfromymine = polyreportright.avgerror;
	_theCurve.blinereport.xfromyrme = polyreportright.rmserror;

	curves.push_back(polycurveinfo);


	///��������
	int info;
	CurveInfo rationalcurve;
	//barycentricinterpolant p1;
	//barycentricfitreport rep;
	barycentricfitfloaterhormann(x,y,x.length(),5,info,_theCurve.xfromybarfloaterhormann,_theCurve.xfromybarfloaterhormannrpt);
	barycentricfitfloaterhormann(y,x,y.length(),5,info,_theCurve.yfromxbarfloaterhormann,_theCurve.yfromxbarfloaterhormannrpt);

	real_1d_array apara1;
	//polynomialbar2pow(p1,apara1);
	int nnum1 = apara1.length();

	rationalcurve.Thecurve.polypara = new double[nnum1];
	rationalcurve.Thecurve.type = nnum1;
	rationalcurve.Thecurve.LINETYPE = LINE_RATIONALCURVE;

	for(int i=0;i<nnum1;i++)
	{
		rationalcurve.Thecurve.polypara[i] = apara1[i];
	}

	///��ö���ʽ��Ϻ�Ĳв�,��report���¾�û��Ҫ��
	//for(int i=0;i<imgpoints1.size();i++)
	//{
	//	double v = barycentriccalc(p1, imgpoints1[i].x);
	//	double errorvalue = v - imgpoints1[i].y;
	//	CurveError er;
	//	er.errorvalue = errorvalue;
	//	er.xpoint = imgpoints1[i].x;
	//	rationalcurve.TheError.push_back(er);
	//}

	//curves.push_back(rationalcurve);

	///�ض�����Ժ��ߵ��������������ר�Ŷ�˫���߽�����ϣ���
	///˫����,���õ��Ǵ����ض�Լ������С�������
	CurveInfo Hyperbola;


	///ʹ�ã����Σ���������
	CurveInfo cubicspline;
	/////������ģ���д������п��ܲ��������㣬��˴˴����ûع�ģ��
	//spline1dinterpolant p2;
	//spline1dfitreport rep2;
	ae_int_t info2;
	ae_int_t Min2 = 100;
	double rho = 0.0;
	spline1dfitpenalized(x,y,Min2,5.0,info2,_theCurve.xfromyspline,_theCurve.xfromysplinerpt);
	spline1dfitpenalized(y,x,Min2,5.0,info2,_theCurve.yfromxspline,_theCurve.yfromxsplinerpt);



	cubicspline.Thecurve.polypara = new double[5];
	cubicspline.Thecurve.type = 5;
	cubicspline.Thecurve.LINETYPE = LINE_SPLINE;

	cubicspline.Thecurve.polypara[0] = _theCurve.xfromysplinerpt.avgerror;
	cubicspline.Thecurve.polypara[1] = _theCurve.xfromysplinerpt.maxerror;
	cubicspline.Thecurve.polypara[2] = _theCurve.xfromysplinerpt.rmserror;

	///�ڴ���û�б�Ҫ������
	//for(int i=0;i<imgpoints1.size();i++)
	//{
	//	double v = spline1dcalc(p2, imgpoints1[i].x);
	//	double errorvalue = v - imgpoints1[i].y;
	//	CurveError er;
	//	er.errorvalue = errorvalue;
	//	er.xpoint = imgpoints1[i].x;
	//	cubicspline.TheError.push_back(er);
	//}

	double test1,test2,test3;
	if(_P_SingleCheckActivated)
	{
		//	test1 = barycentriccalc(p,_P_SinglePoint_X);
		//	test2 = barycentriccalc(p1,_P_SinglePoint_X);
		//	test3 = spline1dcalc(p2,_P_SinglePoint_X);
	}



	///��������ϵĵ㡢�������������������µ��������ϵ������Ϻ�Ĳв������xml�ļ���
	///���գ�����ʵ��ֻ�����ϵĵ㣬����������ϵĲ����϶����������һ��Ҫ��׼�������������Ȳ��ż�д
	XMLHelper xmlhelper;
	char sprintfpath[512];
	sprintf_s(sprintfpath,512,"E:\\ImageProcess_Last\\datas\\TestOfEpiloarEssay\\EpilorLines\\curveinfo%dleft%d.xml",sequence,(int)lefts);


	try{
		///���������ó���õ���Ϸ���
		_f_GetLineType();
		if(fabs(_theCurve.generalreport.yfromxrme)>0.1||fabs(_theCurve.generalreport.xfromyrme)>0.1)
		{
			///���������ͨ��
			string ss = "���󣡣�������Ͼ��Ȳ��ϸ񣬳���0.1������";
			throw  ss;
		}
	}
	catch(string e)
	{
		cout<<e.c_str()<<endl;
	}



	//xmlhelper.OutPutLineInformation(sprintfpath,imgpoints1,curves);
}

void AGLLIBHELPER::_P_F_GetIntersectionPoint(const ImgS& p1,const ImgS& p2,const ImgS& p3,const ImgS& p4, ImgS& intersectionPoint)
{
	//ImgS intersectionPoint;
	intersectionPoint = p1;
	double a1,b1,c1;
	double a2,b2,c2;
	double det_inv;
	double matrix;
	double m1,m2;
	if((p2.x - p1.x)!=0)
		m1 = (p2.y - p1.y)/(p2.x - p1.x);
	else 
		m1 = (double)1e+20;

	if((p4.x - p3.x)!=0)
		m2 = (p4.y - p3.y)/(p4.x - p3.x);
	else 
		m2 = (double)1e+20;

	a1 = m1;
	a2 = m2;
	b1 = -1;
	b2 = -1;
	c1 = (p1.y - m1*p1.x);
	c2 = (p3.y - m2*p3.x);

	det_inv = 1/(a1*b2 - a2*b1);

	intersectionPoint.x = ((b1*c2 - b2*c1)*det_inv);
	intersectionPoint.y = ((a2*c1 - a1*c2)*det_inv);

	return ;
}