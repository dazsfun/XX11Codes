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
	///获取
	try
	{
		polynomialbar2pow(_theCurve.xfromybarycentric,_theCurve.xfromybarycentricPow);
		polynomialbar2pow(_theCurve.yfromxbarycentric,_theCurve.yfromxbarycentricPow);

		polynomialbar2pow(_theCurve.xfromybarfloaterhormann,_theCurve.xfromybarycentricHormannPow);
		polynomialbar2pow(_theCurve.yfromxbarfloaterhormann,_theCurve.yfromxbarycentricHormannPow);


		///对展开式的正确性进行测试,不需要使用正则表达式


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

	///获得法线
	double a,b,c; double tempvalue;
	tempvalue = 0;
	for(int i=0;i<nums;i++)
	{
		tempvalue += i * curvearray[i] * pow(x,i-1);
	}

	double y = barycentriccalc(p,x);
	dline = DirectLine(1,tempvalue,-x - y*tempvalue);				
	///法线获取完毕

	return false;
}

///空间曲线法线的一般式
///如函数名所示
int AGLLIBHELPER::GetNormalVectorFromX(double x,DirectLine& normalline)
{
	try
	{
		switch (_theCurve.typeforlineyfromx)
		{
			//多项式模式
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

///GetNORMALVECTORFROMY和GetNormalVectorFromX所得到的直线对应的坐标系是相反的
int AGLLIBHELPER::GetNormalVectorFromY(double y,DirectLine& normalline)
{
	try
	{
		switch (_theCurve.typeforlineyfromx)
		{
			//多项式模式
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
		///获取系数
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
			string ss = "错误！函数表达式的展开式计算结果替代精度大于千分之一，请检查";
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
	///使用try，catch来捕获异常
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
				///有异常必报，保证效率的同时，宁可错杀三千，不放过一个异常
				string ss = "捕获到错误，直线未被拟合的情况下被内插取值";
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
	///使用try，catch来捕获异常
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
				string  ss = "捕获到错误，直线未被拟合的情况下被内插取值";
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

				///获得直线上的两点
				ImgS p1,p2,p3,p4;
				p1.x = 0; p1.x = lineofnv.GetY(0);
				p2.x = _imgWidth; p2.y = lineofnv.GetY(_imgWidth);

				//获取拟合曲线上的两点
				p3.y = 0; p3.x = GetXFromYValue(0);
				p4.y = _imgHeight; p4.x = GetXFromYValue(_imgHeight);
				///获取拟合曲线上的两点

				///为了保持封装性，使用本类自己封装的相交监测
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


///初始的模块，只是为了检查核线拟合是否满足规律
///正式模块需要做到1. 双向的拟合，X-》Y的拟合和Y->X的拟合。
///2. 良好的封装，从有理曲线，多项式拟合，样条曲线本类值三种方案种选择精度最高的方案
///将拟合结果存储到类成员变量THECURVE中去，需要注意的是，这个成员变量应当是私有的
///保证本类只向外暴露两个函数GETX and GETY， 这与源代码中的Curve很类似，
///重新只是为了让自己更了解精度情况
///便于维护性：中
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

	///这里我们还想找到y方向的最小，最大值，以方便找到核线影像部分的边界
	int p,q;
	_startyvalue = 999999; _endyvalue = -99999;
	for(p =0;p<y.length();p++)
	{
		if(y[p]>_endyvalue) /*_maxvalue = y[p];*/ {_endyvalue = y[p]; _endxvalue = x[p];}
		if(y[p]<_startyvalue) /*_minvalue = y[p];*/ {_startyvalue = y[p]; _startxvalue = x[p];}
	}
	_valuetype = 0;
	vector<CurveInfo> curves;
	/////使用有理多项式拟合
	CurveInfo polycurveinfo;
	//barycentricinterpolant p;
	polynomialfitreport polyreport;
	polynomialfitreport polyreportright;
	int infos;
	//polynomialbuild(x, y, _theCurve.yfromxbarycentric);
	//polynomialbuild(y, x, _theCurve.xfromybarycentric);
	polynomialfit(x,y,5,infos,_theCurve.yfromxbarycentric,polyreport);
	polynomialfit(y,x,5,infos,_theCurve.xfromybarycentric,polyreportright);
	///获得多项式的系数（即获得曲线的一般表达式）
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

	/////获得多项式拟合后的残差
	/////这项工作太需要计算时间，因为多项式拟合没有提供一个最大误差，最小误差和中误差，这里才进行计算
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

	/////与预测相符，两个方向的拟合在核线很接近垂直或者平行于推扫方向时，不可能拟合效果都很好，通常是两者取其1
	//rmerror = sqrt(rmerror);

	///迭代统计过程被取消，太占用时间

	_theCurve.blinereport.yfromxmine = polyreport.avgerror;
	_theCurve.blinereport.yfromxmaxe = polyreport.maxerror;
	_theCurve.blinereport.yfromxrme = polyreport.rmserror;

	_theCurve.blinereport.xfromymaxe = polyreportright.maxerror;
	_theCurve.blinereport.xfromymine = polyreportright.avgerror;
	_theCurve.blinereport.xfromyrme = polyreportright.rmserror;

	curves.push_back(polycurveinfo);


	///有理曲线
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

	///获得多项式拟合后的残差,有report以下就没必要了
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

	///特定（针对核线的特殊情况，我们专门对双曲线进行拟合），
	///双曲线,采用的是带有特定约束的最小二乘拟合
	CurveInfo Hyperbola;


	///使用（三次）样条曲线
	CurveInfo cubicspline;
	/////由于在模型中存在误差，有可能产生噪声点，因此此处采用回归模型
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

	///在此已没有必要计算了
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



	///将用于拟合的点、各种曲线类型拟合情况下的曲线拟合系数，拟合后的残差输出到xml文件中
	///最终，我们实际只输出拟合的点，至于曲线拟合的参数肯定会输出，但一定要标准化，所以现在先不着急写
	XMLHelper xmlhelper;
	char sprintfpath[512];
	sprintf_s(sprintfpath,512,"E:\\ImageProcess_Last\\datas\\TestOfEpiloarEssay\\EpilorLines\\curveinfo%dleft%d.xml",sequence,(int)lefts);


	try{
		///分析，检查得出最好的拟合方案
		_f_GetLineType();
		if(fabs(_theCurve.generalreport.yfromxrme)>0.1||fabs(_theCurve.generalreport.xfromyrme)>0.1)
		{
			///拟合条件不通过
			string ss = "错误！：曲线拟合精度不合格，超过0.1个像素";
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