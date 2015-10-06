#pragma once

#include "ImageModel.h"
#include "ImgPoint.h"
#include "OpenCVHelper.h"
#include "ap.h"
#include "integration.h"
#include "interpolation.h"
#include "safeopencv.h"
#include "safestring.h"
using namespace alglib;


struct SwingCablib
{
	double _start;
	double _end;
	double _avg;
	double _rms;
};

typedef struct Time_Calibstruct
{
	barycentricinterpolant _interkernel;
	double _startvalue;
	double _endvalue;
}TimeCalibStruct;


#pragma  comment(lib,"MatchBundleModel1.lib")
/*
创建这个类型的思路是  一个被检校好的红外（摆扫）模型，其需要走过的路径是
根据辅助数据建立的初始模型 -> 对整个成像刚体校正后的模型（外畸变）->
->对单列摆扫光束的刚体校正(校正摆扫畸变) -> 对单个CCD光束的刚体校正（校正内部畸变）
根据粒度的不同  对应不同粒度的子类
CalibratedModel 对整个成像刚体校正后的模型
一个模型一旦建立 ，它就应该被视为唯一的
每一层成像模型的实例只有一个
因此  所有模型不允许 默认构造
不允许拷贝、复制
只允许从基类对象集成，然后反馈给子类对象，并且不允许越级反馈
*/
///Calibra
class CalibratedModel :
	public JSONOBJECTSerializer
{
public:
	int pownumber;
	//CalibratedModel();
	CalibratedModel(ImageModel* preparedinfo);
	//不允许拷贝 ，不允许重新赋值，只允许传一个指针进来，因为所有上一级的信息在下一级中都是不可以被修改的
	virtual void SetPropertys();
	virtual void SetUTCAPIByTimestring() {}
	void Calibration(vector<ImgPoint> m_points, ImgPoint thepoint);
	vector<TimeCalibStruct> _arryay_inter;
	vector<TimeCalibStruct> _arryay_inter_x;
private:
	double AScan2Cameara0Roll;
	double AScan2Cameara1Pitch;
	double AScan2Cameara2YAW;
	double BCamera2Install0Roll;
	double BCamera2Install1Pitch;
	double BCamera2Install2YAW;
	double CInstall2Body0Roll;
	double CInstall2Body1Pitch;
	double CInstall2Body2YAW;
	double DCalibratedPara0Roll;
	double DCalibratedPara1Pitch;
	double DCalibratedPara2YAW;

	double EDateTimeString;
	double EControlPointErrorX;
	double EControlPointErrorY;


	void SortAndReduct(vector<ImgPoint> hei, vector<ImgPoint>& heinow)
	{
		///对控制点按照CCD探元索引号排列
		std::sort(hei.begin(), hei.end(), myobject2);

		cout << "排序完毕!" << endl;
		///去除其中距离整数索引过远的点
		vector<ImgPoint>::iterator itr = hei.begin();
		do
		{
			if (fabs((*itr).ImgX() - safestring::Findneighbour((*itr).ImgX())) > 0.1)
			{
				itr++;
			}
			else
			{
				heinow.push_back(ImgPoint(*itr, false));
				itr++;//这里迭代器会更新出错
			}
		} while ((itr + 1) != hei.end());

		cout << "筛选完毕!" << endl;
	}

	vector<vector<ImgPoint>> swing_division;
public:
	struct myclass1 {
		bool operator() (ImgPoint pt1, ImgPoint pt2) { return (pt1.ImgY() < pt2.ImgY()); }
	} myobject1;

	struct myclass2 {
		bool operator() (ImgPoint pt1, ImgPoint pt2) { return (pt1.ImgX() < pt2.ImgX()); }
	} myobject2;

	void Divide_Swing( vector<ImgPoint>& heinow,int rooms,int start_line,int end_line,int steps);

	void Link_Swing_Division(int index,vector<ImgPoint>& uppper_line,vector<ImgPoint>& down_line);
	
	void Calibration_Swing(vector<ImgPoint> heinow, ImgPoint base_p,int start_line,int endline,int steps);

	void CalibrateSwing(vector<ImgPoint> heinow, ImgPoint base_p, vector<SwingCablib>& _swingcalib, double& endswing)
	{
		_swingcalib.clear();
		vector<ImgPoint> arraypoint;
		double xmoveavg, ymoveavg;
		double xmoverms, ymoverms;
		xmoveavg = ymoveavg = xmoverms = ymoverms = 0;

		for (int i = 0; i < heinow.size(); i++)
		{
			xmoveavg += (heinow[i].ImgX() - heinow[i]._rotation4) / heinow.size();
			ymoveavg += (heinow[i].ImgY() - heinow[i]._rotation5) / heinow.size();
		}

		for (int i = 0; i < heinow.size(); i++)
		{
			xmoverms += pow(heinow[i].ImgX() - heinow[i]._rotation4 - xmoveavg, 2.0) / (heinow.size() - 1);
			ymoverms += pow(heinow[i].ImgY() - heinow[i]._rotation5 - ymoveavg, 2.0) / (heinow.size() - 1);
		}

		xmoverms = sqrt(xmoverms);
		ymoverms = sqrt(ymoverms);

		for (int i = 0; i < heinow.size(); i++)
		{
			if (fabs(heinow[i].ImgX() - heinow[i]._rotation4) > 12)
				continue;
			if (fabs(heinow[i].ImgY() - heinow[i]._rotation5) > 12)
				continue;
			arraypoint.push_back(ImgPoint(heinow[i], false));
		}

		//arraypoint.insert(arraypoint.begin(), heinow.begin(), heinow.end());
		std::sort(arraypoint.begin(), arraypoint.end(), myobject1);

		ofstream ofs2("c:\\Temp\\match\\摆扫畸变输入.xls");
		for (int i = 0; i < arraypoint.size(); i++)
		{
			ofs2 << arraypoint[i].ImgX() << "\t" << arraypoint[i]._rotation4 << "\t" << arraypoint[i].ImgY() << "\t" << arraypoint[i]._rotation5 << "\t" << arraypoint[i]._rotation4 - arraypoint[i].ImgX() << "\t" <<
				arraypoint[i]._rotation5 - arraypoint[i].ImgY() << endl;
		}


		double ybvalue = arraypoint[0].ImgY();
		double yevalue = arraypoint[arraypoint.size() - 1].ImgY();

		cout << "获取到的观测值，y上限为：" << ybvalue << "y下限位：" << yevalue << endl;
		cout << "X方向畸变平均值：%lf" << xmoveavg << "\t" << "y方向畸变平均值:" << ymoveavg << endl;
		cout << "X方向畸变中误差：%lf" << xmoverms << "\t" << "y方向畸变中误差:" << ymoverms << endl;


		double itbegin = ybvalue;

		int itindex1, itindex2;
		itindex1 = 0;

		//vector<SwingCablib> heiarray;
		do
		{
			cout << itindex1 << endl;
			double yshiftvaluee = 0;
			double itend = itbegin + 5;
			int ncounter = 0;
			//if (itend > ybvalue) break;
			vector<double> shitfarray;
			bool noline = false;
			for (int i = itindex1; i < arraypoint.size(); i++)
			{
				itindex2 = i;
				if (arraypoint[i].ImgY() > itend)
				{
					if (i == itindex1)
					{
						noline = true;
					}
					itindex1 = itindex2;
					break;
				}

				for (int j = 0; j < arraypoint.size(); j++)
				{
					if (j == i) continue;
					if (fabs(arraypoint[j].ImgY() - arraypoint[i].ImgY()) < 10
						/*&&b
						fabs(arraypoint[j].ImgY() - arraypoint[i].ImgY()) > 1&&
						fabs(arraypoint[j].ImgX() - arraypoint[i].ImgX()) > 1&&
						fabs(arraypoint[j].ImgX() - arraypoint[i].ImgX()) > 5*/)
					{
						//对于5个像素的分段来说， 最大的误差是12像素
						double scale = fabs((arraypoint[i].ImgY() - arraypoint[j].ImgY()) / (arraypoint[j]._rotation5 - arraypoint[i]._rotation5));
						if (scale < 1.01 && scale >(1.0 / 1.01))
						{
							ncounter++;
							shitfarray.push_back(fabs((arraypoint[i].ImgY() - arraypoint[j].ImgY()) / (arraypoint[j]._rotation5 - arraypoint[i]._rotation5)));
							yshiftvaluee += fabs((arraypoint[i].ImgY() - arraypoint[j].ImgY()) / (arraypoint[j]._rotation5 - arraypoint[i]._rotation5));
						}
					}

					else if ((arraypoint[j].ImgY() - arraypoint[i].ImgY()) > 10){
						continue;
					}
				}

			}

			double rms_error = 0;

			if (ncounter == 0)
			{
				if (noline)
					;//cout << "区域（" << itbegin << "," << itend << ")" << "内没有匹配点"<<endl;
				else
					cout << "区域（" << itbegin << "," << itend << ")" << "内点不够密集！！" << endl;
				if (_swingcalib.size() == 0)
				{
					SwingCablib heinow;// = _swingcalib[_swingcalib.size() - 1];
					heinow._avg = 1;
					heinow._rms = 0;
					heinow._start = itbegin;
					heinow._end = itend;
					_swingcalib.push_back(heinow);
					itbegin = itend;
					itindex1 = itindex2;
				}
				else
				{
					SwingCablib heinow = _swingcalib[_swingcalib.size() - 1];
					heinow._start = itbegin;
					heinow._end = itend;
					_swingcalib.push_back(heinow);
					itbegin = itend;
					itindex1 = itindex2;
				}
			}
			else
			{
				yshiftvaluee /= ncounter;
				for (int t = 0; t < shitfarray.size(); t++)
				{
					rms_error += pow(shitfarray[t] - yshiftvaluee, 2.0) / (shitfarray.size() - 1);
				}

				rms_error = sqrt(rms_error);
				SwingCablib heinow;
				heinow._start = itbegin;
				heinow._end = itend;
				heinow._avg = yshiftvaluee;
				heinow._rms = rms_error;

				itbegin = itend;
				itindex1 = itindex2;
				_swingcalib.push_back(heinow);
			}

			///检查区域,index1,index2

		} while (itbegin < yevalue);

		cout << "正确的拟合组数应当是" << static_cast<int>((yevalue - ybvalue) / 5.0) + 1 << endl;
		cout << "实际的组数是:" << _swingcalib.size() << endl;

		ofstream ofs("c:\\Temp\\match\\restresult.xls");
		ofs.precision(20);
		for (int i = 0; i < _swingcalib.size(); i++)
		{
			ofs << _swingcalib[i]._start << "\t" << _swingcalib[i]._end << "\t" << _swingcalib[i]._avg << "\t" << _swingcalib[i]._rms << "\t" << yevalue << endl;
		}
		vector<Point3d> comp;
		for (double i = yevalue; i > ybvalue; i--)
		{
			Point3d comp1;
			comp1.x = (yevalue - i);
			int step = (yevalue - i) / 5;
			double re = 0;
			for (int j = _swingcalib.size() - 1, i = 0; i < step; j--, i++)
			{
				re += 5 * _swingcalib[j]._avg;
			}

			re += ((yevalue - i) - 5 * step) * _swingcalib[_swingcalib.size() - 1 - step]._avg;
			comp1.y = re;
			comp1.z = i;
			comp.push_back(comp1);
		}

		ofs.close();

		ofstream ofscureve("c:\\Temp\\match\\mycurve.xls");
		ofscureve.precision(20);
		for (int i = 0; i < comp.size(); i++)
		{
			ofscureve << comp[i].z << "\t" << comp[i].x << "\t" << comp[i].y << "\t" << comp[i].z << "\t" << comp[i].x - comp[i].y << endl;
		}

		ofscureve.close();
	}



	void CalibrateLookAngle(vector<ImgPoint> hei, ImgPoint base_p)
	{
		CSatOrbit helper;
		//对导入的控制点进行排序，筛选，去除不利于进行指向角计算的点
		vector<ImgPoint> heinow;
		cout << "排序筛选前匹配控制点个数" << hei.size() << endl;
		SortAndReduct(hei, heinow);
		cout << "筛选后用于指向角校正的匹配控制点个数 " << heinow.size() << endl;
		//cout << "时间周期解算控制点个数" << toplines.size() << endl;
		//cout << "摆扫成像距离为" << srange << endl;

		int calculationsize = heinow.size();

		//double Kappa = k1;
		//double T0 = Vt;
		////对每个点 ,列误差方程 
		int pointread = 0;
		vector<Point3d> lookerror;		int inputcounter = 0;
		vector<Point3d> hpointarr;
		ofstream ofs_look_distance("c:\\Temp\\match\\lookangle_find_error_distance.xls");
		for (int i = 0; i < heinow.size(); i += inputcounter)
		{
			///确定索引号 
			int x_index = safestring::Findneighbour(heinow[i].ImgX());
			inputcounter = 0;
			do{
				if (safestring::Findneighbour(heinow[i + inputcounter].ImgX()) != x_index)
					break;
				inputcounter++;
				pointread++;
			} while (true);

			if (inputcounter < 2) {
				cout << "!!" << endl; continue;
			}

			cout << "对线阵CCD第" << x_index << "个探元进行校正" << endl;
			cout << "可用于标定的匹配控制点个数为:" << inputcounter << endl;
			double templ[3];
			base_p.ILightVector(x_index, x_index, templ);
			int ncounter = 0;
			double hei222 = 1.0; hei222 /= templ[2];
			cout << hei222 << endl;
		//	cin >> ncounter;

			double lx[2];			double ly[2];
			lx[0] = lx[1] = hei222 * templ[0];
			ly[0] = ly[1] = hei222 * templ[1];
			cout << "原有的指向角大小为: (" << lx[0] << "," << ly[0] << "," << 1 << ")" << endl;

			Point3d hpoint; hpoint.z = x_index;
			ofs_look_distance.precision(20);
			do
			{

				double* coeffarray = new double[inputcounter * 2 * 2];
				double* staticarray = new double[inputcounter * 2]; // 观测值
				double* tarray = new double[inputcounter * 2 * 2];  // 常量矩阵
				double CTC[4] = { 0 };
				double CTL[4] = { 0 };


				int calc = 0;
				for (int j = i; j < i + inputcounter; j++)
				{

					if (fabs(heinow[j].ImgX() - x_index)> 0.5)
					{
						cout << "异常情况，找不到索引号！" << endl;
						cout << heinow[j].ImgX() << "\t" << x_index << endl;
						cout << "是否终止程序" << endl;
						int infonow = 0;
						cin >> infonow;
						if (infonow == 1)
							abort();
					}

					double rotforpoint[9] = { 0 };
			//		cout << heinow[j]._rotation7 << "\t" << heinow[j]._rotation8 << "\t" << heinow[j]._rotation9 << endl;
				//	cout << heinow[j].X() << "\t" << heinow[j].Y() << "\t" << heinow[j].Z() << endl;
				//	cout << heinow[j].OriX() << "\t" << heinow[j].OriY() << "\t" << heinow[j].OriZ() << endl;

				//	cout << lx[0] << "\t" << ly[0] << endl;
					safeopencv::Eulor2Rot(heinow[j]._rotation7, heinow[j]._rotation8,
						heinow[j]._rotation9, rotforpoint);
					//cout << heiyou << endl;
				//	cout << lx[0] << "\t" << ly[0] << endl;
					double A1 = (rotforpoint[0] * lx[0] + rotforpoint[1] * ly[0] + rotforpoint[2]);
					double B1 = (rotforpoint[3] * lx[0] + rotforpoint[4] * ly[0] + rotforpoint[5]);
					double C1 = (rotforpoint[6] * lx[0] + rotforpoint[7] * ly[0] + rotforpoint[8]);
					staticarray[calc * 2 + 0] = (heinow[j].X() - heinow[j]._rotation4) / (heinow[j].Z() - heinow[j]._rotation6)
						- A1/C1;// (rotforpoint[0] * lx[0] + rotforpoint[1] * ly[0] + rotforpoint[2]) / (rotforpoint[6] * lx[0] + rotforpoint[7] * ly[0] + rotforpoint[8]);

					//cout << (heinow[j].X() - heinow[j]._rotation4) / (heinow[j].Z() - heinow[j]._rotation6) << endl;
					staticarray[calc * 2 + 1] = (heinow[j].Y() - heinow[j]._rotation5) / (heinow[j].Z() - heinow[j]._rotation6)
						- B1/C1;// (rotforpoint[3] * lx[0] + rotforpoint[4] * ly[0] + rotforpoint[5]) / (rotforpoint[6] * lx[0] + rotforpoint[7] * ly[0] + rotforpoint[8]);
					//cout << (heinow[j].Y() - heinow[j]._rotation5) / (heinow[j].Z() - heinow[j]._rotation6) << endl;

				//	cout << staticarray[calc * 2 + 0] << endl;
			//		cout << staticarray[calc * 2 + 1] << endl;
	




					coeffarray[calc * 2 * 2 + 0] = (rotforpoint[0] * C1 - A1 * rotforpoint[6]) / (C1 * C1);
					coeffarray[calc * 2 * 2 + 1] = (rotforpoint[1] * C1 - A1 * rotforpoint[7]) / (C1 * C1);
					coeffarray[calc * 2 * 2 + 2] = (rotforpoint[3] * C1 - B1 * rotforpoint[6]) / (C1 * C1);
					coeffarray[calc * 2 * 2 + 3] = (rotforpoint[4] * C1 - B1 * rotforpoint[7]) / (C1 * C1);
					calc++;
				}

				helper.transpose(coeffarray, tarray, inputcounter * 2, 2);
				helper.mult(tarray, coeffarray, CTC, 2, inputcounter * 2, 2);
				helper.mult(tarray, staticarray, CTL, 2, inputcounter * 2, 1);



				//CTL为2X1
				if (!helper.invers_matrix(CTC, 2))		 					//求逆
				{
					cout << "不可逆" << endl;
					;
				}

				double V[2];

				
				helper.mult(CTC, CTL, V, 2, 2, 1);
				lx[0] += V[0];
				ly[0] += V[1];

				delete[] coeffarray;
				delete[] staticarray;
				delete[] tarray;
				ncounter++;
			} while (ncounter < 10);

			templ[0] = lx[0]; templ[1] = ly[0];
			cout << "第" << x_index << "个探元" << endl;
			cout << "校正后的指向角大小为: (" << templ[0] << "," << templ[1] << ")" << endl;
			cout << "修正比例:" << fabs(lx[1] - lx[0]) / lx[1] << "," << fabs(ly[1] - ly[0]) / ly[1] << endl;
			ofs_look_distance << x_index << "\t" << fabs(lx[1] - lx[0]) / lx[1] << "\t" << fabs(ly[1] - ly[0]) / ly[1] << endl;

			int pause_12;
			//cin >> pause_12;
			hpoint.x = lx[0];
			hpoint.y = ly[0];
			hpointarr.push_back(hpoint);
		}

		ofstream ofslookangle("c:\\Temp\\match\\lookangle_find_error.txt");
		ofslookangle.precision(20);
		for (int i = 0; i < hpointarr.size(); i++)
		{
			ofslookangle << static_cast<int>(hpointarr[i].z) << "\t" << hpointarr[i].x << "\t" << hpointarr[i].y << endl;
		}
		ofslookangle.close();
	}


	void CalibrateTime_(vector<ImgPoint> heinow, ImgPoint base_p, vector<SwingCablib>& _swingcalib, double& endswing)
	{
		_swingcalib.clear();
		vector<ImgPoint> arraypoint;
		double xmoveavg, ymoveavg;
		double xmoverms, ymoverms;
		xmoveavg = ymoveavg = xmoverms = ymoverms = 0;

		for (int i = 0; i < heinow.size(); i++)
		{
			xmoveavg += (heinow[i].ImgX() - heinow[i]._rotation4) / heinow.size();
			ymoveavg += (heinow[i].ImgY() - heinow[i]._rotation5) / heinow.size();
		}

		for (int i = 0; i < heinow.size(); i++)
		{
			xmoverms += pow(heinow[i].ImgX() - heinow[i]._rotation4 - xmoveavg, 2.0) / (heinow.size() - 1);
			ymoverms += pow(heinow[i].ImgY() - heinow[i]._rotation5 - ymoveavg, 2.0) / (heinow.size() - 1);
		}

		xmoverms = sqrt(xmoverms);
		ymoverms = sqrt(ymoverms);

		for (int i = 0; i < heinow.size(); i++)
		{
			if (fabs(heinow[i].ImgX() - heinow[i]._rotation4) > 12)
				continue;
			if (fabs(heinow[i].ImgY() - heinow[i]._rotation5) > 12)
				continue;
			arraypoint.push_back(ImgPoint(heinow[i], false));
		}

		//arraypoint.insert(arraypoint.begin(), heinow.begin(), heinow.end());
		std::sort(arraypoint.begin(), arraypoint.end(), myobject2);

		ofstream ofs2("c:\\Temp\\match\\摆扫畸变输入_x.xls");
		for (int i = 0; i < arraypoint.size(); i++)
		{
			ofs2 << arraypoint[i].ImgX() << "\t" << arraypoint[i]._rotation4 << "\t" << arraypoint[i].ImgY() << "\t" << arraypoint[i]._rotation5 << "\t" << arraypoint[i]._rotation4 - arraypoint[i].ImgX() << "\t" <<
				arraypoint[i]._rotation5 - arraypoint[i].ImgY() << endl;
		}


		double xbvalue = arraypoint[0].ImgX();
		double xevalue = arraypoint[arraypoint.size() - 1].ImgX();

		cout << "获取到的观测值，x上限为：" << xbvalue << "x下限位：" << xevalue << endl;
		cout << "X方向畸变平均值：%lf" << xmoveavg << "\t" << "y方向畸变平均值:" << ymoveavg << endl;
		cout << "X方向畸变中误差：%lf" << xmoverms << "\t" << "y方向畸变中误差:" << ymoverms << endl;


		double itbegin = xbvalue;

		int itindex1, itindex2;
		itindex1 = 0;

		//vector<SwingCablib> heiarray;
		do
		{
			cout << itindex1 << endl;
			double Xshiftvaluee = 0;
			double itend = itbegin + 3;
			int ncounter = 0;
			//if (itend > ybvalue) break;
			vector<double> shitfarray;
			bool noline = false;
			for (int i = itindex1; i < arraypoint.size(); i++)
			{
				itindex2 = i;
				if (arraypoint[i].ImgX() > itend)
				{
					if (i == itindex1)
					{
						noline = true;
					}
					itindex1 = itindex2;
					break;
				}

				for (int j = 0; j < arraypoint.size(); j++)
				{
					if (j == i) continue;
					if (fabs(arraypoint[j].ImgX() - arraypoint[i].ImgX()) < 3
						/*&&b
						fabs(arraypoint[j].ImgY() - arraypoint[i].ImgY()) > 1&&
						fabs(arraypoint[j].ImgX() - arraypoint[i].ImgX()) > 1&&
						fabs(arraypoint[j].ImgX() - arraypoint[i].ImgX()) > 5*/)
					{
						//对于5个像素的分段来说， 最大的误差是12像素
						double scale = fabs((arraypoint[i].ImgX() - arraypoint[j].ImgX()) / (arraypoint[j]._rotation4 - arraypoint[i]._rotation4));
						if (scale < 1.01 && scale >(1.0 / 1.01))
						{

							ncounter++;
							shitfarray.push_back(scale);// (fabs((arraypoint[i].ImgX() - arraypoint[j].ImgX()) / (arraypoint[j]._rotation4 - arraypoint[i]._rotation4)));
							Xshiftvaluee += scale;// fabs((arraypoint[i].ImgX() - arraypoint[j].ImgX()) / (arraypoint[j]._rotation4 - arraypoint[i]._rotation4));
						}
					}

					else if ((arraypoint[j].ImgX() - arraypoint[i].ImgX()) > 10){
						continue;
					}
				}

			}

			double rms_error = 0;

			if (ncounter == 0)
			{
				if (noline)
					;//cout << "区域（" << itbegin << "," << itend << ")" << "内没有匹配点"<<endl;
				else
					cout << "区域（" << itbegin << "," << itend << ")" << "内点不够密集！！" << endl;
				if (_swingcalib.size() == 0)
				{
					SwingCablib heinow;// = _swingcalib[_swingcalib.size() - 1];
					heinow._avg = 1;
					heinow._rms = 0;
					heinow._start = itbegin;
					heinow._end = itend;
					_swingcalib.push_back(heinow);
					itbegin = itend;
					itindex1 = itindex2;
				}
				else
				{
					SwingCablib heinow = _swingcalib[_swingcalib.size() - 1];
					heinow._start = itbegin;
					heinow._end = itend;
					_swingcalib.push_back(heinow);
					itbegin = itend;
					itindex1 = itindex2;
				}
			}
			else
			{
				Xshiftvaluee /= ncounter;
				for (int t = 0; t < shitfarray.size(); t++)
				{
					rms_error += pow(shitfarray[t] - Xshiftvaluee, 2.0) / (shitfarray.size() - 1);
				}

				rms_error = sqrt(rms_error);
				SwingCablib heinow;
				heinow._start = itbegin;
				heinow._end = itend;
				heinow._avg = Xshiftvaluee;
				heinow._rms = rms_error;

				itbegin = itend;
				itindex1 = itindex2;
				_swingcalib.push_back(heinow);
			}

			///检查区域,index1,index2

		} while (itbegin < xevalue);

		cout << "正确的拟合组数应当是" << static_cast<int>((xevalue - xbvalue) / 3.0) + 1 << endl;
		cout << "实际的组数是:" << _swingcalib.size() << endl;

		ofstream ofs("c:\\Temp\\match\\restresult_x.xls");
		ofs.precision(20);
		for (int i = 0; i < _swingcalib.size(); i++)
		{
			ofs << _swingcalib[i]._start << "\t" << _swingcalib[i]._end << "\t" << _swingcalib[i]._avg << "\t" << _swingcalib[i]._rms << "\t" << xevalue << endl;
		}
		vector<Point3d> comp;
		for (double i = xevalue; i > xbvalue; i--)
		{
			Point3d comp1;
			comp1.x = (xevalue - i);
			int step = (xevalue - i) / 3;
			double re = 0;
			for (int j = _swingcalib.size() - 1, i = 0; i < step; j--, i++)
			{
				re += 3 * _swingcalib[j]._avg;
			}

			re += ((xevalue - i) - 3 * step) * _swingcalib[_swingcalib.size() - 1 - step]._avg;
			comp1.y = re;
			comp1.z = i;
			comp.push_back(comp1);
		}

		ofs.close();

		ofstream ofscureve("c:\\Temp\\match\\mycurve_x.xls");
		ofscureve.precision(20);
		for (int i = 0; i < comp.size(); i++)
		{
			ofscureve << comp[i].z << "\t" << comp[i].x << "\t" << comp[i].y << "\t" << comp[i].z << "\t" << comp[i].x - comp[i].y << endl;
		}

		ofscureve.close();
	}


	double CalibrateTime(const vector<Point2d>& arraypoint, bool xfactor = false);

public:
	void IAngleCalibration(vector<TimeCalibStruct>& heire, bool bxarray = false)
	{
		if (!bxarray)
		{
			for (int i = 0; i < _arryay_inter.size(); i++)
				heire.push_back(_arryay_inter[i]);
		}
		else
		{
			for (int i = 0; i < _arryay_inter_x.size(); i++)
				heire.push_back(_arryay_inter_x[i]);
		}
	}

public:
	void FindTime(vector<ImgPoint> controlist, ImgPoint basepoint, double Vx, double k1, double k2, double& Vt, double srange)
	{
		///在控制点中寻找符合要求的
		vector<ImgPoint> toplines;
		vector<ImgPoint> alignlines;
		for (int i = 0; i < controlist.size(); i++)
		{
			if (controlist[i].ImgY() < 60)
			{
				for (int j = 0; j < controlist.size(); j++)
				{
					if ((controlist[j].ImgX() - controlist[i].ImgX())<0.999&&controlist[j].ImgY()>10400)
					{
						toplines.push_back(ImgPoint(controlist[i], false));
						alignlines.push_back(ImgPoint(controlist[j], false));
					}
				}
			}
		}

		cout << "时间周期解算控制点个数" << toplines.size() << endl;
		cout << "摆扫成像距离为" << srange << endl;

		int calculationsize = toplines.size();
		double* coeffarray = new double[calculationsize * 2];
		double* staticarray = new double[calculationsize]; // 观测值
		double* tarray = new double[calculationsize * 2];  // 常量矩阵

		double CTC[4] = { 0 };
		double CTL[2] = { 0 };

		double pt[2] = { 0 };
		double V0 = Vx;
		double K0 = k2;
		double Kappa = k1;
		double T0 = Vt;

		CSatOrbit helper;
		int ncounter = 0;
		do
		{
			ncounter++;
			double avgscale = 0;
			for (int i = 0; i < calculationsize; i++)
			{

				///换算出真实距离
				//double  S1 = sqrt(pow(toplines[i].OriX - alignlines[i].OriX, 2.0) + pow(toplines[i].OriY - alignlines[i].OriY, 2.0) + pow(toplines[i].OriZ - alignlines[i].OriZ, 2.0));
				///解算错误
				cout.precision(20);
				double  S1 = sqrt(pow(toplines[i].OriX() - alignlines[i].OriX(), 2.0) + pow(toplines[i].OriY() - alignlines[i].OriY(), 2.0) + pow(toplines[i].OriZ() - alignlines[i].OriZ(), 2.0));

				double  S2 = sqrt(pow(toplines[i].X() - alignlines[i].X(), 2.0) + pow(toplines[i].Y() - alignlines[i].Y(), 2.0) + pow(toplines[i].Z() - alignlines[i].Z(), 2.0));
				cout << toplines[i].X() << "\t" << toplines[i].Y() << "\t" << toplines[i].Z() << endl;

				avgscale += S2 / S1 / calculationsize;


				if (false)
				{
					cout << "检查: 成像距离解算" << endl;

					double S22 = S2 / (alignlines[i].ImgY() - toplines[i].ImgY()) * 10786.0;
					cout << "真实距离：" << S22 << endl;
					cout << "Kappa:!" << Kappa * 180.0 / CV_PI << "\t" << "K0:" << K0 * 180.0 / CV_PI << "V0:" << V0 << endl;
					double S1 = sqrt(pow(fabs(V0 * cos(K0) / (cos(Kappa))) * T0, 2.0) + pow(srange, 2.0));
					cout << "实际模型计算距离:" << S1 << endl;

					cout << "!!!, 是否继续" << endl;

					int pause_for_range;
					cin >> pause_for_range;

					coeffarray[i * 2 + 0] = (-sin(K0)*1.0 / tan(Kappa - K0) + cos(K0) * 1.0 / pow(sin(Kappa - K0), 2.0))*V0 * T0 / cos(Kappa);
					coeffarray[i * 2 + 1] = V0 * cos(K0) / cos(Kappa) * cos(Kappa - K0);

					staticarray[i] = S22 - S1;
				}
				///列出误差方程
			}

			cout << "摆扫比例误差为:" << avgscale << endl;
			cout << "!!!, 是否继续" << endl;

			int pause_for_range;
			cin >> pause_for_range;


			helper.transpose(coeffarray, tarray, calculationsize, 2);
			helper.mult(tarray, coeffarray, CTC, 2, calculationsize, 2);
			helper.mult(tarray, staticarray, CTL, 2, calculationsize, 1);

			if (!helper.invers_matrix(CTC, 8))		 					//求逆
			{
				cout << "不可逆" << endl;
				;
			}

			double V[2];

			helper.mult(CTC, CTL, V, 2, 2, 1);
			//memcpy(x, al, sizeof(double)* 8);
			cout.precision(20);

			cout << "解算的安装角改正值： " << V[0] << "改正幅度" << V[0] / K0 << endl;
			cout << "结算的成像时间改正值" << V[1] << "改正幅度" << V[1] / T0 << endl;
			cout << "- -!!!!" << endl;
			int pause_for_range_1;
			cin >> pause_for_range_1;

			K0 = K0 + V[0];
			T0 = T0 + V[1];

		} while (ncounter < 10);

		Vt = T0;

	}



protected:
	ImageModel* p_InfoCarrier;

public:
	~CalibratedModel();
};

