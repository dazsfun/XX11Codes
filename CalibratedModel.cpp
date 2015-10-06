#include "stdafx.h"
#include "CalibratedModel.h"

void CalibratedModel::Calibration_Swing(vector<ImgPoint> heinow, ImgPoint base_p,int start_line,int end_line,int steps)
{
	CSatOrbit helper;

	///摆扫角非线性校正
	int rooms = (end_line - start_line) / (steps)+1;
	cout << "摆扫非线性校正组个数:" << rooms << endl;
	Divide_Swing(heinow, rooms, start_line, end_line, steps);
	ofstream ofs1("c:\\Temp\\match\\heimatch_swing.xls");
	for (int i = 0; i < rooms - 1; i++)
	{
		vector<ImgPoint> upper_line;
		vector<ImgPoint> down_line;

		Link_Swing_Division(i, upper_line, down_line);

		int size_observation = upper_line.size();
		cout << "第" << i << "行所在分段用于最小二乘的观测值个数为" << endl;
		cout << size_observation << endl;

		if (size_observation < 3)
		{
			cout << "第" << i << "行所在分段由于观测数不足未进行最小二乘结算解算" << endl;
			continue;
		}

		double* coeffarray = new double[size_observation  * 2];
		double* staticarray = new double[size_observation ]; // 观测值
		double* tarray = new double[size_observation * 2];
		double CTC[4] = { 0 };
		double CTL[4] = { 0 };


		int ncounter = 0; 
		double avg_scale = 0;
		double scale_rms = 0;
		for (int j = 0; j < size_observation; j++)
		{
			double distance_left = upper_line[j].ImgY() - down_line[j].ImgY();
			double distance_right = upper_line[j]._rotation5 - down_line[j]._rotation5; ///观测值 
			avg_scale += distance_right / distance_left / size_observation;
		}

		double max_scale = -100;
		double min_scale = 100;
		for (int j = 0; j < size_observation; j++)
		{
			double distance_left = upper_line[j].ImgY() - down_line[j].ImgY();
			double distance_right = upper_line[j]._rotation5 - down_line[j]._rotation5; ///观测值 
			max_scale = ((distance_right / distance_left) > max_scale) ? distance_right / distance_left : max_scale;
			min_scale = ((distance_right / distance_left) < min_scale) ? distance_right / distance_left : min_scale;

			scale_rms += ((distance_right / distance_left - avg_scale)*(distance_right / distance_left - avg_scale)) / (size_observation - 1);
		}

		scale_rms = sqrt(scale_rms);
		int count_valid = 0;
		double avg_error_valid = 0;
		for (int j = 0; j < size_observation; j++)
		{
			double distance_left = upper_line[j].ImgY() - down_line[j].ImgY();
			double distance_right = upper_line[j]._rotation5 - down_line[j]._rotation5; ///观测值 
			max_scale = ((distance_right / distance_left) > max_scale) ? distance_right / distance_left : max_scale;
			min_scale = ((distance_right / distance_left) < min_scale) ? distance_right / distance_left : min_scale;

			if (fabs(distance_right / distance_left - avg_scale) < 1.5 * scale_rms)
			{

				count_valid++;
				avg_error_valid += distance_right / distance_left;
			}
			else if (fabs(distance_right / distance_left - avg_scale) > 1.5 * scale_rms)
			{
				cout << "upper_line:(" << upper_line[j].ImgX() << "," << upper_line[j].ImgY() << ")"
					<< "=>" << "(" << upper_line[j]._rotation4 << ", " << upper_line[j]._rotation5 << ")"<<endl;
				cout << "down_line:(" << down_line[j].ImgX() << "," << down_line[j].ImgY() << ")"
					<< "=>" << "(" << down_line[j]._rotation4 << ", " << down_line[j]._rotation5 << ")" << endl;
				cout << "distance_left" << distance_left << endl;
				cout << "distance_right" << distance_right << endl;
				cout << distance_right / distance_left << "\t" << avg_scale << "\t" << scale_rms << endl;
				cout << "正确？" << endl;
				int correction_1;
				cin >> correction_1;
			}
			
		}

		avg_error_valid /= count_valid;

		double kvalue[2] = { avg_scale, avg_scale };
		double kappa[2] = { 1.4 * CV_PI / 180.0, 1.4 * CV_PI / 180.0 };
		cout.precision(20);
		double avg_error = 0;

		do
		{
			ncounter++;

			avg_error = 0;

			for (int j = 0; j < size_observation; j++)
			{
				double distance_left = upper_line[j].ImgY() - down_line[j].ImgY();
				double distance_right = upper_line[j]._rotation5 - down_line[j]._rotation5; ///观测值 
				//avg_scale += distance_right / distance_left / size_observation;
				coeffarray[j * 2 + 0] = (distance_left * cos(kappa[0]));
				coeffarray[j * 2 + 1] = (distance_left)* kvalue[0] * (-sin(kappa[0]));

				staticarray[j] = distance_right -  distance_left * cos(kappa[0]) * kvalue[0];
				avg_error += fabs(staticarray[j]) /static_cast<double>( size_observation);
			}

			helper.transpose(coeffarray, tarray, size_observation, 2);
			helper.mult(tarray, coeffarray, CTC, 2, size_observation , 2);
			helper.mult(tarray, staticarray, CTL, 2, size_observation, 1);

			//CTL为2X1
			if (!helper.invers_matrix(CTC, 2))		 					//求逆
			{
				cout << "第" << ncounter << "次迭代不可逆" << endl;
				break;
			}

			double V[2];


			helper.mult(CTC, CTL, V, 2, 2, 1);
			kvalue[0] += V[0];
			kappa[0] += V[1];

			if (avg_error < 0.000000000001)
			{
				cout << "满足收敛条件：本次解算迭代次数为:" << ncounter << endl;
				cout << "距离为：" << avg_error << endl;
				cout << "平均比例为：" << avg_scale << endl;
				break;
			}
		} while (ncounter < 20);

		delete[] coeffarray;
		delete[] staticarray;
		delete[] tarray;

		cout << "摆扫行" << safestring::Findneighbour(upper_line[0].ImgY()) << "所在行的摆扫比率为:" << endl;
		cout << kvalue[0] << "," << kappa[0] << endl;
		cout << "平均距离:" << avg_error << endl;
		cout << "平均比例为：" << avg_scale << endl;
		cout << "比例中误差为:" << scale_rms << endl;
		cout << "最大比例：" << max_scale << "最小比例:" << min_scale << endl;
		cout << "良好的观测个数为:" << count_valid << endl;
		cout << "真实平均比例：" << endl;
		cout << avg_error_valid << endl;
		int pause_for_debug = 0;
		int i_index = (static_cast<int>(safestring::Findneighbour(upper_line[0].ImgY()) + 2.5) / 5) * 5;
		ofs1 << i_index << avg_error << "\t" << avg_scale << "\t" << scale_rms << "\t" << count_valid << "\t" << avg_error_valid << endl;
		cin >> pause_for_debug;

	}

	ofs1.close();
}

void CalibratedModel::Link_Swing_Division(int index, vector<ImgPoint>& uppper_line, vector<ImgPoint>& down_line)
{
	if (index < 0 || index > swing_division.size() - 2)
	{
		cout << "摆扫角分组超出索引范围,请修正!" << endl;
		return;
	}

	int sizeofall = swing_division[index].size();

	for (int i = 0; i < sizeofall; i++)
	{
		int sizeofnext = swing_division[index + 1].size();

		for (int j = 0; j < sizeofnext; j++)
		{
			double d1 = safestring::Findneighbour(swing_division[index + 1][j].ImgX());
			double d2 = safestring::Findneighbour(swing_division[index][i].ImgX());
			double d_1 = swing_division[index + 1][j].ImgY() - swing_division[index][i].ImgY();


			if (fabs(d1 - d2) < 2&& fabs(d_1) >3)
			{
				uppper_line.push_back(ImgPoint(swing_division[index][i], false));
				down_line.push_back(ImgPoint(swing_division[index + 1][j], false));
			}
		}
	}
}

void CalibratedModel::Divide_Swing( vector<ImgPoint>& heinow, int rooms,int start_line,int end_line,int steps)
{
	swing_division.reserve(rooms);
	swing_division.resize(rooms);

	for (int i = 0; i < heinow.size(); i++)
	{
		///计算分组公式    取整 ( ((匹配点行数 - 起始行编号) + (浮点型)(步长)/2) ) /(步长 ))
		///举例        ((22.5 - 20) + 2.5) / (5) = 1
		///            ((22.4 - 20) + 2.5 ) / (5) = 0
		///            ((17.6 - 20) + 2.5) /(5) = 0
		///            ((17.4 - 20) + 2.5) /(5) = -1  (注 : 小于0和大于rooms -1 时舍弃)
		int division_index = static_cast<int>((heinow[i].ImgY() - start_line + static_cast<double>(steps) / 2.0) / static_cast<double>(steps));
		if (division_index < 0 || division_index > rooms - 1)
		{
			cout << "匹配点被舍弃,超出分类索引组的范畴" << endl;
		}
		else
		{
			swing_division[division_index].push_back(ImgPoint(heinow[i],false));
		}
	}
}


CalibratedModel::CalibratedModel(ImageModel* preparedinfo) : p_InfoCarrier(preparedinfo)
{
	SetPropertys();
	AScan2Cameara0Roll = AScan2Cameara1Pitch = AScan2Cameara2YAW
		= BCamera2Install0Roll = BCamera2Install1Pitch = BCamera2Install2YAW
		= CInstall2Body0Roll = CInstall2Body1Pitch = CInstall2Body2YAW
		= DCalibratedPara0Roll = DCalibratedPara1Pitch = DCalibratedPara2YAW
		= EDateTimeString = EControlPointErrorX = EControlPointErrorY = 0;
}


CalibratedModel::~CalibratedModel()
{
}

void CalibratedModel::SetPropertys()
{
	/*
		SetProperty(0 + 15, asDouble, &GPSUTC);
		*/
	objectname = "校正外畸变成像模型";
	/*
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
	*/
	SetProperty(0 + 45, asDouble, &AScan2Cameara0Roll);
	SetProperty(1 + 45, asDouble, &AScan2Cameara1Pitch);
	SetProperty(2 + 45, asDouble, &AScan2Cameara2YAW);
	SetProperty(3 + 45, asDouble, &BCamera2Install0Roll);
	SetProperty(4 + 45, asDouble, &BCamera2Install1Pitch);
	SetProperty(5 + 45, asDouble, &BCamera2Install2YAW);
	SetProperty(6 + 45, asDouble, &CInstall2Body0Roll);
	SetProperty(7 + 45, asDouble, &CInstall2Body1Pitch);
	SetProperty(8 + 45, asDouble, &CInstall2Body2YAW);
	SetProperty(9 + 45, asDouble, &DCalibratedPara0Roll);
	SetProperty(10 + 45, asDouble, &DCalibratedPara1Pitch);
	SetProperty(11 + 45, asDouble, &DCalibratedPara2YAW);
	SetProperty(12 + 45, asDouble, &EDateTimeString);
	SetProperty(13 + 45, asDouble, &EControlPointErrorX);
	SetProperty(14 + 45, asDouble, &EControlPointErrorY);
}


#include "ap.h"
#include "interpolation.h"

void ClearVectorSame_(vector<Point2d>& _content)
{
	vector<Point2d>::iterator itr = _content.begin();
	do
	{
		if ((*itr).y == (*(itr + 1)).y)
		{
			itr = _content.erase(itr);
		}
		else
			itr++;//这里迭代器会更新出错
	} while ((itr + 1) != _content.end());
}

struct myclass {
	bool operator() (cv::Point2d pt1, cv::Point2d pt2) { return (pt1.x < pt2.x); }
} myobject;



double CalibratedModel::CalibrateTime(const vector<Point2d>& arraypoint,bool bxarray)
{
	///获得均值与中误差
	int calibsize = arraypoint.size();
	double avg_error = 0;
	for (int i = 0; i < calibsize; i++)
	{
		avg_error += arraypoint[i].y / calibsize;
	}

	double rms_error = 0;
	for (int i = 0; i < calibsize; i++)
	{
		rms_error += pow(arraypoint[i].y - avg_error, 2.0) / (calibsize - 1);
	}

	rms_error = sqrt(rms_error);

	cout << "平均值为：" << avg_error << endl;
	cout << "中误差为:" << rms_error << endl;
	vector<Point2d> calipoint;
	ofstream ofs;
	ofs.open("c:\\Temp\\match\\hei_wrong.txt");
	for (int i = 0; i < calibsize; i++)
	{
		if (fabs(arraypoint[i].y - avg_error) > 3 * rms_error)
		{
			///

			ofs << arraypoint[i].x << "\t " << arraypoint[i].y << endl;
			continue;
		}
		else
			calipoint.push_back(arraypoint[i]);
	}

	std::sort(calipoint.begin(), calipoint.end(), myobject);
	ClearVectorSame_(calipoint);

	cout << "拟合总点数::" << calipoint.size() << endl;

	int step_i = calipoint.size();
	ofstream ofsLsmRe("c:\\Temp\\match\\拟合结果.xls");

	vector<TimeCalibStruct> my_inter;
	for (int j = 0; j < step_i; j += 0)
	{
		cout << j << "\\" << step_i << "\r";
		int startindex = j;
		int endindex = startindex + 5;
		double error_avg[2], error_rms[2], error_max[2];
		int error_type[2] = { -1, -1 };
		barycentricinterpolant by_inter[2];
		do
		{
			vector<Point2d> _calitemp;
			for (int c = startindex; c < endindex; c++)
			{
				_calitemp.push_back(calipoint[c]);
			}

			char inputvalue[1024];
			int casize = _calitemp.size();
			string xinitial = "[";
			string yinitial = "[";
			for (int i = 0; i < casize; i++)
			{
				sprintf_s(inputvalue, 512, "%lf,", _calitemp[i].x);

				xinitial += string(inputvalue);
				sprintf_s(inputvalue, 512, "%lf,", _calitemp[i].y);
				yinitial += string(inputvalue);
			}
			xinitial = xinitial.substr(0, xinitial.size() - 1);
			yinitial = yinitial.substr(0, yinitial.size() - 1);
			xinitial += "]";
			yinitial += "]";

			real_1d_array x = xinitial.c_str();
			real_1d_array y = yinitial.c_str();	ae_int_t info2;


			double rho = 0.0;
			barycentricinterpolant x2y;
			polynomialfitreport polyreport;
			int infos;
			ae_int_t hei_pl;
			//polynomialbuild(x, y, _theCurve.yfromxbarycentric);
			//polynomialbuild(y, x, _theCurve.xfromybarycentric);


			barycentricinterpolant byx2y;
			barycentricfitreport byreport;
			ae_int_t byinfo;
			barycentricfitfloaterhormann(x, y, x.length(), 5, byinfo, byx2y, byreport);
			

			polynomialfit(x, y, x.length(), pownumber, hei_pl, x2y, polyreport);



			double  cho_value = 5;
			ae_int_t hei_ae;
			spline1dinterpolant x2yspline;
			spline1dfitreport x2ysplinereport;
			spline1dfitpenalized(x, y, x.length(), cho_value, hei_ae, x2yspline, x2ysplinereport);



			double minerror = (byreport.rmserror < polyreport.rmserror) ? byreport.rmserror : polyreport.rmserror;
			minerror = (minerror < x2ysplinereport.rmserror) ? minerror : x2ysplinereport.rmserror;

			if (true) //
			{
				if (byreport.rmserror < polyreport.rmserror && byreport.rmserror < x2ysplinereport.rmserror)
				{
					error_type[1] = 0;
					error_rms[1] = byreport.rmserror;
					error_avg[1] = byreport.avgerror;
					error_max[1] = byreport.maxerror;
				}
				else if (polyreport.rmserror < byreport.rmserror && polyreport.rmserror < x2ysplinereport.rmserror)
				{
					error_type[1] = 1;
					error_rms[1] = polyreport.rmserror;
					error_avg[1] = polyreport.avgerror;
					error_max[1] = polyreport.maxerror;
				}
				else if (x2ysplinereport.rmserror < polyreport.rmserror && x2ysplinereport.rmserror < byreport.rmserror)
				{
					error_type[1] = 2;
					error_rms[1] = x2ysplinereport.rmserror;
					error_avg[1] = x2ysplinereport.avgerror;
					error_max[1] = x2ysplinereport.maxerror;
				}
			}

			if (minerror < 1 && error_avg[1] < 1)
			{
				error_type[0] = error_type[1];
				error_rms[0] = error_rms[1];
				error_max[0] = error_max[1];
				error_avg[0] = error_avg[1];
			}

			else
			{
				TimeCalibStruct calib;
				calib._startvalue = calipoint[startindex].x;
				calib._endvalue = calipoint[endindex - 1].x;
				calib._interkernel = x2y;
				my_inter.push_back(calib);
				break;
			}

			endindex++;
		} while (endindex < calipoint.size());

		ofsLsmRe.precision(15);
		ofsLsmRe << error_type[0] << "\t" << error_rms[0] << "\t" << error_avg[0] << "\t" << error_max[0] << "\t" << endindex - 1 - startindex << endl;

		j = endindex - 1;
	}

	int groupsize = my_inter.size();
	int startcheckindex = my_inter[0]._startvalue;
	int endcheckindex = my_inter[groupsize - 1]._endvalue;

	ofstream curveinfo;
	if(bxarray)
		curveinfo.open("c:\\Temp\\match\\畸变曲线_x.xls");
	else
		curveinfo.open("c:\\Temp\\match\\畸变曲线.xls");

	curveinfo.precision(15);
	for (int i = startcheckindex; i < endcheckindex; i++)
	{
		double getvalue = 0;
		for (int j = 0; j < groupsize; j++)
		{
			if (i > my_inter[j]._startvalue && i < my_inter[j]._endvalue)
			{
				double intervalue = i + 0.01;
				getvalue = barycentriccalc(my_inter[j]._interkernel, i + 0.01);
				if (fabs(fabs(getvalue) - avg_error) > 3 * rms_error)
				{
					cout << i << "\t" << getvalue << "\t" << my_inter[j]._startvalue << "\t" << my_inter[j]._endvalue << endl;
					cout << barycentriccalc(my_inter[j]._interkernel, my_inter[j]._startvalue) << endl;
					cout << barycentriccalc(my_inter[j]._interkernel, my_inter[j]._endvalue) << endl;


					int pause_new;
					cin >> pause_new;
				}

				break;
			}

			else if (j < groupsize - 1 && i > my_inter[j]._endvalue && i < my_inter[j + 1]._startvalue)
			{
				double calleft = barycentriccalc(my_inter[j]._interkernel, my_inter[j]._endvalue);
				double calright = barycentriccalc(my_inter[j + 1]._interkernel, my_inter[j + 1]._startvalue);

				double distance = my_inter[j + 1]._startvalue - my_inter[j]._endvalue;
				getvalue = calleft * (my_inter[j + 1]._startvalue - i) / distance + calright * (i - my_inter[j]._endvalue) / distance;
				break;
			}
		}

		curveinfo << i << "\t" << getvalue << endl;
	}

	cout << "畸变拟合结束" << endl;

	if (!bxarray)
	{
		for (int i = 0; i < my_inter.size(); i++)
		{
			_arryay_inter.push_back(my_inter[i]);
		}
	}

	else
	{
		for (int i = 0; i < my_inter.size(); i++)
		{
			_arryay_inter_x.push_back(my_inter[i]);

		}
	}

	cout << "畸变拷贝结束" << endl;

	curveinfo.close();
	ofsLsmRe.close();
}






#include "GeoBase.h"
void CalibratedModel::Calibration(vector<ImgPoint> m_points_1, ImgPoint hei)
{
	vector<ImgPoint> m_points;
	ImgPoint heinow = hei;
	cout << "开始选取外检校控制点" << endl;
	for (int i = 48; i < 480 -48; i+=20)
	{
		for(int j = 50; j < 10000; j += 300)
		{
			for (int q = 0; q < m_points_1.size(); q++)
			{
				if (m_points_1[q].ImgX()> i && m_points_1[q].ImgX() < i + 20 && m_points_1[q].ImgY()> j
					&&m_points_1[q].ImgY() < j + 300)
				{
					heinow / m_points_1[q];
					m_points.push_back(heinow);
					break;
				}
			}
		}
	}

	int pointsize = m_points.size();
	CSatOrbit helper;
	double matrix1[9], matrix2[9], matrixtemp1[9], matrix3[9], matrires[9];
	helper.rot(AScan2Cameara0Roll, AScan2Cameara1Pitch, AScan2Cameara2YAW, matrix1);
	helper.rot(BCamera2Install0Roll, BCamera2Install1Pitch, BCamera2Install2YAW, matrix2);
	helper.rot(CInstall2Body0Roll, CInstall2Body1Pitch, CInstall2Body2YAW, matrix3);
	helper.mult(matrix1, matrix2, matrixtemp1, 3, 3, 3);
	helper.mult(matrixtemp1, matrix3, matrires, 3, 3, 3);
	//pointsize = pointsize - 2;
	vector<double*> lightq; lightq.reserve(pointsize); lightq.resize(pointsize);
	vector<double*> lightp; lightp.reserve(pointsize); lightp.resize(pointsize);
	//vector* matrixarray = new double*[(m_points).size()];

	double avgp[3] = { 0 }; double avgq[3] = { 0 };
	double avg_dist_ = 0;
	for (int i = 0; i < pointsize; i++)
	{
		double rotates[9] = { 0 };
		lightq[i] = new double[3];
		lightp[i] = new double[3];
		double timenow = m_points[i].utcapi.GetUTC();
		double position[3];
		double body2WGS84Rotation[9] = { 0 };
		p_InfoCarrier->IBody2WGS84Test(timenow, body2WGS84Rotation, position);

		double actGroundPos[3] = { m_points[i].OriX(), m_points[i].OriY(), m_points[i].OriZ() };
		double distance = sqrt(pow(m_points[i].OriX() - m_points[i].X(), 2.0) + pow(m_points[i].OriY() - m_points[i].Y(), 2.0) + pow(m_points[i].OriZ() - m_points[i].Z(), 2.0));
		avg_dist_ += fabs(distance) / pointsize;
		//cout << "distance is:" << distance << endl;
		m_points[i].TrueVector(position, lightq[i]);
		m_points[i].ActualVector(position, actGroundPos, lightp[i]);
		///实际光线则是指有姿轨数据计算出的光线
		///得到指向角
		/*
		cout << "lightp is:" << endl;
		for (int j = 0; j < 3; j++)
		{
			cout << lightp[i][j] << endl;
		}

		cout << "lightq is:" << endl;
		for (int j = 0; j < 3; j++)
		{
			cout << lightq[i][j] << endl;
		}
		*/
		for (int j = 0; j < 3; j++)
		{
			avgp[j] += lightp[i][j] / pointsize;
			avgq[j] += lightq[i][j] / pointsize;
		}
	}



	///接下来要做的工作是
	///Least-Squares Rigid Motion Using SVD
	///通过 最小二乘来估算真实光线与实际光线间的偏转
	///算法原理来源:
	///http://igl.ethz.ch/projects/ARAP/svd_rot.pdf
	///算法封装在动态苦苦OpenCVHelper
	OpenCVHelper cvhelper;

	double transCompensate[3] = { 0 };
	double compensateVe[3] = { 0 };
	double rotCompensate[9] = { 0 };
	cout << "lightp size:" << lightp.size() << endl;
	cout << "lightq size is:" << lightq.size() << endl;
	cvhelper.FindRT(lightp, lightq, avgp, avgq, rotCompensate, transCompensate);


	Mat Rot(3, 3, CV_64F, &rotCompensate);
	ofstream ofs1("fs.txt");
	ofs1 << Rot << endl;
	ofs1.close();
	CvMat value1test = Rot;
	Mat rotationVe(3, 1, CV_64F, &compensateVe);
	CvMat value2test = rotationVe;
	cvRodrigues2(&value1test, &value2test);

	cout << avgp[0] << "\t" << avgp[1] << "\t" << avgp[2] << endl;
	cout << avgq[0] << "\t" << avgq[1] << "\t" << avgq[2] << endl;

	cout << "Avg dist is" << avg_dist_ << endl;
	cout << "Transcompensate is:" << transCompensate[0] << "\t" << transCompensate[1] << "\t" << transCompensate[2] << endl;
	cout << "ROLLCompensate:" << compensateVe[0] * 180.0 / CV_PI << endl << "PitchCompensate" << compensateVe[1] * 180.0 / CV_PI
		<< endl
		<< "YAWCompensate:" << compensateVe[2] * 180.0 / CV_PI << endl;
	//helper.rot2eulor(rotCompensate, DCalibratedPara0Roll, DCalibratedPara1Pitch, DCalibratedPara2YAW);

}