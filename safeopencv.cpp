#include "stdafx.h"
#include "safeopencv.h"
#include <opencv2/core/core.hpp>        // Basic OpenCV structures


using namespace cv;

#ifndef FLT_EPSILON
#define FLT_EPSILON     1.192092896e-07F        /* smallest such that 1.0+FLT_EPSILON != 1.0 */
#endif

#ifndef _ASSERT
#define _ASSERT
#endif


void SubPixFitParabola(cv::Point2d* Result, cv::Point2d& P1, cv::Point2d& P2, cv::Point2d& P3);


double  safeopencv::MinMaxloc_subpixel(CV_OUT cv::Point2d* SubPixLoc,
	cv::InputArray src,
	CV_IN_OUT cv::Point* LocIn,
	CV_IN_OUT const int Method
	)
{/*

 in
 src = Heat map. It must be single-channel 32-bit floating-point
 LocIn = X,y location you want to start with. This is normally  from MinMaxLoc, however you can explore other peaks.
 Method = -1 to disable math for testing, 0 for the best we have, 1 for this parabola curve fit version
 out
 SubPixLoc = location adjusted with improved precision.

 notes:
 If you rescale the heat map after match template... you must also adjust LocIn to match.

 */

	// set default result in case we bail
	SubPixLoc->x = (float)LocIn->x;
	SubPixLoc->y = (float)LocIn->y;

	if (Method == -1) { // disable any changes for testing		
		return 0;
	};

	// At this time we don't have anything other than Parabola math so we can ignore "Method".

	{ // Parabola math	

		/*
		The values returned from MatchTemplate are not linear past the point where it just starts to drop.
		The reason is that we can not assume that the template is the same on all sides. Imagine a sloped edge on one side and a sharp edge on the other.

		We can also get several values at the top that are all the same since the source and the template are integers.

		We also have to protect against the situation where the source is a solid white or solid black. The result is a constant value heat map.
		*/
		Mat HeatMap = src.getMat();

		// pick some limiting values
		// It's not expected that the template peak values will span more than 1/16 of the source size. This also limits the processing time used when looking at a blank image.
		Size MaxScan; MaxScan.width = HeatMap.cols >> 4;  MaxScan.height = HeatMap.rows >> 4;


		Point ScanRectMin;  // I used two Points instead of a Rect to prevent having Rect compute right/left values in each loop below
		Point ScanRectMax;
		ScanRectMin.x = LocIn->x - MaxScan.width; if (ScanRectMin.x < 0) ScanRectMin.x = 0;
		ScanRectMin.y = LocIn->y - MaxScan.height; if (ScanRectMin.y < 0) ScanRectMin.y = 0;
		ScanRectMax.x = LocIn->x + MaxScan.width; if (ScanRectMax.x >= HeatMap.cols) ScanRectMax.x = HeatMap.cols - 1;
		ScanRectMax.y = LocIn->y + MaxScan.height; if (ScanRectMax.y >= HeatMap.rows) ScanRectMax.y = HeatMap.rows - 1;

		// were we are starting at
		const float FloatValueChange = FLT_EPSILON * 10.0f; // smallest change that we can do math on with some meaningful result.

		// scan to find area to use. this can get complicated since we may be given a point near any of the edges of the blob we want to use.		
		float SrcStartingPoint = HeatMap.at<float>(LocIn->y, LocIn->x);

		Point Center = *LocIn;

		// results
		Point ScanRight;
		Point ScanLeft;
		Point ScanUp;
		Point ScanDown;

		//for (int rescan = 0; rescan < 2; ++rescan){

		ScanRight = Center;
		while (true){
			++ScanRight.x; // no point checking the passed location. so inc first
			if (ScanRight.x > ScanRectMax.x){
				//					_ASSERT(0);
				return 1; // ran out of room to scan
			};
			float Val = HeatMap.at<float>(ScanRight.y, ScanRight.x);
			if (abs(Val - SrcStartingPoint) > FloatValueChange){
				break;
			};
		};

		ScanLeft = Center;
		while (true){
			--ScanLeft.x; // no point checking the passed location. so inc first
			if (ScanLeft.x < ScanRectMin.x){
				//					_ASSERT(0);
				return 1; // ran out of room to scan
			};
			if (abs(HeatMap.at<float>(ScanLeft.y, ScanLeft.x) - SrcStartingPoint) > FloatValueChange){
				break;
			};
		};

		ScanUp = Center;
		while (true){
			++ScanUp.y; // assume G cords. The actual direction of Up in the image is not important since the math is symmetrical
			if (ScanUp.y > ScanRectMax.y){
				//					_ASSERT(0);
				return 1; // ran out of room to scan
			};
			if (abs(HeatMap.at<float>(ScanUp.y, ScanUp.x) - SrcStartingPoint) > FloatValueChange){
				break;
			};
		};

		ScanDown = Center;
		while (true){
			--ScanDown.y; // assume G cords. The actual direction of Up in the image is not important since the math is symmetrical
			if (ScanDown.y < ScanRectMin.y){
				//					_ASSERT(0);
				return 1; // ran out of room to scan
			};
			if (abs(HeatMap.at<float>(ScanDown.y, ScanDown.x) - SrcStartingPoint) > FloatValueChange){
				break;
			};
		};

		// At this point we have a good starting point on the blob area, but our initial scan may be stuck on one side so center and rescan once more

		//Center.x =  ((ScanRight.x - ScanLeft.x) >> 1) + ScanLeft.x;
		//Center.y =  ((ScanUp.y    - ScanDown.y) >> 1) + ScanDown.y;

		// did center change?
		//if ((Center.x == LocIn->x) && (Center.y == LocIn->y)) break; // done early

		//}; // for rescan

		// measure polarity if needed



		// At this point we have a center of a blob with some extents to use

		// for each axis we now do a triangulation math.


		// imagine the match numbers as height and the pixel numbers as horizontal.

		//B is highest, A and C are on the sides


		double ErrorVal = 0;

		{// X axis

			Point2d A;
			A.x = ScanLeft.x; // The pixel cords
			A.y = HeatMap.at<float>(ScanLeft.y, ScanLeft.x); // the Heat map value

			Point2d B; // center
			B.x = Center.x; // The pixel cords
			B.y = HeatMap.at<float>(Center.y, Center.x); // the Heat map value

			Point2d C;
			C.x = ScanRight.x; // The pixel cords
			C.y = HeatMap.at<float>(ScanRight.y, ScanRight.x); // the Heat map value

			Point2d Result;
			SubPixFitParabola(&Result, A, B, C);
			// we throw away the y and use the x

			// clip and set error
			if (Result.x < ScanLeft.x){
				_ASSERT(0);
				Result.x = ScanLeft.x;
				ErrorVal = 1;
			};
			if (Result.x > ScanRight.x){
				_ASSERT(0);
				Result.x = ScanRight.x;
				ErrorVal = 1;
			};
			SubPixLoc->x = Result.x;
		}; // X axis



		{// Y axis

			// this time we swap x and y since the parabola is always found in the x
			Point2d A;
			A.x = ScanDown.y; // The pixel cords
			A.y = HeatMap.at<float>(ScanDown.y, ScanDown.x); // the Heat map value

			Point2d B; // center
			B.x = Center.y; // The pixel cords
			B.y = HeatMap.at<float>(Center.y, Center.x); // the Heat map value

			Point2d C;
			C.x = ScanUp.y; // The pixel cords
			C.y = HeatMap.at<float>(ScanUp.y, ScanUp.x); // the Heat map value

			Point2d Result;
			SubPixFitParabola(&Result, A, B, C);
			// we throw away the y and use the x
			Result.y = Result.x;

			// clip and set error
			if (Result.y < ScanDown.y){
				_ASSERT(0);
				Result.y = ScanDown.y;
				ErrorVal = 1;
			};
			if (Result.y > ScanUp.y){
				_ASSERT(0);
				Result.y = ScanUp.y;
				ErrorVal = 1;
			};
			SubPixLoc->y = Result.y;
		}; // X axis


		return ErrorVal;


	}; // Bill's Tilt math

	return 0;

};

// Parabolic fit
void SubPixFitParabola(cv::Point2d* Result, cv::Point2d& P1, cv::Point2d& P2, cv::Point2d& P3)
{/*
 Parabola fit and resulting peak

 The parabola is aligned along the X axis with the peak being in the Y.

 in
 P1 = a point on one side
 P2 = the center point
 P3 = a point on the other side
 out
 Result = the peak point in the center of the parabola
 */

	Result->x = P2.x; // default in case of an error
	Result->y = P2.y;


	/* from http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	This is really just a simple linear algebra problem, so you can do the calculation symbolically. When you substitute in the x and y values of your three points, you'll get three linear equations in three unknowns.

	A x1^2 + B x1 + C = y1
	A x2^2 + B x2 + C = y2
	A x3^2 + B x3 + C = y3

	The straightforward way to solve this is to invert the matrix

	x1^2  x1  1
	x2^2  x2  1
	x3^2  x2  1

	and multiply it by the vector

	y1
	y2
	y3

	The result of this is... okay, not exactly all that simple ;-) I did it in Mathematica, and here are the formulas in pseudocode:
	*/

	double denom = (P1.x - P2.x) * (P1.x - P3.x) * (P2.x - P3.x); // can't be zero since X is from pixel locations.
	double A = (P3.x * (P2.y - P1.y) + P2.x * (P1.y - P3.y) + P1.x * (P3.y - P2.y)) / denom;
	double B = ((P3.x * P3.x) * (P1.y - P2.y) + (P2.x * P2.x) * (P3.y - P1.y) + (P1.x * P1.x) * (P2.y - P3.y)) / denom;
	double C = (P2.x * P3.x * (P2.x - P3.x) * P1.y + P3.x * P1.x * (P3.x - P1.x) * P2.y + P1.x * P2.x * (P1.x - P2.x) * P3.y) / denom;



	// y = A * x^2 + B * x + C 

	//now find the center

	double xv = -B / (2 * A);
	double yv = C - (B*B) / (4 * A);


	Result->x = xv;
	Result->y = yv;
};

void ClearVectorSame(vector<Point2d>& _content)
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

///只支持小数
void safeopencv::Mat_Fill(vector<Point2f>& modelarray, int windowsize, Mat_Fill_Type filltype)
{
	if (filltype == Enum_Mat_Fill_Symmetric)
	{
		int width_ = windowsize;
		int height_ = windowsize;
		modelarray.reserve(windowsize * windowsize);
		modelarray.resize(windowsize * windowsize);
		if (width_ != height_)
		{
			cout << "窗口宽高不一致，无法中心化" << endl;
			return;
		}

		else if (width_ % 2 == 0)
		{
			cout << "窗口宽高为偶数，无法中心化" << endl;
		}

		int halfsize = width_ / 2;

		for (int j = 0; j < width_; j++)
		{
			for (int i = 0; i < width_; i++)
			{
				modelarray[j * width_ + i] = Point2f(-halfsize + i, -halfsize + j);
			}
		}

		return;
	}
}


////没错了,04/19
double safeopencv::Curve_MinMax(vector<Point2d> curveinfo, double& maxvaluereturn, double& minvaluereturn, double& maxloc, double& minloc)
{
	//ClearVectorSame(curveinfo);

	double samevalue = 0;

	string xinitial = "[";
	string yinitial = "[";
	double minvalue = 99999;
	double maxvalue = -99999;
	double minyvalue, maxyvalue;
	minyvalue = 99999;
	maxyvalue = -99999;
	int minyvalueindex, maxyvalueindex;
	minyvalueindex = maxyvalueindex = 0;
	int minyvectorindex = 0;
	minvalue = 0; maxvalue = 4085;
	for (int i = 0; i < curveinfo.size(); i++)
	{
		char inputvalue[512];

		if (curveinfo[i].y < minyvalue)
		{
			minyvalue = curveinfo[i].y;
			minyvalueindex = curveinfo[i].x;
			minyvectorindex = i;
		}

	}

	char inputvalue[512];

	int upindex = (minyvectorindex + 25) > curveinfo.size() ? curveinfo.size() : (minyvectorindex + 25);
	int downindex = (minyvectorindex - 25) < 0 ? 0 : (minyvectorindex - 25);
	for (int i = downindex; i < upindex; i++)
	{
		sprintf_s(inputvalue, 512, "%lf,", curveinfo[i].x);
		xinitial += string(inputvalue);
		sprintf_s(inputvalue, 512, "%lf,", curveinfo[i].y * 10000);
		yinitial += string(inputvalue);
	}
	xinitial = xinitial.substr(0, xinitial.size() - 1);
	yinitial = yinitial.substr(0, yinitial.size() - 1);
	xinitial += "]";
	yinitial += "]";

	real_1d_array x = xinitial.c_str();
	real_1d_array y = yinitial.c_str();

	///这里我们还想找到y方向的最小，最大值，以方便找到核线影像部分的边界
	ae_int_t info2;

	double rho = 0.0;
	barycentricinterpolant x2y;
	polynomialfitreport polyreport;
	int infos;
	//polynomialbuild(x, y, _theCurve.yfromxbarycentric);
	//polynomialbuild(y, x, _theCurve.xfromybarycentric);
	alglib::polynomialfit(x, y, upindex - downindex + 1, infos, x2y, polyreport);



	double summitup, summitdown;
	double summitupindex, summitdownindex;
	summitup = -999;
	summitdown = 99999;
	for (double i = minyvalueindex - 2; i < minyvalueindex + 2; i += 0.05)
	{
		double valuenow = barycentriccalc(x2y, i);
		if (valuenow < summitdown)
		{
			summitdown = valuenow;
			summitdownindex = i;
		}
	}



	minvaluereturn = summitdown;


	minloc = summitdownindex;


	return 0;
}

bool safeopencv::Inside(Rect parent, Rect child)
{
	if (child.x < parent.x || child.y < parent.y)
	{
		return false;
	}

	if (child.x + child.width > parent.x + parent.width || child.y + child.height > parent.y + parent.height)
	{
		return false;
	}

	return true;
}

void safeopencv::Pixel_Warp_Affine(const vector<Point2f>& srcMat, vector<Point2f>& dstMat, const Point2f& pixelstart, double* affine_data)
{
	int sizepoint = srcMat.size();
	dstMat.reserve(sizepoint);
	dstMat.resize(sizepoint);

	for (int i = 0; i < sizepoint; i++)
	{
		dstMat[i].x = affine_data[0] * srcMat[i].x + affine_data[1] * srcMat[i].y + affine_data[2] + pixelstart.x;
		dstMat[i].y = affine_data[3] * srcMat[i].x + affine_data[4] * srcMat[i].y + affine_data[5] + pixelstart.y;
	}
}

double Bilinear_Cubic(Mat& srcMat_, const Point2f& pixel_loc)
{
	int cubic_index_x = pixel_loc.x;
	int cubic_index_y = pixel_loc.y;
	double cubic_shift_x = pixel_loc.x - cubic_index_x;
	double cubic_shift_y = pixel_loc.y - cubic_index_y;
	//cout << "cubic_index_x :" << cubic_index_x << endl;
	//cout << "cubic_index_y :" << cubic_index_y << endl;

	if (cubic_index_x < 0 || cubic_index_y < 0 || cubic_index_x > srcMat_.cols - 2 || cubic_index_y > srcMat_.rows - 2)
	{
		return 0;
	}

	Mat _Cubic = srcMat_(Rect(cubic_index_x, cubic_index_y, 2, 2)).clone();

	void* pvoiddata = (void*)_Cubic.data;
	double* pdata = static_cast<double*>(pvoiddata);

	return pdata[0] * (1 - cubic_shift_x) * (1 - cubic_shift_y) +
		pdata[1] * cubic_shift_x * (1 - cubic_shift_x) +
		pdata[2] * (1 - cubic_shift_x) * cubic_shift_y +
		pdata[3] * cubic_shift_x * cubic_shift_y;
}

const double safeopencv::Resample(Mat& srcMat, Point2f pixel_loc, ResampleKernels* resmaplekernel)
{
	double pixelvalue__ = 0;
	int kernelsize = resmaplekernel->GetCurKernelSize();
	int halfsize = kernelsize / 2;
	Rect imgrange(0, 0, srcMat.cols, srcMat.rows);
	if (!((pixel_loc.x >= 0) && (pixel_loc.y >= 0) && (pixel_loc.x <= srcMat.cols - 1) && (pixel_loc.y <= srcMat.rows - 1)))
	{
		return false;
	}



	if (pixel_loc.x == 0 || pixel_loc.x == srcMat.cols - 1)
	{
		if (pixel_loc.y == 0 || pixel_loc.y == srcMat.rows - 1)
		{
			///角点
			cout << "采样到角点！！" << endl;
			return srcMat.at<double>(pixel_loc.x, pixel_loc.y);  ///在多个两个点的情况下都不要使用opencv的at		                                                   
		}

		///列方向线性
		int yindex = pixel_loc.y;
		double yoffset = pixel_loc.y - yindex;
		return srcMat.at<double>(pixel_loc.x, pixel_loc.y) * (1 - yoffset) + srcMat.at<double>(pixel_loc.x, pixel_loc.y + 1) * yoffset;
	}

	if (pixel_loc.y == 0 || pixel_loc.y == srcMat.rows - 1)
	{
		int xindex = pixel_loc.x; double xoffset = pixel_loc.x - xindex;
		///行方向线性
		return srcMat.at<double>(pixel_loc.x, pixel_loc.y) * (1 - xoffset) + srcMat.at<double>(pixel_loc.x + 1, pixel_loc.y) * xoffset;
	}

	Rect subrange(pixel_loc.x - halfsize + 1, pixel_loc.y - halfsize + 1, kernelsize, kernelsize);
	if (!safeopencv::Inside(imgrange, subrange))
	{
		///双线性
		//cout << "异常，！！ 重采样窗口超出" << endl;
		//cout << pixel_loc.x << "\t" << pixel_loc.y << endl;
	
		int pause_mode = 0;
	//	cin >> pause_mode;
		if ((pixel_loc.x > 0) && (pixel_loc.y > 0) && (pixel_loc.x < srcMat.cols - 1) && (pixel_loc.y < srcMat.rows - 1))
		{
			//cout << "继续？" << endl;
			return  Bilinear_Cubic(srcMat, pixel_loc);
		}

		return 0;

	}

	else
	{
		///使用指定采样核
		Mat subMat = srcMat(subrange).clone();
		void* pvoidsub = static_cast<void*>(subMat.data);
		double* pdoublesub = static_cast<double*>(pvoidsub);
		int windowall = kernelsize * kernelsize;
		double returnvalue = 0;
		double* kerneldata = new double[kernelsize * kernelsize];
		memset(kerneldata, 0, sizeof(double)* kernelsize * kernelsize);
		resmaplekernel->GetCurPixelResampleKernels(pixel_loc.x, pixel_loc.y, kerneldata);

		returnvalue = 0;
		for (int j = 0; j < kernelsize; j++)
		{
			cout.precision(10);
			for (int i = 0; i < kernelsize; i++)
			{
				//cout << pdoublesub[j * kernelsize + i] << "\t" << subMat.at<double>(j, i) <<endl;
				returnvalue += pdoublesub[i * kernelsize + j] * kerneldata[j * kernelsize + i];
			}
		}
		delete[] kerneldata;
		//cout << returnvalue << endl;
		//cout << "像素差异:" << pdoublesub[halfsize * (halfsize -1) + halfsize  -1 ] - subMat.at<double>(static_cast<int>(pixel_loc.x), static_cast<int>(pixel_loc.y));
		return returnvalue;
		//return srcMat.at<double>(static_cast<int>(pixel_loc.x), static_cast<int>(pixel_loc.y));
	}
}

void safeopencv::Rot2Eulor(const double* rotvalue, double& rot1, double& rot2, double& rot3)
{
	double rotcomvalue[3];// = { rot1 / 180.0 * CV_PI, rot2 / 180.0 * CV_PI, rot3 / 180.0 * CV_PI };
	Mat rotatevector(3, 1, CV_64F, &rotcomvalue);

	double rotcom1value[9];
	memcpy(rotcom1value, rotvalue, sizeof(double)* 9);
	CvMat valuecom1 = rotatevector;
	Mat rotationcomMatrix(3, 3, CV_64F, &rotcom1value);
	CvMat valuecom2 = rotationcomMatrix;
	cvRodrigues2(&valuecom2, &valuecom1);

	rot1 = rotcomvalue[0] * 180.0 / CV_PI;
	rot2 = rotcomvalue[1] * 180.0 / CV_PI;
	rot3 = rotcomvalue[2] * 180.0 / CV_PI;
}

void safeopencv::Eulor2Rot(const double& rot1, const double& rot2, const double& rot3, double* rotvalue)
{
	double rotcomvalue[3] = { rot1 / 180.0 * CV_PI, rot2 / 180.0 * CV_PI, rot3 / 180.0 * CV_PI };
	Mat rotatevector(3, 1, CV_64F, &rotcomvalue);

	double rotcom1value[9];
	CvMat valuecom1 = rotatevector;
	Mat rotationcomMatrix(3, 3, CV_64F, &rotcom1value);
	CvMat valuecom2 = rotationcomMatrix;
	cvRodrigues2(&valuecom1, &valuecom2);

	memcpy(rotvalue, rotcom1value, sizeof(double)* 9);
}

///拟合匹配曲线，显示结果 
void safeopencv::Match_Curve(vector<Point2d>& matchleft, vector<Point2d>& matchright, string img_path_left, string img_path_right, barycentricinterpolant& left2right_x, barycentricinterpolant& left2right_y)
{
	if (matchleft.size() != matchright.size())
	{
		cout << "匹配曲线拟合不成功，点集大小不一致" << endl;
		return;
	}
	double samevalue = 0;

	cout << "拟合前匹配点数目:" << endl;
	cout<< matchleft.size() << endl;

	barycentricinterpolant x2x_info, y2y_info;
	do
	{

		string xinitial = "[";
		string yinitial = "[";

		string xinitial_1 = "[";
		string yinitial_1 = "[";
		char inputvalue[512];


		for (int i = 0; i < matchleft.size(); i++)
		{
			sprintf_s(inputvalue, "%lf,", matchleft[i].x);
			xinitial += string(inputvalue);
			sprintf_s(inputvalue, "%lf,", matchright[i].x);
			yinitial += string(inputvalue);

			sprintf_s(inputvalue, "%lf,", matchleft[i].y);
			xinitial_1 += string(inputvalue);
			sprintf_s(inputvalue, "%lf,", matchright[i].y);
			yinitial_1 += string(inputvalue);
			//cout << matchleft[i].x << "\t" << matchright[i].x << "\t" << matchleft[i].y << "\t" << matchright[i].y << endl;
		}
		xinitial = xinitial.substr(0, xinitial.size() - 1);
		yinitial = yinitial.substr(0, yinitial.size() - 1);
		xinitial_1 = xinitial_1.substr(0, xinitial_1.size() - 1);
		yinitial_1 = yinitial_1.substr(0, yinitial_1.size() - 1);
		xinitial += "]";
		yinitial += "]";
		xinitial_1 += "]";
		yinitial_1 += "]";

		real_1d_array x = xinitial.c_str();
		real_1d_array y = yinitial.c_str();
		real_1d_array x_1 = xinitial_1.c_str();
		real_1d_array y_1 = yinitial_1.c_str();

		ae_int_t info2;

		double rho = 0.0;

		cout << "开始进行多项式拟合" << endl;

		///输出拟合信息

		int info_1;
		barycentricinterpolant x2y_3, x2y_4;
		barycentricfitreport report3, report4;
		barycentricfitfloaterhormann(x, y, x.length(), 4, info_1, x2y_3, report3);
		barycentricfitfloaterhormann(x_1, y_1, x_1.length(), 4, info_1, x2y_4, report4);

		cout << "x方向坐标有理曲线拟合结果:" << endl;
		cout << "残差中误差:" << report3.rmserror << ".最大误差:" << report3.maxerror << ".平均误差" << report3.avgerror << endl;

		cout << "y方向坐标有理曲线拟合结果:" << endl;
		cout << "残差中误差:" << report4.rmserror << ".最大误差:" << report4.maxerror << ".平均误差" << report4.avgerror << endl;

		vector<Point2d> _letemp;
		vector<Point2d> _rtTemp;
		for (int i = 0; i < matchleft.size(); i++)
		{
			double xvalue = barycentriccalc(x2y_3, matchleft[i].x);
			if (fabs(xvalue - matchright[i].x) > 3 * report3.rmserror)
				continue;
			double yvalue = barycentriccalc(x2y_4, matchleft[i].y);
			if (fabs(yvalue - matchright[i].y) > 3 * report4.rmserror)
				continue;

			_letemp.push_back(matchleft[i]);
			_rtTemp.push_back(matchright[i]);
		}

		matchleft.clear();
		matchright.clear();

		for (int i = 0; i < _letemp.size(); i++)
		{
			matchleft.push_back(_letemp[i]);
			matchright.push_back(_rtTemp[i]);
		}

		int pause_iteration = 0;
		cout << "是否继续?" << endl;
		cin >> pause_iteration;
		if (pause_iteration)
		{
			x2x_info = x2y_3;
			y2y_info = x2y_4;
			break;
		}
	} while (true);

	cout << "剩余匹配点数目:" << endl;
	cout << matchleft.size() << endl;

	ofstream ofs("c:\\Temp\\match\\heimatchcurve.xls");
	ofs.precision(20);
	for (int i = 0; i < matchleft.size(); i++)
	{
		double xvalue = barycentriccalc(x2x_info, matchleft[i].x);
		double yvalue = barycentriccalc(y2y_info, matchleft[i].y);

		ofs << matchleft[i].x << "\t" << matchleft[i].y << "\t" << matchright[i].x << "\t" <<matchright[i].y << "\t" << xvalue << "\t" << yvalue << "\t" << xvalue - matchleft[i].x << "\t" << yvalue - matchleft[i].y << endl;
	}
	
	ofs.close();

	left2right_x = x2x_info;
	left2right_y = y2y_info;
}
