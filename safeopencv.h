#pragma once
#include <opencv2\opencv.hpp>
#include "ap.h"
#include "interpolation.h"
#include "log.h"
#include "ResampleKernels.h"
using namespace alglib;
using namespace cv;

enum Mat_Fill_Type
{
	Enum_Mat_Fill_Symmetric, ///中心化，
	Enum_Mat_Fill_Sequence， ///顺序
};

///这个Namespace主要包含对OPENCV中算法的改进
namespace safeopencv
{
	//From minmaxloc to MinMaxLoc_Subpixel
	///用途: 匹配
	///使用要求: 输入数据src为 单通道 float型(不为float类型肯定会出错),使用前对数据作一次转换即可, 函数内部代码不做格式转换
	///         输入数据LocIn为最大最小值的初值,可以任意设定,也可以使用MinMaxLoc的计算结果 
	///			结果输出SubPixLoc,得到的最大最小值位置
	double  MinMaxloc_subpixel(CV_OUT cv::Point2d* SubPixLoc,
		cv::InputArray src,
		CV_IN_OUT cv::Point* LocIn,
		CV_IN_OUT const int Method
		);

	///用途：曲线的高精度拟合，直方图匹配（辐射校正），一维匹配
	///依赖： opencv， alglib
	///找到一条曲线最佳（而非离散点）的峰值，谷值
	///说明： 输入数据为double型点，因为这个是最方便的数据类型
	///功能上类似于MinMaxloc_subpixel，但是是针对1维的，不需要初值
	double Curve_MinMax(vector<Point2d> curveinfo, double& maxvalue, double& minvalue, double& maxloc, double& minloc);
	bool Inside(Rect parent, Rect child);

	///像素模版填充
	///
	void Mat_Fill(vector<Point2f>& srcMat_,int windowsize,Mat_Fill_Type filltype);
	void Pixel_Warp_Affine(const vector<Point2f>& srcMat, vector<Point2f>& dstMat,const Point2f& pixestart, double*  affine_data);

	void Rot2Eulor(const double* rotvalue, double& rot1, double& rot2, double& rot3);

	void Eulor2Rot(const double& rot1, const double& rot2, const double& rot3, double* rotvalue);

	const double Resample( Mat& srcMat, Point2f pixel_loc, ResampleKernels* resmaplekernel);

	void Match_Curve( vector<Point2d>& matchleft,  vector<Point2d>& matchright,string img_path_left,string img_path_right,barycentricinterpolant& left2right_x,barycentricinterpolant& left2right_y);

	
};

