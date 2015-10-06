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
	Enum_Mat_Fill_Symmetric, ///���Ļ���
	Enum_Mat_Fill_Sequence�� ///˳��
};

///���Namespace��Ҫ������OPENCV���㷨�ĸĽ�
namespace safeopencv
{
	//From minmaxloc to MinMaxLoc_Subpixel
	///��;: ƥ��
	///ʹ��Ҫ��: ��������srcΪ ��ͨ�� float��(��Ϊfloat���Ϳ϶������),ʹ��ǰ��������һ��ת������, �����ڲ����벻����ʽת��
	///         ��������LocInΪ�����Сֵ�ĳ�ֵ,���������趨,Ҳ����ʹ��MinMaxLoc�ļ����� 
	///			������SubPixLoc,�õ��������Сֵλ��
	double  MinMaxloc_subpixel(CV_OUT cv::Point2d* SubPixLoc,
		cv::InputArray src,
		CV_IN_OUT cv::Point* LocIn,
		CV_IN_OUT const int Method
		);

	///��;�����ߵĸ߾�����ϣ�ֱ��ͼƥ�䣨����У������һάƥ��
	///������ opencv�� alglib
	///�ҵ�һ��������ѣ�������ɢ�㣩�ķ�ֵ����ֵ
	///˵���� ��������Ϊdouble�͵㣬��Ϊ�����������������
	///������������MinMaxloc_subpixel�����������1ά�ģ�����Ҫ��ֵ
	double Curve_MinMax(vector<Point2d> curveinfo, double& maxvalue, double& minvalue, double& maxloc, double& minloc);
	bool Inside(Rect parent, Rect child);

	///����ģ�����
	///
	void Mat_Fill(vector<Point2f>& srcMat_,int windowsize,Mat_Fill_Type filltype);
	void Pixel_Warp_Affine(const vector<Point2f>& srcMat, vector<Point2f>& dstMat,const Point2f& pixestart, double*  affine_data);

	void Rot2Eulor(const double* rotvalue, double& rot1, double& rot2, double& rot3);

	void Eulor2Rot(const double& rot1, const double& rot2, const double& rot3, double* rotvalue);

	const double Resample( Mat& srcMat, Point2f pixel_loc, ResampleKernels* resmaplekernel);

	void Match_Curve( vector<Point2d>& matchleft,  vector<Point2d>& matchright,string img_path_left,string img_path_right,barycentricinterpolant& left2right_x,barycentricinterpolant& left2right_y);

	
};

