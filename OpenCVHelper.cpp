#include "stdafx.h"
#include "OpenCVHelper.h"
#include <opencv2\calib3d\calib3d.hpp>

OpenCVHelper::OpenCVHelper()
{
}


OpenCVHelper::~OpenCVHelper()
{
}


bool OpenCVHelper::FindRT(vector<double*> pvalue, vector<double*> qvalue, double* pmiddle, double* qmidde, double* rot, double* trans)
{
	if (pvalue.size() != qvalue.size())
	{
		cout << "点集大小不一致" << endl;
		return false;
	}


	Mat middleq(3, 1, CV_64F, &qmidde);
	Mat middlep(3, 1, CV_64F, &pmiddle);
	/*
	std::vector<Point3f> _pvector;
	std::vector<Point3f> _qvector;
	for (int i = 0; i < pvalue.size(); i++)
	{
	double a[3][1];
	Point3d tempvecp(pvalue[i][0], pvalue[i][1], pvalue[i][2]);
	Point3d tempvecq(qvalue[i][0], qvalue[i][1], qvalue[i][2]);
	_pvector.push_back(tempvecp);
	_qvector.push_back(tempvecq);
	}

	std::vector<Point3f> inlier;
	inlier.resize(10);
	inlier.reserve(10);
	double rottest[3][4] = { 0 };
	Mat H;
	//estimateAffine3D(_pvector, _qvector, H, inlier);
	//memcpy(rot, rottest, sizeof(double)* 9);
	double checkpp[3] = { pvalue[0][0],pvalue[0][1],pvalue[0][2] };
	double checkqq[3] = { qvalue[0][0], qvalue[0][1], qvalue[0][2] };
	Mat checkmatpp(1, 3, CV_64F, &checkpp);
	Mat checkmatqq(1, 3, CV_64F, &checkqq);
	Mat calcmatqq(1, 3, CV_64F,&trans);
	calcmatqq = checkmatpp * H - checkmatqq;
	*/

	Mat Xmatrix(3, pvalue.size(), CV_64F);
	Mat YmatrixT(pvalue.size(), 3, CV_64F);
	for (int i = 0; i < pvalue.size(); i++)
	{
		Mat mx = Xmatrix.col(i); Mat mxnew(3, 1, CV_64F);
		for (int j = 0; j < 3; j++)
			mxnew.at<double>(j, 0) = pvalue[i][j];
		Mat my = YmatrixT.row(i); Mat mynew(1, 3, CV_64F);
		for (int j = 0; j < 3; j++)
			mynew.at<double>(0, j) = qvalue[i][j];
		mxnew.copyTo(mx);
		mynew.copyTo(my);
	}

	Mat_<double> S = Xmatrix * YmatrixT;
	SVD svd(S);
	Mat svd_u = svd.u;
	Mat svd_vt = svd.vt;
	Mat svd_w = svd.w;
	Matx33d W(0, -1, 0, 1, 0, 0, 0, 0, 1);//HZ 9.13
	Mat_<double> R = svd_u * Mat(W) * svd_vt; //
	double rottest[3][3] = { 0 };
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			rottest[i][j] = R[i][j];
	}
	memcpy(rot, rottest, sizeof(double)* 9);
	Mat_<double> T = svd_u.col(2); //u3
	for (int i = 0; i < 3; i++)
	{
		trans[i] = T[i][0];
	}
	double det = determinant(svd_vt.t() * svd_u.t());
	double kk[3][3] = { 1, 0, 0, 0, 1, 0, 0, 0, det };//HZ 9.13
	Mat WNew(3, 3, CV_64F, &kk);
	Mat_<double> Rnew = svd_vt.t() * WNew * svd_u.t();
	Mat hei = middleq - (svd_vt.t() * WNew * svd_u.t()).inv() * middlep;
	Mat_<double> RnewInv = Rnew.inv();
	Vec3d heit = hei.col(0);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			rot[i * 3 + j] = RnewInv[i][j];
		}
	}

	for (int i = 0; i < 3; i++)
	{
		trans[i] = hei.at<double>(i, 0);
	}

}