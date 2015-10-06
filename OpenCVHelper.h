#include <opencv2\opencv.hpp>
#pragma once
#ifndef _MY_DLL_MANAGE_
#endif
#include <iostream>
using namespace std;
using namespace cv;
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define random() (static_cast<double>(rand()%1000) / (1000))
class OpenCVHelper
{
public:
	OpenCVHelper();

public:
	bool SelfCheck()
	{
		double roll = 10;
		double pitch = 6;
		double yaw = 10;
		double rotvecvalue[3] = { roll / 180 * CV_PI, pitch / 180 * CV_PI, yaw / 180 * CV_PI };
		Mat rotatevector(3,1,CV_64F,&rotvecvalue);

		double rotvalue[9];
		CvMat value1 = rotatevector;
		Mat rotationMatrix(3, 3, CV_64F,&rotvalue);
		CvMat value2 = rotationMatrix;
		cvRodrigues2(&value1, &value2);
		double avgP[3] = { 0 }; double  avgQ[3] = {0};
		vector<double*> pointP, pointQ;
		for (int x = 0; x < 10; x++)
		{

			double* tempvalue1 = new double[3];
			double* tempvalue2 = new double[3];
			srand((int)time(0)*1000000);
			tempvalue1[0] =  random();
			//srand((int)time(0));
			srand((int)time(0)*100);
			tempvalue1[1] = (1 - tempvalue1[0]) * random();
			tempvalue1[2] = sqrt(1 - pow(tempvalue1[0], 2) - pow(tempvalue1[1], 2));
			Mat vecp(1,3,CV_64F,tempvalue1);
			Mat resulthere(1, 3, CV_64F,tempvalue2);
			resulthere = vecp * rotationMatrix;
			for (int i = 0; i < 3; i++)
			{
				avgP[i] += tempvalue1[i] / 10;
				tempvalue2[i] += 0.01;
		//		if(i%2==0)tempvalue2[i] +=  0.03;
			//	else tempvalue2[i] -= 0.03;
				avgQ[i] += tempvalue2[i] / 10;
			}
			pointP.push_back(tempvalue1);
			pointQ.push_back(tempvalue2);
		}

		for (int i = 0; i < 10; i++)
		{
			for (int x = 0; x < 3; x++)
			{
				pointP[i][x] -= avgP[x];
				pointQ[i][x] -= avgQ[x];
			}
		}

		double rotCaled[9] = { 0 };
		double trans[3] = { 0 };
		FindRT(pointP, pointQ, avgP, avgQ, rotCaled, trans);
		Mat rotresult(3, 3, CV_64F, &rotCaled);
		CvMat result1 = rotresult;
		double rollpitchyaw[3] = { 0 };
		Mat rotvecresult(3, 1, CV_64F, &rollpitchyaw);
		CvMat result2 = rotvecresult;
		cvRodrigues2(&result1, &result2);
		cout.precision(20);
		cout << rollpitchyaw[0] * 180 / CV_PI << "\t" << rollpitchyaw[1] * 180 / CV_PI << "\t" << rollpitchyaw[0] * 180 / CV_PI<<endl;
		cout << trans[0] << "\t" << trans[1] << "\t" << trans[2] << endl;
		for (int i = 0; i < 9; i++)
		{
			cout << rotCaled[i] << endl;
			cout << rotvalue[i] << endl;
		}
		return true;
	}

	bool FindRT(vector<double*> pvalue, vector<double*> qvalue, double* pmiddle, double* qmidde, double* rot, double* trans);


	~OpenCVHelper();
};

