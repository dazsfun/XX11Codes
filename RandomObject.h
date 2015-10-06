#pragma once
#include <opencv2\opencv.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include <iostream>
#include "SatOrbitNew.h"
#include <stdlib.h> 
#include <time.h>  

using namespace cv;
using namespace std;

class RandomObject
{
private:
	int width_;
	int height_;

public:
	RandomObject(int _width, int _height)
	{
		width_ = _width;
		height_ = _height;
		_rightx = new double[_width * _height];
		_righty = new double[_width * _height];
		for (int j = 0; j < height_; j++)
		{
			for (int i = 0; i < width_; i++)
			{
				_rightx[j * width_ + i] = i;
			}
		}

		for (int j = 0; j < height_; j++)
		{
			for (int i = 0; i < width_; i++)
			{
				_righty[j * width_ + i] = j;
			}
		}

		cout << "是否加入随机误差" << endl;
		int heish0;
		cin >> heish0;
		if (heish0)
		{
			for (int j = 0; j < height_; j++)
			{
				for (int i = 0; i < width_; i++)
				{
					time_t t;
					srand((unsigned)i * 10000 + j * 1000 + (unsigned)time(&t));

					_rightx[j * width_ + i] += rand()/(RAND_MAX + 1.0)* 6.7 *  1.3 * ((8500 - fabs((double)j - 8500)) / 8500 + 0.4);

					srand((unsigned)i * 10000 + j * 1000 + (unsigned)time(&t));
					_righty[j * width_ + i] += rand() / (RAND_MAX + 1.0) * 1.5;

				}
			}

		}

		cout << "是否拉伸高度" << endl;
		int heish;
		cin >> heish;
		if (heish)
		{
			HeightStretch("hei1", 44);
		}
		int heish2;
		cin >> heish2;
		if (heish2)
		{
			WidthStretch("he2", 25);
		}


		int heish3;

		cout << "是否加仿射变换" << endl;
		cin >> heish3;
		double angle = 0.2 * CV_PI / 180.0;
		if (heish3)
		{ 
			double rotm[9] = { cos(angle), -sin(angle), 15.25
				, sin(angle), cos(angle), 94.27,
				0, 0, 1 };

			for (int j = 0; j < height_; j++)
			{
				for (int i = 0; i < width_; i++)
				{
					double hei[3] = { _rightx[j * width_ + i], _righty[j * width_ + i], 1 };
					double info[3];
					CSatOrbit helper;
					helper.mult(rotm, hei, info, 3, 3, 1);
					_rightx[j * width_ + i] = info[0];
					_righty[j * width_ + i] = info[1];
				}
			}
		}

		ifstream reademydata("c:\\Temp\\match\\heimatch.xls");
		ofstream reademydata1("c:\\Temp\\match\\heimatch1.xls");

		double a[9];
		char readernow[1024];
		
		reademydata1.precision(20);
		do
		{
			reademydata.getline(readernow, 1024);
			sscanf_s(readernow, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &a[0], &a[1], &a[2], &a[3], &a[4], &a[5],
				&a[6], &a[7], &a[8]);

			int indexX = a[0];
			int indexY = a[1];
			double xvalue = _rightx[indexY * width_ + indexX] - a[0];
			double yvalue = _righty[indexY * width_ + indexX] - a[1];
			reademydata1 << a[0] << "\t" << a[1] << "\t" << xvalue << "\t" << yvalue << "\t" << "\t" << _rightx[indexY * width_ + indexX] << "\t" << _righty[indexY * width_ + indexX] << endl;
		} while (reademydata.peek() != EOF);

		reademydata.close();
		reademydata1.close();

	}

	void HeightStretch(string curveinfo, int startindex)
	{
		ifstream curve("c:\\Temp\\match\\mycurve.xls");
		int ncouter = 0;
		for (int i = startindex; i < height_; i++)
		{
			char reader[1024];
			double a[5];
			curve.getline(reader, 1024);
			sscanf_s(reader, "%lf %lf %lf %lf %lf", &a[0], &a[1], &a[2], &a[3], &a[4]);
			if (i == startindex) continue;
			for (int j = 0; j < width_; j++)
			{
				_righty[i * width_ + j] = a[2];
			}

			if (curve.peek() == EOF)
			{
				break;
			}
		}

	}

	void WidthStretch(string curveinfo, int startindex)
	{
		ifstream curve("c:\\Temp\\match\\mycurve_x.xls");
		int ncouter = 0;
		for (int i = startindex; i < width_; i++)
		{
			char reader[1024];
			double a[5];
			curve.getline(reader, 1024);
			sscanf_s(reader, "%lf %lf %lf %lf %lf", &a[0], &a[1], &a[2], &a[3], &a[4]);
			if (i == startindex) continue;
			for (int j = 0; j < height_; j++)
			{
				_rightx[i * width_ + j] = a[2];
			}

			if (curve.peek() == EOF)
			{
				break;
			}
		}
	}

	Mat rightx;
	Mat righty;

	double* _rightx;
	double* _righty;



	~RandomObject();
};

