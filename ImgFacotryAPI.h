#pragma once 
#ifndef _CLASS_HEADER_IMGFACTORY_API_HEADER
#define _CLASS_HEADER_IMGFACTORY_API_HEADER
#include <gdal_priv.h>
#include <ogrsf_frmts.h> 
#include <opencv2\opencv.hpp>
#include <opencv2\features2d\features2d.hpp>
#include <fftw3.h>
#include "safeopencv.h"
#include "ResampleKernels.h"
#include <omp.h>

using namespace alglib;

/*
Author:Fang Chen. Version: T0.01
2015/02/19
*/

/*
影像处理类
将使用两个开源库
OPENCV和GDAL
*/#include <opencv2\core\core.hpp>
#include <opencv2\features2d\features2d.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2/video/tracking.hpp>
#include <set>


#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2\nonfree\features2d.hpp>
#include <gdal_priv.h>
#include <string>
#include "log.h"


using namespace cv;
using namespace std;
#include "VirtualMachine.h"


//#include "RCProcess.h"
enum MatchFactoryType
{
	Enum_Match_RigurousModel,   //影像本身具有严密模型,从指针获取数据,或者反序列化
	Enum_Match_GEOTIFF,  /// 
	Enum_Match_IMGHEADER,
};

typedef struct MatchResult
{
	double _lx;
	double _ly;
	double _rx;
	double _ry;
	double _lpx;
	double _lpy;
	double _rpx;
	double _rpy;

	double corre;
	int type;
};


typedef struct MatchToolType
{
	string filepath;
	void* _ptool;
	string outputpath;
	bool xy2ground;
	bool SwitchDirction()
	{
		xy2ground = !xy2ground;
	}

	void Translate(bool fromxy2ground, Point3d* container_result, int size_calc, int index_id,bool mosaic_mode = false)
	{
		if (fromxy2ground)
		{
			for (int i = 0; i < size_calc; i++)
			{
				double loninfo, latinfo;
				static_cast<VirtualMachine*>(_ptool)->FromXY2LonLatTest(container_result[i].x, container_result[i].y, 100,
					loninfo, latinfo, index_id/*static_cast<VirtualMachine*>(_ptool)->Current()*/);
				container_result[i].x = loninfo;
				container_result[i].y = latinfo;
				container_result[i].z = 100;
			}
		}
		else
		{
			int overlap_count = 0;
			for (int i = 0; i < size_calc; i++)
			{
				//cout << "i: " << i << "\\" << overlap_count << "\\" << size_calc << "\r";
				double xinfo,  yinfo;
				static_cast<VirtualMachine*>(_ptool)->FromlatlonH2xy(container_result[i].x, container_result[i].y,
					 100, xinfo, yinfo, index_id);
				if (mosaic_mode)
				{
					double xinfo_2, yinfo_2;
					 static_cast<VirtualMachine*>(_ptool)->FromlatlonH2xy(container_result[i].x, container_result[i].y,
						100, xinfo_2, yinfo_2, index_id + 1);
					

				
					container_result[i].x = (fabs(xinfo - 479) < fabs(xinfo_2)) ? xinfo_2 : xinfo;
					container_result[i].y = (fabs(xinfo - 479) < fabs(xinfo_2)) ? yinfo_2 : yinfo;
					container_result[i].z = (fabs(xinfo - 479) < fabs(xinfo_2)) ? 1 : 0;

				}
				else
				{
					container_result[i].x = xinfo;
					container_result[i].y = yinfo;
					container_result[i].z = 0;
					continue;
				}
			}
		}
	}

	void FromXY2LONLAT(const double& imgx, const double& imgy, const double& height, double& loninfo, double& latinfo) const
	{
		if (_typeofFac == Enum_Match_RigurousModel)
		{
			static_cast<VirtualMachine*>(_ptool)->FromXY2LonLatTest(imgx, imgy, height, loninfo, latinfo, 67/*static_cast<VirtualMachine*>(_ptool)->Current()*/);
		}
		else if (xy2ground)
		{
			double transforminfo[6];
			memcpy(transforminfo, (static_cast<double*>(_ptool)), sizeof(double)* 6);
			loninfo = transforminfo[0] + imgx *transforminfo[1] + imgy * transforminfo[2];
			latinfo = transforminfo[3] + imgx * transforminfo[4] + imgy * transforminfo[5];
		}
	};
	void FromXY2Rect(const double& imgx, const double& imgy, const double& height, double& xinfo, double& yinfo, double& zinfo)
	{
		double loninfo, latinfo;
		FromXY2LONLAT(imgx, imgy, height, loninfo, latinfo);
		static_cast<VirtualMachine*>(_ptool)->Geo2Rect(loninfo, latinfo, height, xinfo, yinfo, zinfo);
	}
	void FromRect2XY(const double& lon, const double& lat, const double& height, double& imgx, double& imgy) const
	{
		double geoinfo[6];
		memcpy(geoinfo, static_cast<double*>(_ptool), sizeof(double)* 6);
		imgx = geoinfo[0] * lon + geoinfo[1] * lat + geoinfo[2];
		imgy = geoinfo[3] * lon + geoinfo[4] * lat + geoinfo[5];
	}
	MatchFactoryType _typeofFac;
}ImgFactoryExtend;

//#pragma comment(lib,"D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\RCProcess.lib")
using namespace std;

#pragma comment(lib,"gdal_i.lib")
class ImgFacotryAPI
{
private:
	int heivaluecount;
	unsigned int _width;
	unsigned int _height;
	unsigned int _singlebyte;
	string _path;
	unsigned int* _imgcontent;
	OGRCoordinateTransformation * poTransform;
	void DrawHist(const MatND& b_hist);
	bool Histo_Match_Image(MatchToolType& leftinfi, MatchToolType& rightinfi);
	void Translate(const MatchToolType& _lefttool, const MatchToolType& _rightTool, const double& height, double* transx, double* transy, int sizeall);
private:
	void Match_Pixel(MatchToolType& _lefttool, MatchToolType& _rightool);
	void Match_Phase_Correlation(MatchToolType& _left, MatchToolType& _right, vector<DMatch>);
	double CalcCorrelation(double * pl, double *pr, int wm, int wn);
	vector<Point2d> _leftkey_;
	vector<Point2d> _rightkey_;
	vector<TimeCalibStruct> p_mycalib;
	vector<Point2d> _heitemp;
	vector<Point2d> _heitemp_x;

	barycentricinterpolant x2xinfo;
	barycentricinterpolant y2yinfo;


	void SetCoordinates()
	{
		OGRSpatialReference oSRS, *poLatLong;
		oSRS.SetProjCS("UTM30(WGS84) IN NORTH");
		oSRS.SetWellKnownGeogCS("WGS84");
		oSRS.SetUTM(30, TRUE);
		poLatLong = oSRS.CloneGeogCS();
		poTransform = OGRCreateCoordinateTransformation(poLatLong, &oSRS);
	}

public:
	void GetNeightContent(const vector<Point2d>& LEFTMATCH, const vector<Point2d>& RIGHTMATCH)
	{
		_leftkey_.clear();
		_rightkey_.clear();
		if (LEFTMATCH.size() != RIGHTMATCH.size())
		{
			cout << "邻近景匹配点个数不一致" << endl;
			return;
		}

		else
		{
			for (int i = 0; i < LEFTMATCH.size(); i++)
			{
				_leftkey_.push_back(LEFTMATCH[i]);
				_rightkey_.push_back(RIGHTMATCH[i]);
			}
		}
	}

	void GetClibContent(vector<Point2d>& container_, bool bxarray = false)
	{
		if (!bxarray)
		{
			for (int i = 0; i < _heitemp.size(); i++)
			{
				container_.push_back(_heitemp[i]);
			}
		}
		else
		{
			for (int i = 0; i < _heitemp.size(); i++)
			{
				container_.push_back(_heitemp_x[i]);
			}
		}
	}
	void Begin_Match_WorkFlow(string srcImg, string refimg, VirtualMachine* info);
	ImgFacotryAPI(unsigned int width, unsigned int height, unsigned int singlebyte, string path);

	void Mosaic(vector<string> imagenames, int width, int height, bool align);
	~ImgFacotryAPI();

	bool  StitchFile(vector<string> iamgenames, string output, int eachwidth, int eachheight);
	void LonLat2MapXY(const double& lon, const double& lat, const double& height, double& mapx_, double& mapy_);
public:
	void MatchAndCompare(string rootdir, Mat* heimat);
	template<typename T>void DetectorHistoMatchCalibration(string imagename, int singleCCDwidth, int singleccdHeight, int singlebytenumber, GDALDataType etype)
	{
		SplitByCol(imagename, singleCCDwidth);  //Classify image by detector id
		Mat src;
		src = imread(imagename, CV_LOAD_IMAGE_ANYDEPTH);

		int bins = static_cast<int>(pow(2, singlebytenumber));
		int hist_size[] = { bins };
		float range[] = { 0, static_cast<float>(bins) };
		const float* ranges[] = { range };
		int wholewidth = src.cols;
		int wholehei = src.rows;
		int sizecal = wholewidth / singleCCDwidth;

		MatND hiswhole; int channels[] = { 0 };


		calcHist(&src, 1, channels, Mat(), hiswhole, 1, hist_size, ranges, true, false);

		cout << "begin histo matching" << endl;
		vector<float*> hiscomponents;
		cout << endl;
		hiscomponents.resize(singleCCDwidth);
		//#pragma omp parallel for
#pragma omp parallel for num_threads(4)
		for (int i = 0; i < singleCCDwidth; i++)
		{
			//获取单探元的影像
			///得到单探元的查找表
			float* looktable = new float[bins];
			char reader[1024];
			sprintf_s(reader, "%d_new.tif", i);


			HistoMatching(imagename + static_cast<string>(reader), hiswhole, singlebytenumber, looktable);
			hiscomponents[i] = (looktable);

			ofstream ofstemp;
			char readfuck[1024];
			sprintf_s(readfuck, "counter_%d.txt", i);
			ofstemp.open(readfuck, ios::out);
			ofstemp << i << endl;
			ofstemp.close();
		}

		cout << "Histo matching finished" << endl;
		///永远 ，不要用OPENCV读取影像，尤其是大块的
		float* infobytes = new float[wholewidth * wholehei];
		unsigned short* infobyteshortint = new unsigned short[wholewidth * wholehei];
		memset(infobytes, 0, sizeof(float)* wholewidth * wholehei);
		memset(infobyteshortint, 0, sizeof(unsigned short)* wholewidth * wholehei);
		GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");


		GDALDataset* poDataOri = (GDALDataset*)GDALOpen(imagename.c_str(), GA_ReadOnly);
		GDALRasterBand* poBandOri = poDataOri->GetRasterBand(1);

		poBandOri->RasterIO(GF_Read, 0, 0, wholewidth, wholehei, infobyteshortint, wholewidth, wholehei, GDT_UInt16, 0, 0);
		cout << "got image" << endl;
		for (int j = 0; j < wholehei; j++)
		{
			for (int i = 0; i < wholewidth; i++)
			{
				infobytes[j * wholewidth + i] = static_cast<float>(infobyteshortint[j * wholewidth + i]);
			}
		}
		delete[] infobyteshortint;

		cout << "begin allocation" << endl;
		for (int j = 0; j < wholehei; j++)
		for (int i = 0; i < wholewidth; i++)
		{
			int heitempint = static_cast<int>(infobytes[j * wholewidth + i]);
			heitempint = (heitempint > bins - 1) ? bins - 1 : heitempint;
			infobytes[j * wholewidth + i] = hiscomponents[i%singleCCDwidth][heitempint];
		}

		cout << "allocation finished" << endl;
		vector<string>  mosaicname;
		for (int i = 0; i < sizecal; i++)
		{
			char readernow[1024];
			float* infweget = new float[singleCCDwidth * wholehei];
			//unsigned char* infwegetchar = new unsigned char[singleCCDwidth * wholehei];

			sprintf_s(readernow, "%d_.tif", i);

			GDALDataset* poDataset = poDriver->Create((imagename + static_cast<string>(readernow)).c_str(), singleCCDwidth, wholehei, 1, GDT_Float32, NULL);
			GDALRasterBand* poBand = poDataset->GetRasterBand(1);
			for (int j = 0; j < wholehei; j++)
			{
				memcpy(infweget + j * singleCCDwidth, infobytes + j * wholewidth + i * singleCCDwidth, sizeof(float)* singleCCDwidth);
			}




			poBand->RasterIO(GF_Write, 0, 0, singleCCDwidth, wholehei, infweget, singleCCDwidth, wholehei, GDT_Float32, 0, 0);
			GDALClose(poDataset);
			delete[] infweget;
			//			delete[] infwegetchar;

			//delete[] infweget;
			//delete[] infwegetchar;
			cout << i << "/" << singleCCDwidth << "\r";
			mosaicname.push_back(imagename + static_cast<string>(readernow));
		}


		GDALClose(poDataOri);

		delete[] infobytes;
		for (int i = 0; i < hiscomponents.size(); i++)
		{
			if (hiscomponents[i] != nullptr)
			{
				delete[] hiscomponents[i];
				hiscomponents[i] = nullptr;
			}
		}

		//StitchFile(mosaicname, "whatevernew.tif", 480, 10786);
	};
	unsigned int* IRow(unsigned int rowindex);
	bool IOutPut(bool radiometrix = false)
	{
		GDALDataType dataType;
		if (_singlebyte <= 8)
		{
			dataType = GDT_Byte;
		}

		else if (_singlebyte > 8 && _singlebyte <= 16)
		{
			dataType = GDT_UInt32;
		}

		else if (_singlebyte > 16 && _singlebyte <= 32)
		{
			dataType = GDT_UInt32;
		}

		else if (_singlebyte > 32 && _singlebyte <= 64)
		{
			dataType = GDT_Float64;
		}

		else
		{
			return  false;
		}

		GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset* poDataset = poDriver->Create(_path.c_str(), _width, _height, 1, dataType, NULL);

		GDALRasterBand* poBand = poDataset->GetRasterBand(1);

		poBand->RasterIO(GF_Write, 0, 0, _width, _height, _imgcontent, _width, _height, dataType, 0, 0);

		GDALClose(poDataset);

		//bool __declspec(dllexport) X6RCProcessed(const char *str_FileName, const char *str_RCType, int OverLapPixelNum[]);
		if (radiometrix)
		{
			int  overlap[5] = { 0 };
			//X6RCProcessed(_path.c_str(), "Histo", overlap);
		}
		//	delete poBand;
		//	delete poDriver;
	}
	template<typename T>bool Mosaic(vector<string> images, vector<Point2f> positions, int eachwidth, int eachheight)
	{
		int maxwidth, maxheight;
		bool changealigh = false;
		int offset = fabs(positions[positions.size() - 1].y) + eachheight - 1;
		if (positions[positions.size() - 1].y < 0)
		{
			changealigh = true;
			maxheight = eachheight + fabs(positions[positions.size() - 1].y) + eachheight;
			maxwidth = positions[positions.size() - 1].x + eachwidth;
		}

		else
		{
			maxwidth = positions[positions.size() - 1].x + eachwidth;
			maxheight = positions[positions.size() - 1].y + eachheight;
		}

		GDALDataset* poDataset1 = (GDALDataset*)GDALOpen(images[0].c_str(), GA_ReadOnly);
		GDALDataType typeofdata = poDataset1->GetRasterBand(1)->GetRasterDataType();
		GDALClose(poDataset1);
		GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset* poDataset = poDriver->Create("whateverits.tif", maxwidth, maxheight, 1, typeofdata, NULL);
		GDALRasterBand* poBand = poDataset->GetRasterBand(1);

		for (int i = 0; i < images.size(); i++)
		{
			GDALDataset* poDatasetTemp = (GDALDataset*)GDALOpen(images[i].c_str(), GA_Update);
			GDALRasterBand* poBandTemp = poDatasetTemp->GetRasterBand(1);
			assert(poDatasetTemp->GetRasterBand(1)->GetRasterDataType() == typeofdata);
			unsigned short* infoarea = new unsigned short[eachheight * eachwidth];
			memset(infoarea, 0, sizeof(unsigned short)*eachheight * eachwidth);
			poBandTemp->RasterIO(GF_Read, 0, 0, eachwidth, eachheight, infoarea, eachwidth, eachheight, typeofdata, 0, 0);
			poBand->RasterIO(GF_Write, positions[i].x, offset + positions[i].y, eachwidth, eachheight, infoarea, eachwidth, eachheight, typeofdata, 0, 0);
			GDALClose(poDatasetTemp);
			delete[] infoarea;
		}

		GDALClose(poDataset);
		return true;
	}

public:
	void ImageScale(string path, double topx, double topy, double width, double height, double scalex, double scaley);
	bool PreMatch(MatchToolType& _lefttool, MatchToolType& _right, VirtualMachine* hei = nullptr, VirtualMachine* heiright = nullptr);
	bool PreMatch(MatchToolType& _lefttool, MatchToolType& _right, VirtualMachine* hei = nullptr, int nonsense = 0);
	bool Phase_Correction(string _leftpath, string rightpath);

	///在进行相位匹配之前，请务必保证初值误差在4个像素以内，可以先使用特征匹配，或者严密模型，使得初值误差达到像素级
	///相位匹配将能提供尽可能密集的匹配结果，用于内部畸变校正
	bool Phase_Correction(Mat matleft, Mat matright, double& maxvalue, Point2d& pointmove);

	///直接用成像模型拼接影像
	bool Mosaic_RigurouModel_Edge(MatchToolType& _lefttool,MatchToolType& _right,VirtualMachine* hei,int index_id);

	///Left and right should all be double
	bool Match_Phase_Correction(MatchToolType& _lefttool, MatchToolType& _righttool, int windowsize, int windowmove, double thresold);



	bool LSM_Match(MatchToolType leftinfi, MatchToolType rightinfi, int windowsize);

	///Template matching
	bool Template_PreMatch(string match_path);
	bool Match_Template(MatchToolType leftinfi, MatchToolType rightinfi, int windowsize);
	bool Form_Template();

	bool lsm_multi_times;
private:
	bool LSM_Match(MatchToolType leftinfi, MatchToolType rightinfi, const vector<Point2d>& leftori, const vector<Point2d>& rightori,
		vector<Point2d>& _leftnewkey, vector<Point2d>& _rightnewkey,
		vector<MatchResult>& _lsmresult, int windowsize);
	bool is_shifted;
	void HistoMatching(string matinfo, MatND hisinfo, int singlebytenumber, float* T);
	void HistoMatching(Mat* info, MatND hisinfo, int singlebytenumber, unsigned short* T);
	void SplitByCol(string imagename, int singleccdbyte);
	void MatchAndCompare(string _left, string _right, cv::Ptr<cv::Mat> matinfo, cv::Rect rectleft, cv::Rect rectright);
};
#endif
