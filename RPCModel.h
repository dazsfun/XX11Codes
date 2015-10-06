// RPCModel.h: interface for the RPCModel class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RPCMODEL_H__91231331_BBBB_4527_B780_D8F0BEA9DDBB__INCLUDED_)
#define AFX_RPCMODEL_H__91231331_BBBB_4527_B780_D8F0BEA9DDBB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "RSGeometryModel.h"
#include "Curve.h"
#include "GeoRSDataStruct.h"
#include "GridInterpolate.h"

static ImageUnchangeModel ImageUnchangeSinglelon;
static ObjectUnchangeModel ObjectUnchangeSinglelon;

/************************************************************************/
/* RPC(rational polynomial coefficients) Model                         */
/* Derived from the RSGeometryModel                                     */
/* FromGroundToImage                                                    */
/* FromImageToGround && GetModelRange && GetGolbalAffineModel           */
/* should be redefined.                                                 */
/* There are two ways to initialize the model: read rpc file or         */
/* computer the rpc from other RSGeometryModel                          */
/*                                                                      */
/*                                                         _            */
/*  Author: 潘红播                                   |_|  |_)           */
/*  Email:  hongbo640@gmail.com                      | |  |_)           */
/*  Edit time: 2010 3 04                                                */
/*  Last modified: 2011 2 27                                            */
/************************************************************************/

/************************************************************************/
/* actually, the rpc model usually provided the the corresponding from  */
/* ground to image space, so should we provide a switcher to change the */
/* corresponding.                                                       */
/*  add the function of compute the coordinate from image to ground     */
/*                                                           2011 12 27 */
/************************************************************************/

/* another question: if it's from (x,y,h) -> (lat, lon)                  */
/* so a switcher is needed in this section                               */
/* four part of this class should modified                               */
/* init, creat, from ground to image and from image to ground            */


// copy construction and copy assignment don't contain the information of solving the RPC parameters
#include "VirtualMachine.h"
#include <opencv2\opencv.hpp>
class VirtualMachine;
class  RSRPCModel : public RSGeometryModel
{

public:
	double MaxHeight()
	{
		return HEIGHT_OFF + HEIGHT_SCALE/2;
	}

	double MinHeight()
	{
		return HEIGHT_OFF - HEIGHT_SCALE/2;
	}

	double Height()
	{
		return HEIGHT_OFF;
	}

	double ImgHeight()
	{
		return SAMP_SCALE;
	}

	double ImgWidth()
	{
		return LINE_SCALE;
	}

private:

/*************************************************************************/
	double LNUM[20],LDEN[20],SNUM[20],SDEN[20];
	double LAT_OFF,LAT_SCALE,LONG_OFF,LONG_SCALE,HEIGHT_OFF,HEIGHT_SCALE;
	double LINE_OFF,LINE_SCALE,SAMP_OFF,SAMP_SCALE;

	// affine transform parameter
	double px[3],py[3];
	double invpx[3],invpy[3];
/*************************************************************************/

	//RPC扩展模型
	//##ModelId=4E4A89B20283
	RPCGCPSTRUCT *m_RPCgcp;	
	//##ModelId=4E4A89B20292
	long m_nRPCgcpNum;

	bool m_bflag;

	//RPC 控制点
	double *m_Y,*m_X;
	double *m_P,*m_L,*m_H;

	// the range of the model
	//##ModelId=4E4A89B20294
	DPolygon m_ModelRange;

	// 计算RPC模型
	double   m_MaxElevation, m_MinElevation;
	//##ModelId=4E4A89B20299
	DRect    m_ImgRange;
	int		 m_nx,m_ny,m_nz;	
	//##ModelId=4E4A89B2029F
	const void*    m_pImageTrans;
	//##ModelId=4E4A89B202A0
	const void*    m_pGroundTrans;
	//##ModelId=4E4A89B202AF
	const void*    m_pOrignalModel;

	// for accelerate the backward transformation, and it should be initialized from the init function
	// and it should needed in the method of reading from the RPCfile
	//##ModelId=4E4A89B202B1
	AffineTranslationModel m_GlobalAffineModel;

public:
	//##ModelId=4E4A89B202BF
	RSRPCModel();
	//##ModelId=4E4A89B202DE
	virtual ~RSRPCModel();

	void GetHeightOff(double& heightoffer)
	{
		heightoffer = HEIGHT_OFF;
	}

	//##ModelId=4E4A89B202FD
	RSRPCModel( const RSRPCModel& rhs);
	//##ModelId=4E4A89B2031C
	RSRPCModel& operator=( const RSRPCModel& rhs);

// method
public:
	// east-lon, north-lat
	//##ModelId=4E4A89B2034B
	virtual bool FromImageToGround(double x, double y, double h, double& east, double& north)const;

	virtual bool FromGroundToImage(double east, double north, double h, double& x, double& y)const;
//	virtual bool FromImageToGround(double x, double y, double h, double& east, double& north)const;
	//##ModelId=4E4A89B2037A
	virtual DPolygon GetModelRange()const;
private:
	//##ModelId=4E4A89B20399
	virtual AffineTranslationModel GetGolbalAffineModel()const;

public:
	//##ModelId=4E4A89B203B9
	// the flag define the model direction, from ground to image should be true
	// from image to ground should be false
	bool Init( const char* strRPCFilePath, bool flag = true);

	//##ModelId=4E4A89B30000
	bool WriteRPCFile( const char* strRPCFileName);

	// set px py, invpx and invpy
	// this function is only for the true model
	// i don't think it should be the from ground to the image model
	bool  SetAffineParameter(double AffineTrans[12]);

	// create the RPC by other geometry model
	// valid range of the rpc model should be the translate from the  original model
	// and more, I need to know the translation between the original image and the transformed image
	// and the Image range, not just for the Image
	// and get from the parameter.
	// flag could be true
	// ImageSize just record the image height and width, but the range
	// the image range should be record in the ImageMappingModel
	bool CreateRPCModel(const DRect& ImageSize,	int dx, int dy, int nz,
						double MaxElevation, double MinElevation, const RSGeometryModel* OrignalModel, bool flag,
						const ImageMappingModel* ImageTransModel = &ImageUnchangeSinglelon, const ObjectMappingModel* GroundTransModel = &ObjectUnchangeSinglelon);
	bool CreateRPCModel(VirtualMachine* OriginModel,bool flag,int* ImageSize,int index,int dx,int dy,int nz);

	// the path to store the error of the model 
	//##ModelId=4E4A89B3004E
	void RPCError( const char* pPath );

private:

	//##ModelId=4E4A89B30050
	bool ReadRPCFile( const char* strRPCFilePath, const char* strImgRangeFile = NULL);

private:
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////It's not necessary, actually////////////////////////////////////////////////
	// just for the forward RFM and backward RFM based forward 
	double distanceOnePixel(double latitude,double longitude)const;

	void PreciseBasedAffine(double &latitude,double &longitude,double dlat,double dlon,double x0,double y0,double h0)const;

	void LATLONGHEIGHT2LineSample(double latitude,double longitude,double height,double &x,double &y)const;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void LineSample2LATLONGHEIGHT(double line,double sample,double height,double &latitude,double &longitude)const;

	// computer the RPC parameter
	void RPCPLH2a(double P,double L,double H,double *a);


	void DoRPC1();//分母不同，阶数为1	
	void DoRPC2();//分母不同，阶数为2	
	void DoRPC3();//分母不同，阶数为3	
	void DoRPC4();//分母相同不为1，阶数为1	
	void DoRPC5();//分母不同不为1，阶数为2	
	void DoRPC6();//分母不同不为1，阶数为3	
	void DoRPC7();//无分母，阶数为1	
	void DoRPC8();//无分母，阶数为2	
	void DoRPC9();//无分母，阶数为3
	void IDoRPC3();//迭代求解
	
};

 
#endif // !defined(AFX_RPCMODEL_H__91231331_BBBB_4527_B780_D8F0BEA9DDBB__INCLUDED_)
