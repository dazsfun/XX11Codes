// RSGeometryModel.h: interface for the RSGeometryModel class.
//
//////////////////////////////////////////////////////////////////////

#if !defined RSGEOMETRYMODEL_H_
#define RSGEOMETRYMODEL_H_


#include "Curve.h"
class GDALRasterBand;


/******************************************************************************/
/* this is affine translation model for local region or inaccuracy model      */
/* it could be set up by other models                                         */
/* the rpc model is a more complex model                                      */
/* 还需要实现的功能有 必须初始化才能调用相关函数                              */
/*                                                                            */
/*                                                         _                  */
/*  Author: 潘红播                                   |_|  |_)                 */
/*  Email:  hongbo640@gmail.com                      | |  |_)                 */
/*  Edit time: 2011 2 26                                                      */
/*  Last modified: 2011 2 27                                                  */
/******************************************************************************/

class RSGeometryModel;           // declaration 

//##ModelId=4E4A89B20128
class __declspec(dllexport) AffineTranslationModel
{
public:
	CSatOrbit helper;
	//##ModelId=4E4A89B20129
	AffineTranslationModel(){}                  // inline
	//##ModelId=4E4A89B2012A
	virtual ~AffineTranslationModel(){}
	
	//##ModelId=4E4A89B20138
	bool fromImageToGround(double x, double y, double& lat,double& lon) const;
	//##ModelId=4E4A89B2013E
	bool fromGroundToImage(double lat, double lon, double& x, double& y) const;
	
	//##ModelId=4E4A89B20144
	bool createAffineTranslationModel(const RSGeometryModel* pModel, double x, double y, double h, double dx, double dy);
	
	bool createAffineTranslationModel(int nSize, double* lat, double* lon, double* x, double* y);
private:
	double m_LatOff, m_LatScale;
	double m_LonOff, m_LonScale;
	double m_SampOff, m_SampScale;
	double m_LineOff, m_LineScale;
	//##ModelId=4E4A89B2014B
	double m_H;                     // 仿射变换成立的高程平面
	//##ModelId=4E4A89B2014C
	double m_ps[3];
	//##ModelId=4E4A89B2014D
	double m_pl[3];
};



/************************************************************************/
/* Abstract Base Class                                                  */
/* Here Supply three operation:                                         */
/* Get The Model validing range                                         */
/* From image coordinate to the ground coordinate                       */
/* From ground coordinate to the image coordinate                       */
/*                                                                      */
/*                                                         _            */
/*  Author: 潘红播                                   |_|  |_)           */
/*  Email:  hongbo640@gmail.com                      | |  |_)           */
/*  Edit time: 2010 3 04                                                */
/*  Last modified: 2011 2 27                                            */
/************************************************************************/

//##ModelId=4E4A89B20157
class __declspec(dllexport) RSGeometryModel  
{
public:

	CSatOrbit helper;
	//##ModelId=4E4A89B20158
	RSGeometryModel(){}
	//##ModelId=4E4A89B20159
	virtual ~RSGeometryModel(){}

	// the backward transformation is based on iteration of the forward transformation 
	// the starting value comes from the Global affine translation model
	// rpc model has been tested, it's ok 2011 03 01
	//##ModelId=4E4A89B2015B
	// 为什么这里是const, 为了利用const RSGeometryModel*来进行某些计算
	virtual bool FromGroundToImage(double east, double north, double h, double& x, double& y)const;

	// different model means different forward transformation
	//##ModelId=4E4A89B20167
	virtual bool FromImageToGround(double x, double y, double h, double& east, double& north)const=0;

	// here we supply the vertex in the image coordinate.
	// return the number of the vertex,if the polygon is invalid , return -1.
	//##ModelId=4E4A89B2016E
	virtual DPolygon GetModelRange()const=0;

private:
	//##ModelId=4E4A89B20170
	virtual AffineTranslationModel GetGolbalAffineModel()const=0;
};


bool __declspec(dllexport) ReadHeightCoarse(GDALRasterBand* poBand,double adfGeoTransform[6],double lat,double lon,double &H);
bool __declspec(dllexport) GetRegionDEMInfo(RSGeometryModel* ImgModel, const char* strDEMFile, 
											double posx, double posy, double offx, double offy, 
											double& MaxElevation, double& MinElevation, double& AveElevation);


// update again for const


// Abstract Base Class
// To build the relationship between the original images and  the transformed image
//##ModelId=4E4A89B20177
class __declspec(dllexport) ImageMappingModel  
{
public:
	//##ModelId=4E4A89B20178
	virtual ~ImageMappingModel(){}
	
	//##ModelId=4E4A89B20186
	virtual void OriginalToTransformed(double ox, double oy, double& tx, double& ty)const = 0;
	//##ModelId=4E4A89B2018C
	virtual void TransformedToOriginal(double tx, double ty, double& ox, double& oy)const = 0;
};

//##ModelId=4E4A89B20197
class __declspec(dllexport) ObjectMappingModel
{
public:
	//##ModelId=4E4A89B201A5
	virtual ~ObjectMappingModel(){}
	
	//##ModelId=4E4A89B201A7
	virtual void OriginalToTransformed(double oeast, double onorth, double oh, double& teast, double& tnorth, double& th)const = 0;
	//##ModelId=4E4A89B201B6
	virtual void TransformedToOriginal(double teast, double tnorth, double th, double& oeast, double& onorth, double& oh)const = 0;
};

//##ModelId=4E4A89B201C5
class __declspec(dllexport) ImageUnchangeModel : public ImageMappingModel
{
public:
	//##ModelId=4E4A89B201D5
	ImageUnchangeModel(){}
	//##ModelId=4E4A89B201D6
	~ImageUnchangeModel(){}
	//##ModelId=4E4A89B201D7
	void OriginalToTransformed(double ox, double oy, double& tx, double& ty)const
	{
		tx = ox;
		ty = oy;
	}
	//##ModelId=4E4A89B201DC
	void TransformedToOriginal(double tx, double ty, double& ox, double& oy)const
	{
		ox = tx;
		oy = ty;
	}
};

//##ModelId=4E4A89B201E5
class __declspec(dllexport) ObjectUnchangeModel : public ObjectMappingModel
{
public:
	//##ModelId=4E4A89B201E7
	ObjectUnchangeModel(){}
	//##ModelId=4E4A89B201E8
	~ObjectUnchangeModel(){}
	
	//##ModelId=4E4A89B201E9
	void OriginalToTransformed(double oeast, double onorth, double oh, double& teast, double& tnorth, double& th) const
	{
		teast  = oeast;
		tnorth = onorth;
		th     = oh;
	}
	//##ModelId=4E4A89B201F9
	void TransformedToOriginal(double teast, double tnorth, double th, double& oeast, double& onorth, double& oh) const
	{
		oeast  = teast;
		onorth = tnorth;
		oh     = th;
	}
};


//##ModelId=4E4A89B20203
class __declspec(dllexport) ImageAffineTransModel : public ImageMappingModel
{
public:
	//##ModelId=4E4A89B20205
	ImageAffineTransModel();
	//##ModelId=4E4A89B20213
	~ImageAffineTransModel();
	//##ModelId=4E4A89B20214
	ImageAffineTransModel(double coe[12]);
	//##ModelId=4E4A89B20216
	ImageAffineTransModel(const ImageAffineTransModel& rhs);
	//##ModelId=4E4A89B20218
	ImageAffineTransModel& operator=(const ImageAffineTransModel& rhs);
	//##ModelId=4E4A89B2021A
	void OriginalToTransformed(double ox, double oy, double& tx, double& ty) const
	{
		tx = px[0]+px[1]*ox+px[2]*oy;
		ty = py[0]+py[1]*ox+py[2]*oy;
	}
	//##ModelId=4E4A89B20222
	void TransformedToOriginal(double tx, double ty, double& ox, double& oy) const
	{
		ox = invpx[0]+invpx[1]*tx+invpx[2]*ty;
		oy = invpy[0]+invpy[1]*tx+invpy[2]*ty;
	}
private:
	//##ModelId=4E4A89B20227
	double px[3];
	//##ModelId=4E4A89B20228
	double py[3];
	//##ModelId=4E4A89B20229
	double invpx[3];
	//##ModelId=4E4A89B2022A
	double invpy[3];
};

#endif 
