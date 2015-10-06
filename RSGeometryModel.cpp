// RSGeometryModel.cpp: implementation of the RSGeometryModel class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
//#include "RSGeometry.h"
#include "RSGeometryModel.h"
#include <gdal_priv.h>
#include "SatOrbitNew.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// ##ModelId=4E4A89B2015B
// It looks seem that, the dx and dy should be the same 
bool RSGeometryModel::FromGroundToImage(double east, double north, double h, double& x, double& y) const
{
	double& lat = north;
	double& lon = east;
	double latp,lat2,lonp,lon2;
	GetGolbalAffineModel().fromGroundToImage(lat, lon, x, y); // 
	
	FromImageToGround(x,y,h,lonp,latp);                 // 1 1
	FromImageToGround(x+10,y+10,h,lon2,lat2);           // 1 2
	double pixelsize2 = ((latp-lat2)*(latp-lat2)+(lonp-lon2)*(lonp-lon2))/(200.0);
	double pixelsize = sqrt(pixelsize2 );
	
	// 
	AffineTranslationModel LocalAffineModel;
	double dy = fabs((latp-lat)/pixelsize);
	double dx = fabs((lonp-lon)/pixelsize);
	if(dx < 0.001 && dy < 0.001)
		return true;
	double dd;   
	// because of lower accuracy, so this should be ok
//	double dd = dx > dy ? dx : dy;
//	if( dd < 10) dd=10;
	LocalAffineModel.createAffineTranslationModel(this, x, y, h, dx, dy);  
	
	double e,e0;
	e=((latp-lat)*(latp-lat)+(lonp-lon)*(lonp-lon))/pixelsize2;
    int num=0;
	do {
		e0=e;

		// compute the difference between affine model and FromImageToGround 
		LocalAffineModel.fromGroundToImage(lat, lon, x, y);
		FromImageToGround(x,y,h,lonp,latp);
		e=((latp-lat)*(latp-lat)+(lonp-lon)*(lonp-lon))/pixelsize2;

		// this is accuracy, should be limited in 0.001 pixel
		if(e<0.000001) break;
		// if the model couldn't be converge, the local affine model should be recalculated
		// here can't be 10 pixel, should be 10, may be very small

		// 重新设定范围
		// because the pixel size shouldn't be recompute like this
//		pixelsize = sqrt(((latp-lat2)*(latp-lat2)+(lonp-lon2)*(lonp-lon2))/(dd*dd) );
		dy = fabs((latp-lat)/pixelsize);
		dx = fabs((lonp-lon)/pixelsize);
		dd = dx > dy ? dx : dy;
		// this shouldn't be limited, the accuracy will be lower
		// the gird of 1 pixel should couldn't achieve enough accuracy
//		if( dd < 1) dd=1;  
		LocalAffineModel.createAffineTranslationModel(this, x, y, h, dd, dd);
				
		// the maximum calculation times should be limited to 10
		num ++;
	} while(num<10);
	
	return true;
}


void ConvertGeoReference2pixel(double adfGeoTransform[6],double X,double Y,double &I,double &J)
{
	I=(X-adfGeoTransform[0])/adfGeoTransform[1];
	J=(Y-adfGeoTransform[3])/adfGeoTransform[5];
}

// read the height from the image
bool __declspec(dllexport) ReadHeightCoarse(GDALRasterBand* poBand,double adfGeoTransform[6],double lat,double lon,double &H)
{
	double I,J;
	double dx,dy;
	lon=lon;///PI*180;
	lat=lat;///PI*180;
	ConvertGeoReference2pixel(adfGeoTransform,lon,lat,I,J);
	int i=INT(I);
	int j=INT(J);
	dx=I-i;dy=J-j;
	double m_NoData=0;
	short bandgray[4];

	if ( i+1<0 || i>poBand->GetXSize()-1 ||
		 j+1<0 || j>poBand->GetYSize()-1 )
	{
		H = -9999;
		return false;
	}

	GDALDataType  eType = poBand->GetRasterDataType();
	poBand->RasterIO(GF_Read, i, j, 2, 2, bandgray, 2, 2, eType, 0, 0);
	for(int nn=0;nn<4;nn++)
	{
		if ( fabs(double(bandgray[nn]+9999))<0.1 )
		{
			bandgray[nn]=(short)m_NoData;
		}
	}	
				
	H=(1-dx)*(1-dy)*bandgray[0]+dx*(1-dy)*bandgray[1]+
		(1-dx)*dy*bandgray[2]+dx*dy*bandgray[3];
	return true;
}


// 修改增加读取平均高程 add by pan 2011 02 22
bool __declspec(dllexport) GetRegionDEMInfo(RSGeometryModel* ImgModel, const char* strDEMFile,
											double posx, double posy, double offx, double offy,
											double& MaxElevation, double& MinElevation, double& AveElevation)
{
	GDALDataset* poDataset;
	GDALAllRegister();
	poDataset = (GDALDataset *)GDALOpen( strDEMFile, GA_ReadOnly );
	
	if ( NULL == poDataset)
	{
		return false;
	}

	int nBandNum = poDataset->GetRasterCount();
	GDALRasterBand* poBand = poDataset->GetRasterBand(1);

	// geo information 
	double adfGeoTransform[6];
	GDALGetGeoTransform( poDataset, adfGeoTransform);
	
	// 
	int xstep = 200;
	int ystep = 200;
	int m_nx=(int)(offx/xstep)+1;
	int m_ny=(int)(offy/ystep)+1;

	// image range
// 	DPoint2D vPo[100];
// 	int nn = .GetPolygon( vPo);
	DPolygon ImgRange = ImgModel->GetModelRange();

    //获得高程范围
	double lat,lon,DemH,x,y;
	MaxElevation=-10000,MinElevation=10000;
	AveElevation=0;
	double nNum = m_nx*m_ny;
	int i,j;
	for(i=0;i<m_ny;i++)//Height
	{
		for(j=0;j<m_nx;j++)//Width
		{
			x = posx + xstep*j;
			y = posy + ystep*i;
			if ( !ImgRange.InPolygon( DPoint2D( x, y)) )
			{
				continue;
			}
			if ( ImgModel->FromImageToGround(x, y, 0, lon, lat))
			{
				ReadHeightCoarse(poBand, adfGeoTransform, lat, lon, DemH);
				AveElevation+=(DemH/nNum);
				if(DemH>MaxElevation)MaxElevation=DemH;
				if(DemH<MinElevation)MinElevation=DemH;
			}
		}
	}

	return true;
}

/************************************************************************/
/* the affine translation model for image space                         */
/* haven't be normalized                                                */
/************************************************************************/


//##ModelId=4E4A89B20144
bool AffineTranslationModel::createAffineTranslationModel(const RSGeometryModel* pModel, double x, double y, double h, double dx, double dy)
{
	double s[5],l[5],lat[5],lon[5];
	s[0]=x;                      l[0]=y;
	s[1]=x-dx;                   l[1]=y-dy;
	s[2]=x+dx;                   l[2]=y-dy;
	s[3]=x+dx;                   l[3]=y+dy;
	s[4]=x-dx;                   l[4]=y+dy;

	m_H = h;
	
	for(int i=0;i<5;i++)
		pModel->FromImageToGround(s[i],l[i],m_H,lon[i],lat[i]);

	return createAffineTranslationModel(5, lat, lon, s, l);
/*	
    // replace by new method
	Compute_avAnddx(lat,5,m_LatOff,m_LatScale);
	Compute_avAnddx(lon,5,m_LonOff,m_LonScale);	
	Compute_avAnddx(l,5,m_LineOff,m_LineScale);
	Compute_avAnddx(s,5,m_SampOff,m_SampScale);			
	
    //归一化
	for(i=0;i<5;i++)
	{
		l[i]=(l[i]-m_LineOff)/m_LineScale;
		s[i]=(s[i]-m_SampOff)/m_SampScale;
		
		lat[i]=(lat[i]-m_LatOff)/m_LatScale;
		lon[i]=(lon[i]-m_LonOff)/m_LonScale;	
	}
    double aas[9],aal[9],a[3];
	
	memset(aas,0,sizeof(double)*9);
	memset(aal,0,sizeof(double)*9);
	memset(m_ps,0,sizeof(double)*3);
	memset(m_pl,0,sizeof(double)*3);
	for(i=0;i<5;i++)
	{
		a[0]=1,a[1]=lat[i],a[2]=lon[i];
		pNormal(a,3,s[i],aas,m_ps,1);
		pNormal(a,3,l[i],aal,m_pl,1);
	}

	solve33(m_ps, aas);
	solve33(m_pl, aal);

	return true;
*/
}

bool AffineTranslationModel::createAffineTranslationModel(int nSize, double* lat, double* lon, double* x, double* y)
{
	int i = 0;
	if (nSize < 3) return false;

	helper.Compute_avAnddx(lat, nSize, m_LatOff, m_LatScale);
	helper.Compute_avAnddx(lon, nSize, m_LonOff, m_LonScale);
	helper.Compute_avAnddx(y, nSize, m_LineOff, m_LineScale);
	helper.Compute_avAnddx(x, nSize, m_SampOff, m_SampScale);
	
    //归一化
	for(int i=0;i<5;i++)
	{
		y[i]=(y[i]-m_LineOff)/m_LineScale;
		x[i]=(x[i]-m_SampOff)/m_SampScale;
		
		lat[i]=(lat[i]-m_LatOff)/m_LatScale;
		lon[i]=(lon[i]-m_LonOff)/m_LonScale;	
	}
    double aas[9],aal[9],a[3];
	
	memset(aas,0,sizeof(double)*9);
	memset(aal,0,sizeof(double)*9);
	memset(m_ps,0,sizeof(double)*3);
	memset(m_pl,0,sizeof(double)*3);
	for(i=0;i<5;i++)
	{
		a[0]=1,a[1]=lat[i],a[2]=lon[i];
		helper.pNormal(a, 3, x[i], aas, m_ps, 1);
		helper.pNormal(a, 3, y[i], aal, m_pl, 1);
	}
	
	helper.solve33(m_ps, aas);
	helper.solve33(m_pl, aal);
	
	return true;
}

// haven't be debuged 
//##ModelId=4E4A89B20138
bool AffineTranslationModel::fromImageToGround(double x, double y, double& lat,double& lon) const
{
	//
	double xx = (x-m_SampOff)/m_SampScale - m_ps[0];
	double yy = (y-m_LineOff)/m_LineScale - m_pl[0];

	double deno = m_ps[1]*m_pl[2] - m_ps[2]*m_pl[1];
	lat = (m_ps[1]*yy - m_pl[1]*xx)/deno;
	lon = (m_pl[2]*xx-m_ps[2]*yy)/deno;

	lat = lat*m_LatScale+m_LatOff;
	lon = lon*m_LonScale+m_LonOff;
	return true;
}

//##ModelId=4E4A89B2013E
bool AffineTranslationModel::fromGroundToImage(double lat, double lon, double& x, double& y) const
{
	lat = (lat-m_LatOff)/m_LatScale;
	lon = (lon-m_LonOff)/m_LonScale;
	x=m_ps[0]+m_ps[1]*lat+m_ps[2]*lon;
	y=m_pl[0]+m_pl[1]*lat+m_pl[2]*lon;
	
	y=y*m_LineScale+m_LineOff;
	x=x*m_SampScale+m_SampOff;

	return true;
}



/************************************************************************/
/* the affine translation model for image space                         */
/* haven't be normalized                                                */
/************************************************************************/

//##ModelId=4E4A89B20205
ImageAffineTransModel::ImageAffineTransModel()
{
	memset( px, 0, sizeof(double)*3);
	memset( py, 0, sizeof(double)*3);
	memset( invpx, 0, sizeof(double)*3);
	memset( invpy, 0, sizeof(double)*3);
	px[1] = py[2] = 1;
	invpx[1] = invpy[2] = 1;
}

//##ModelId=4E4A89B20213
ImageAffineTransModel::~ImageAffineTransModel()
{
	
}

//##ModelId=4E4A89B20214
ImageAffineTransModel::ImageAffineTransModel(double coe[12])
{
	memcpy( px, coe, sizeof(double)*3);
	memcpy( py, coe+3, sizeof(double)*3);
	memcpy( invpx, coe+6, sizeof(double)*3);
	memcpy( invpy, coe+9, sizeof(double)*3);
}

//##ModelId=4E4A89B20216
ImageAffineTransModel::ImageAffineTransModel(const ImageAffineTransModel& rhs)
{
	memcpy( px, rhs.px, sizeof(double)*3);
	memcpy( py, rhs.py, sizeof(double)*3);
	memcpy( invpx, rhs.invpx, sizeof(double)*3);
	memcpy( invpy, rhs.invpy, sizeof(double)*3);
}

//##ModelId=4E4A89B20218
ImageAffineTransModel& ImageAffineTransModel::operator=(const ImageAffineTransModel& rhs)
{
	if( this == &rhs) return *this;
	else 
	{
		memcpy( px, rhs.px, sizeof(double)*3);
		memcpy( py, rhs.py, sizeof(double)*3);
		memcpy( invpx, rhs.invpx, sizeof(double)*3);
		memcpy( invpy, rhs.invpy, sizeof(double)*3);
	}
	return *this;
}

