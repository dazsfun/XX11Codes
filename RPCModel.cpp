// RSRPCModel.cpp: implementation of the RSRPCModel class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
//#include "RSGeometry.h"
#include "RPCModel.h"
#include "SatOrbitNew.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//##ModelId=4E4A89B202BF
RSRPCModel::RSRPCModel()
{

	// read form rpc file
	// default construction
	// m_ModelRange
	memset(LNUM, 0, sizeof(double)*20);
	memset(LDEN, 0, sizeof(double)*20);
	memset(SNUM, 0, sizeof(double)*20);
	memset(SDEN, 0, sizeof(double)*20);
	memset(px,0,sizeof(double)*3);
	memset(py,0,sizeof(double)*3);
	memset(invpx,0,sizeof(double)*3);
	memset(invpy,0,sizeof(double)*3);
	SDEN[19] = LDEN[19] = 1.0;
	px[1] = py[2] = invpx[1] = invpy[2] = 1;
	HEIGHT_SCALE = LINE_SCALE = SAMP_SCALE = LAT_SCALE = LONG_SCALE = 1;
	HEIGHT_OFF = LINE_OFF = SAMP_OFF = LAT_OFF = LONG_OFF = 0;

	// compute the rpc model
	// default construction
	// m_ImgRange
	m_pOrignalModel = NULL;
	m_pImageTrans   = NULL;
	m_pGroundTrans  = NULL;
	m_MinElevation  = 0;
	m_MaxElevation  = 0;
	m_nx            = 0;
	m_ny            = 0;
	m_nz            = 0;

	m_nRPCgcpNum = 0;
	m_RPCgcp = NULL;
    m_Y = NULL;
	m_X = NULL;
	m_P = NULL;
	m_L = NULL;
	m_H = NULL;

	m_bflag = true;
}

//##ModelId=4E4A89B202DE
RSRPCModel::~RSRPCModel()
{
	if(m_Y)delete []m_Y;m_Y=NULL;
	if(m_X)delete []m_X;m_X=NULL;
	if(m_P)delete []m_P;m_P=NULL;
	if(m_L)delete []m_L;m_L=NULL;
	if(m_H)delete []m_H;m_H=NULL;
	
	for(int i=0;i<m_nz;i++)
	{
		if(m_RPCgcp[i].ddx)delete []m_RPCgcp[i].ddx;m_RPCgcp[i].ddx=NULL;
		if(m_RPCgcp[i].ddy)delete []m_RPCgcp[i].ddy;m_RPCgcp[i].ddy=NULL;
	}
	if(m_RPCgcp)delete []m_RPCgcp;m_RPCgcp=NULL;
}


//##ModelId=4E4A89B202FD
RSRPCModel::RSRPCModel(const RSRPCModel& rhs):m_ModelRange( rhs.m_ModelRange),m_ImgRange(rhs.m_ImgRange),
					 HEIGHT_SCALE(rhs.HEIGHT_SCALE), HEIGHT_OFF(rhs.HEIGHT_OFF),
					 LINE_SCALE(rhs.LINE_SCALE), LINE_OFF(rhs.LINE_OFF),
					 SAMP_SCALE(rhs.SAMP_SCALE), SAMP_OFF(rhs.SAMP_OFF),
					 LAT_SCALE(rhs.LAT_SCALE), LAT_OFF(rhs.LAT_OFF),
					 LONG_SCALE(rhs.LONG_SCALE),LONG_OFF(rhs.LONG_OFF)
{
	memcpy(LNUM, rhs.LNUM, sizeof(double)*20);
	memcpy(LDEN, rhs.LDEN, sizeof(double)*20);
	memcpy(SNUM, rhs.SNUM, sizeof(double)*20);
	memcpy(SDEN, rhs.SDEN, sizeof(double)*20);
	memcpy(px,rhs.px,sizeof(double)*3);
	memcpy(py,rhs.py,sizeof(double)*3);
	memcpy(invpx,rhs.invpx,sizeof(double)*3);
	memcpy(invpy,rhs.invpy,sizeof(double)*3);

	// this and image range is not need
	m_MinElevation  = rhs.m_MinElevation;
	m_MaxElevation  = rhs.m_MaxElevation;
	// these should be zero
	// 
	m_pOrignalModel = NULL;
	m_pImageTrans   = NULL;
	m_pGroundTrans  = NULL;
	m_nx            = 0;
	m_ny            = 0;
	m_nz            = 0;

	m_nRPCgcpNum = 0;
	m_RPCgcp = NULL;
    m_Y = NULL;
	m_X = NULL;
	m_P = NULL;
	m_L = NULL;
	m_H = NULL;
}

//##ModelId=4E4A89B2031C
RSRPCModel& RSRPCModel::operator=( const RSRPCModel& rhs)
{
	if( this == &rhs) return *this;
	else 
	{
		m_ModelRange =  rhs.m_ModelRange;
		HEIGHT_SCALE =  rhs.HEIGHT_SCALE;
		HEIGHT_OFF   =  rhs.HEIGHT_OFF;
		LINE_SCALE   =  rhs.LINE_SCALE;
		LINE_OFF     =  rhs.LINE_OFF;
		SAMP_SCALE   =  rhs.SAMP_SCALE;
		SAMP_OFF     =  rhs.SAMP_OFF;
		LAT_SCALE    =  rhs.LAT_SCALE;
		LAT_OFF      =  rhs.LAT_OFF;
		LONG_SCALE   =  rhs.LONG_SCALE;
		LONG_OFF     =  rhs.LONG_OFF;
		
		memcpy(LNUM, rhs.LNUM, sizeof(double)*20);
		memcpy(LDEN, rhs.LDEN, sizeof(double)*20);
		memcpy(SNUM, rhs.SNUM, sizeof(double)*20);
		memcpy(SDEN, rhs.SDEN, sizeof(double)*20);
		memcpy(px,rhs.px,sizeof(double)*3);
		memcpy(py,rhs.py,sizeof(double)*3);
		memcpy(invpx,rhs.invpx,sizeof(double)*3);
		memcpy(invpy,rhs.invpy,sizeof(double)*3);
		
		// these should be zero
		// 
		m_ImgRange   =  rhs.m_ImgRange;
		m_MinElevation  = rhs.m_MinElevation;
		m_MaxElevation  = rhs.m_MaxElevation;
	
		m_bflag = rhs.m_bflag;

		m_nx            = 0;
		m_ny            = 0;
		m_nz            = 0;
		m_nRPCgcpNum    = 0;
		
		if(m_Y)delete []m_Y;m_Y=NULL;
		if(m_X)delete []m_X;m_X=NULL;
		if(m_P)delete []m_P;m_P=NULL;
		if(m_L)delete []m_L;m_L=NULL;
		if(m_H)delete []m_H;m_H=NULL;
		
		for(int i=0;i<m_nz;i++)
		{
			if(m_RPCgcp[i].ddx)delete []m_RPCgcp[i].ddx;m_RPCgcp[i].ddx=NULL;
			if(m_RPCgcp[i].ddy)delete []m_RPCgcp[i].ddy;m_RPCgcp[i].ddy=NULL;
		}
		if(m_RPCgcp)delete []m_RPCgcp;m_RPCgcp=NULL;
		
		m_pOrignalModel = NULL;
		m_pImageTrans   = NULL;
		m_pGroundTrans  = NULL;
	}
	return *this;
}

// 利用RPC文件读取RPC参数
//##ModelId=4E4A89B30050
bool RSRPCModel::ReadRPCFile( const char* strRPCFilePath, const char* strImgRangeFile)
{
	FILE *fp;
	fopen_s(&fp,strRPCFilePath,"r");
	if (NULL == fp) {
		return false;
	}
	
	px[0]=0; px[1]=1; px[2]=0;
	py[0]=0; py[1]=0; py[2]=1;
				
	int i;
	char line[256];

	if ( m_bflag == true)
	{
		fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256); 
		LINE_OFF=atof(line);
		fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); SAMP_OFF=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LAT_OFF=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LONG_OFF=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); HEIGHT_OFF=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LINE_SCALE=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); SAMP_SCALE=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LAT_SCALE=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LONG_SCALE=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); HEIGHT_SCALE=atof(line);fscanf_s(fp,"%s",line,256);
	}
	else
	{
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LAT_OFF=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LONG_OFF=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LINE_OFF=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); SAMP_OFF=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); HEIGHT_OFF=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LAT_SCALE=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LONG_SCALE=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); LINE_SCALE=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); SAMP_SCALE=atof(line);fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);fscanf_s(fp,"%s",line,256); HEIGHT_SCALE=atof(line);fscanf_s(fp,"%s",line,256);
	}


	for(i=0;i<20;i++)
	{
		fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);
		 LNUM[i]=atof(line);
	}
	
	for(i=0;i<20;i++)
	{
		fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);
		 LDEN[i]=atof(line);
	}
	
	for(i=0;i<20;i++)
	{
		fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);
		 SNUM[i]=atof(line);
	}
	
	for(i=0;i<20;i++)
	{
		fscanf_s(fp,"%s",line,256);
		fscanf_s(fp,"%s",line,256);
		 SDEN[i]=atof(line);
	}
	fclose(fp);
	
	if ( NULL == strImgRangeFile )
	{
		DPoint2D ver[4];
		ver[0].x = ver[3].x = SAMP_OFF - SAMP_SCALE;
		ver[1].x = ver[2].x = SAMP_OFF + SAMP_SCALE;
		ver[0].y = ver[1].y = LINE_OFF - LINE_SCALE;
		ver[2].y = ver[3].y = LINE_OFF + LINE_SCALE;
		m_ModelRange = DPolygon(ver, 4);
	} 
	else
	{
		FILE* fp;
		fopen_s(&fp,strImgRangeFile, "r");
		unsigned num;
		fscanf_s(fp, "%d" , &num );
		
		DPoint2D* po = new DPoint2D[num];
		int x,y;
		for (unsigned i=0; i<num; ++i)
		{
			fscanf_s(fp, "%d %d", &x, &y);
			po[i].x = x, po[i].y = y;
		}
		fclose (fp);
		DPolygon DPo(po, static_cast<unsigned> (num) );
		delete []po;
		m_ModelRange = DPo;
	}

	return true;
}

//##ModelId=4E4A89B203B9
bool RSRPCModel::Init( const char* strRPCFilePath, bool flag)
{
	// find range file
	m_bflag = flag;
	if( !ReadRPCFile(strRPCFilePath)) return false;

	if ( !m_GlobalAffineModel.createAffineTranslationModel(this, SAMP_OFF, LINE_OFF, 0, SAMP_SCALE, LINE_SCALE)) return false;
	
	return true;
}


// 将RPC参数写入RPC文件
//##ModelId=4E4A89B30000
bool RSRPCModel::WriteRPCFile( const char* strRPCFileName)
{
	FILE *fp;
	fopen_s(&fp, strRPCFileName,"w");
	if(fp==NULL) return false;
	if ( m_bflag == true)
	{
		fprintf(fp,"LINE_OFF: %+010.2lf pixels\n", LINE_OFF);
		fprintf(fp,"SAMP_OFF: %+010.2lf pixels\n", SAMP_OFF);
		fprintf(fp,"LAT_OFF: %+011.8lf degrees\n", LAT_OFF);
		fprintf(fp,"LONG_OFF: %+012.8lf degrees\n", LONG_OFF);
		fprintf(fp,"HEIGHT_OFF: %+08.3lf meters\n", HEIGHT_OFF);
		
		fprintf(fp,"LINE_SCALE: %+010.2lf pixels\n", LINE_SCALE);
		fprintf(fp,"SAMP_SCALE: %+010.2lf pixels\n", SAMP_SCALE);
		fprintf(fp,"LAT_SCALE: %0+12.8lf degrees\n", LAT_SCALE);
		fprintf(fp,"LONG_SCALE: %+012.8lf degrees\n", LONG_SCALE);
		fprintf(fp,"HEIGHT_SCALE: %+08.3lf meters\n", HEIGHT_SCALE);
	}
	else
	{
		fprintf(fp,"LINE_OFF: %+010.2lf pixels\n", LAT_OFF);
		fprintf(fp,"SAMP_OFF: %+010.2lf pixels\n", LONG_OFF);
		fprintf(fp,"LAT_OFF: %+011.8lf degrees\n", LINE_OFF);
		fprintf(fp,"LONG_OFF: %+012.8lf degrees\n", SAMP_OFF);
		fprintf(fp,"HEIGHT_OFF: %+08.3lf meters\n", HEIGHT_OFF);
		
		fprintf(fp,"LINE_SCALE: %+010.2lf pixels\n", LAT_OFF);
		fprintf(fp,"SAMP_SCALE: %+010.2lf pixels\n", LONG_SCALE);
		fprintf(fp,"LAT_SCALE: %0+12.8lf degrees\n", LINE_SCALE);
		fprintf(fp,"LONG_SCALE: %+012.8lf degrees\n", SAMP_SCALE);
		fprintf(fp,"HEIGHT_SCALE: %+08.3lf meters\n", HEIGHT_SCALE);
	}

	int i;
	for( i=0;i<20;i++)
	{
		fprintf(fp,"LINE_NUM_COEFF_%d:  %+50.40le\n",i+1, LNUM[i]);
	}
	
	for( i=0;i<20;i++)
	{
		fprintf(fp,"LINE_DEN_COEFF_%d:  %+50.40le\n",i+1, LDEN[i]);
	}
	
	for( i=0;i<20;i++)
	{
		fprintf(fp,"SAMP_NUM_COEFF_%d:  %+50.40le\n",i+1, SNUM[i]);
	}
	for( i=0;i<20;i++)
	{
		fprintf(fp,"SAMP_DEN_COEFF_%d:  %+50.40le\n",i+1, SDEN[i]);
	}
	fclose(fp);
	return true;
}

// 计算某个经纬度上像元大小
//##ModelId=4E4A89B3007D
double RSRPCModel::distanceOnePixel(double latitude,double longitude)const
{
	double x1,y1,x2,y2;
	LATLONGHEIGHT2LineSample(latitude,longitude, HEIGHT_OFF,x1,y1);
	LATLONGHEIGHT2LineSample(latitude+0.01* LAT_SCALE,longitude, HEIGHT_OFF,x2,y2);

	return  0.01*LAT_SCALE/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) ;
}

// 基于部分区域计算出对应的放射变换参数
//##ModelId=4E4A89B300AB
void RSRPCModel::PreciseBasedAffine(double &latitude,double &longitude,double dlat,double dlon,double x0,double y0,double h0)const
{
	int i =0;
	double s[5],l[5],lat[5],lon[5];
	lat[0]=latitude     ;lon[0]=longitude;
	lat[1]=latitude+dlat;lon[1]=longitude+dlon;
	lat[2]=latitude+dlat;lon[2]=longitude-dlon;
	lat[3]=latitude-dlat;lon[3]=longitude+dlon;
	lat[4]=latitude-dlat;lon[4]=longitude-dlon;

	for(int i=0;i<5;i++)
	{
		LATLONGHEIGHT2LineSample( lat[i],lon[i],h0,s[i],l[i]);
	}
    double LAT_OFF,LONG_OFF,LINE_OFF,SAMP_OFF,LAT_SCALE,LONG_SCALE,LINE_SCALE,SAMP_SCALE;
	helper.Compute_avAnddx(lat,5,LAT_OFF,LAT_SCALE);
	helper.Compute_avAnddx(lon,5,LONG_OFF,LONG_SCALE);	
	helper.Compute_avAnddx(l,5,LINE_OFF,LINE_SCALE);
	helper.Compute_avAnddx(s,5,SAMP_OFF,SAMP_SCALE);

	//归一化
	for(i=0;i<5;i++)
	{
		l[i]=(l[i]-LINE_OFF)/LINE_SCALE;
		s[i]=(s[i]-SAMP_OFF)/SAMP_SCALE;
		
		lat[i]=(lat[i]-LAT_OFF)/LAT_SCALE;
		lon[i]=(lon[i]-LONG_OFF)/LONG_SCALE;		
	}

	double ps[3],pl[3];
	//P=ps[0]+ps[1]*s+ps[2]*l;
	//L=pl[0]+pl[1]*s+pl[2]*l;
    double aas[9],aal[9],a[3];
	memset(aas,0,sizeof(double)*9);
	memset(aal,0,sizeof(double)*9);
	memset(ps,0,sizeof(double)*3);
	memset(pl,0,sizeof(double)*3);
	for(i=0;i<5;i++)
	{
		a[0]=1,a[1]=s[i],a[2]=l[i];
		helper.pNormal(a,3,lat[i],aas,ps,1);
		helper.pNormal(a,3,lon[i],aal,pl,1);
	}
	helper.Gauss(aas,ps,3);
	helper.Gauss(aal,pl,3);
	
	
	y0=(y0-LINE_OFF)/LINE_SCALE;
	x0=(x0-SAMP_OFF)/SAMP_SCALE;

	double P=ps[0]+ps[1]*x0+ps[2]*y0;
	double L=pl[0]+pl[1]*x0+pl[2]*y0;

	latitude=P*LAT_SCALE+LAT_OFF;
	longitude=L*LONG_SCALE+LONG_OFF;
}

// RPC模型正解
//##ModelId=4E4A89B300EA
void RSRPCModel::LATLONGHEIGHT2LineSample(double latitude,double longitude,double height,double &x,double &y)const
{
	double P,L,H;
	P=(latitude-LAT_OFF)/LAT_SCALE;
	L=(longitude-LONG_OFF)/LONG_SCALE;
	H=(height-HEIGHT_OFF)/HEIGHT_SCALE;

	double a[20];
	a[0]=1;a[1]=L;a[2]=P;a[3]=H;a[4]=L*P;a[5]=L*H;a[6]=P*H;a[7]=L*L;a[8]=P*P;a[9]=H*H;
	a[10]=P*L*H;a[11]=L*L*L;a[12]=L*P*P;a[13]=L*H*H;a[14]=L*L*P;a[15]=P*P*P;a[16]=P*H*H;
	a[17]=L*L*H;a[18]=P*P*H;a[19]=H*H*H;
    
	double LineNum=0,LineDen=0,SampNum=0,SampDen=0;
	
	for(int i=0;i<20;i++)
	{
		LineNum+=LNUM[i]*a[i];
		LineDen+=LDEN[i]*a[i];

		SampNum+=SNUM[i]*a[i];
		SampDen+=SDEN[i]*a[i];		
	}
	double Y,X;
	Y=LineNum/LineDen;
	X=SampNum/SampDen;

	y=Y*LINE_SCALE+LINE_OFF;
	x=X*SAMP_SCALE+SAMP_OFF;
}

// 基于RPC正解的反解法
//##ModelId=4E4A89B300FA
void RSRPCModel::LineSample2LATLONGHEIGHT(double line,double sample,double height,double &latitude,double &longitude)const
{
	latitude  = LAT_OFF;
	longitude = LONG_OFF;
	PreciseBasedAffine( latitude, longitude, LAT_SCALE*0.5, LONG_SCALE*0.5, sample,line,height);

	double xp,yp;
	LATLONGHEIGHT2LineSample(latitude,longitude,height,xp,yp);
	double e=(xp-sample)*(xp-sample)+(yp-line)*(yp-line);
	double e0;
	do {
		e0=e;
		double pixelsize=distanceOnePixel(latitude,longitude);
		
		LATLONGHEIGHT2LineSample( latitude,longitude,height,xp,yp);

		PreciseBasedAffine( latitude,longitude,fabs(yp-line)*pixelsize,fabs(xp-sample)*pixelsize,
			sample,line,height);
		LATLONGHEIGHT2LineSample( latitude,longitude,height,xp,yp);
		
		e=(xp-sample)*(xp-sample)+(yp-line)*(yp-line);

		if(e<0.000001) break;
		
	} while(e<e0);
}

// 
//##ModelId=4E4A89B2034B
bool RSRPCModel::FromImageToGround(double x, double y, double h, double& east, double& north) const
{
	if (false)
	{
		// 迭代处理
		LineSample2LATLONGHEIGHT(y, x, h, north, east);
	}
	else
	{
		// 直接计算
		LATLONGHEIGHT2LineSample(y,x,h,east,north);
	}

	return true;
}

// 因为那个多项式模型建立的是影像与地面的对应关系，所以有差别
bool RSRPCModel::FromGroundToImage(double east, double north, double h, double& x, double& y) const
{
	if ( true)
	{
		LATLONGHEIGHT2LineSample( north, east, h, x, y);
	}
	else
	{
		RSGeometryModel::FromGroundToImage(east, north, h, x, y);
//		LineSample2LATLONGHEIGHT(north,east,h,y,x);
	}

	return true;
}

// 
//##ModelId=4E4A89B2037A
DPolygon RSRPCModel::GetModelRange() const
{
	return m_ModelRange;
}

//##ModelId=4E4A89B20399
AffineTranslationModel RSRPCModel::GetGolbalAffineModel()const
{
	return m_GlobalAffineModel;
}

//////////////////////////////////////////////////////////////////////////
// 求解RPC参数
//////////////////////////////////////////////////////////////////////////

// 
// all should be const
bool RSRPCModel::CreateRPCModel(VirtualMachine* OriginModel, bool flag, int* Imageize, int index, int dx, int dy, int nz)
{
	m_bflag = flag;
	double max_Height, min_Height;
	max_Height = 200;
	min_Height = 0;
	m_MaxElevation = max_Height + 30;
	m_MinElevation = min_Height - 30;

	int min_X_Loc = Imageize[0];
	int min_Y_Loc = Imageize[1];
	int max_X_Loc = Imageize[2] + Imageize[0] - 1;
	int max_Y_Loc = Imageize[3] + Imageize[1] - 1;

	int m_Xsize = static_cast<int>(Imageize[2]) / dx ;
	int m_Ysize = static_cast<int>(Imageize[3]) / dy;
	int m_Hsize = nz;
	double m_hstep = (m_MaxElevation - m_MinElevation) / (m_Hsize - 1);


	m_RPCgcp = new RPCGCPSTRUCT[m_Hsize];
	if (m_RPCgcp == NULL)return false;
	for (int i = 0; i<nz; i++)
	{
		m_RPCgcp[i].h = m_MinElevation + i*nz;
		m_RPCgcp[i].x0 = min_X_Loc;
		m_RPCgcp[i].y0 = min_Y_Loc;
		m_RPCgcp[i].dx = dx;
		m_RPCgcp[i].dy = dy;
		m_RPCgcp[i].nx = m_Xsize;
		m_RPCgcp[i].ny = m_Ysize;
		m_RPCgcp[i].ddx = new double[m_RPCgcp[i].nx*m_RPCgcp[i].ny];//lon
		if (m_RPCgcp[i].ddx == NULL) return false;
		m_RPCgcp[i].ddy = new double[m_RPCgcp[i].nx*m_RPCgcp[i].ny];//lat 
		if (m_RPCgcp[i].ddy == NULL) return false;
	}
	m_nRPCgcpNum = m_Xsize*m_Ysize*m_Hsize;

	m_P = new double[m_nRPCgcpNum];
	if (m_P == NULL)return false;
	m_L = new double[m_nRPCgcpNum];
	if (m_L == NULL)return false;
	m_H = new double[m_nRPCgcpNum];
	if (m_H == NULL)return false;
	m_Y = new double[m_nRPCgcpNum];
	if (m_Y == NULL)return false;
	m_X = new double[m_nRPCgcpNum];
	if (m_X == NULL)return false;

	m_P = new double[m_nRPCgcpNum];
	if (m_P == NULL)return false;
	m_L = new double[m_nRPCgcpNum];
	if (m_L == NULL)return false;
	m_H = new double[m_nRPCgcpNum];
	if (m_H == NULL)return false;
	m_Y = new double[m_nRPCgcpNum];
	if (m_Y == NULL)return false;
	m_X = new double[m_nRPCgcpNum];
	if (m_X == NULL)return false;

	m_nRPCgcpNum = 0;
	double h = m_MinElevation;
	double x, y, lat, lon;
	DPoint3D gpo;
	DPoint2D orpo;
	int ii, j;

	for (ii = 0; ii<m_Hsize; ii++)
	{
		for (int i = 0; i<m_Ysize; i++)
		{
			cout << i << "//" << m_Ysize << "\r";
			y = min_Y_Loc + i*dy;
			for (j = 0; j<m_Xsize; j++)
			{
				x = min_X_Loc + j*dx;

				DPoint2D trpo(x, y);
				DPoint3D gpo;
				OriginModel->FromXY2LonLatTest(x, y, h, gpo.lon, gpo.lat, index);
				gpo.h = h;
				m_RPCgcp[ii].ddy[i*m_Xsize + j] = gpo.lat;
				m_RPCgcp[ii].ddx[i*m_Xsize + j] = gpo.lon;
				if (m_bflag == true)
				{
					m_P[m_nRPCgcpNum] = gpo.lat;
					m_L[m_nRPCgcpNum] = gpo.lon;
					m_Y[m_nRPCgcpNum] = y;
					m_X[m_nRPCgcpNum] = x;
				}
				else
				{
					// 反求rpc模型，此时m_L,m_P表示x,y
					// m_Y,m_X表示lat,lon
					// 那么之后的基准参数什么的都是反过来对应的
					m_P[m_nRPCgcpNum] = y;
					m_L[m_nRPCgcpNum] = x;
					m_Y[m_nRPCgcpNum] = gpo.lat;
					m_X[m_nRPCgcpNum] = gpo.lon;
				}

				m_H[m_nRPCgcpNum] = gpo.h;
				m_nRPCgcpNum++;
			}
		}
		h += m_hstep;
	}

	//计算归一化参数
	helper.Compute_avAnddx(m_P, m_nRPCgcpNum, LAT_OFF, LAT_SCALE);
	helper.Compute_avAnddx(m_L, m_nRPCgcpNum, LONG_OFF, LONG_SCALE);
	helper.Compute_avAnddx(m_H, m_nRPCgcpNum, HEIGHT_OFF, HEIGHT_SCALE);
	helper.Compute_avAnddx(m_Y, m_nRPCgcpNum, LINE_OFF, LINE_SCALE);
	helper.Compute_avAnddx(m_X, m_nRPCgcpNum, SAMP_OFF, SAMP_SCALE);

	// 归一化
	// 
	for (int i = 0; i<m_nRPCgcpNum; i++)
	{
		m_Y[i] = (m_Y[i] - LINE_OFF) / LINE_SCALE;
		m_X[i] = (m_X[i] - SAMP_OFF) / SAMP_SCALE;

		m_P[i] = (m_P[i] - LAT_OFF) / LAT_SCALE;
		m_L[i] = (m_L[i] - LONG_OFF) / LONG_SCALE;
		m_H[i] = (m_H[i] - HEIGHT_OFF) / HEIGHT_SCALE;
	}


	// computer the value of the rpc model
	IDoRPC3();

	return true;
}


bool RSRPCModel::CreateRPCModel(const DRect& ImageSize,	int dx, int dy, int nz,
								double MaxElevation, double MinElevation, const RSGeometryModel* OrignalModel, bool flag,
								const ImageMappingModel* ImageTransModel, const ObjectMappingModel* GroundTransModel)
{
	// 获取原始模型在新的影像上的范围
	int i;
	DPoint2D OrVer[100], TrVer[100];
	m_ModelRange = OrignalModel->GetModelRange();
	int  VerNum = m_ModelRange.GetPolygon( OrVer);
	for ( i=0; i<VerNum; ++i)
	{
		ImageTransModel->OriginalToTransformed( OrVer[i].x, OrVer[i].y, TrVer[i].x, TrVer[i].y);
	}

	m_ImgRange   = ImageSize;
	m_pImageTrans = ImageTransModel;
	m_pGroundTrans= GroundTransModel;
	m_pOrignalModel= OrignalModel;
	m_bflag = flag;

	//获取核线影像起始位置和范围
	double m_Minx         ,m_Miny;
	double m_ImgHeight    ,m_ImgWidth;

	// 计算RPC模型的像素间隔以及高程分层
	m_MaxElevation = MaxElevation + 30;
	m_MinElevation = MinElevation - 30;
	m_ImgHeight =m_ImgRange.bottom - m_ImgRange.top;
	m_ImgWidth  =m_ImgRange.right  - m_ImgRange.left;
	\
	// 格网的实际范围是超过原来的计划的
	m_nx = (int)( m_ImgWidth/dx )+1;
	if( (m_nx-1)*dx<m_ImgWidth )  m_nx++; 
	m_ny = (int)( m_ImgHeight/dy )+1;
	if( (m_ny-1)*dy<m_ImgHeight ) m_ny++;
	m_nz = nz;
	double dz = (m_MaxElevation-m_MinElevation)/(m_nz-1);
	m_Minx = m_ImgRange.left;
	m_Miny = m_ImgRange.top;


	m_RPCgcp=new RPCGCPSTRUCT[m_nz];
	if(m_RPCgcp==NULL)return false;	
	for( i=0;i<nz;i++)
	{
		m_RPCgcp[i].h=m_MinElevation+i*dz;
		m_RPCgcp[i].x0=m_Minx;
		m_RPCgcp[i].y0=m_Miny;
		m_RPCgcp[i].dx=dx;
		m_RPCgcp[i].dy=dy;
		m_RPCgcp[i].nx=m_nx;
		m_RPCgcp[i].ny=m_ny;
		m_RPCgcp[i].ddx=new double[m_RPCgcp[i].nx*m_RPCgcp[i].ny];//lon
		if(m_RPCgcp[i].ddx==NULL) return false;
		m_RPCgcp[i].ddy=new double[m_RPCgcp[i].nx*m_RPCgcp[i].ny];//lat 
		if(m_RPCgcp[i].ddy==NULL) return false;
	}
	m_nRPCgcpNum=m_nx*m_ny*m_nz;

	m_P=new double[m_nRPCgcpNum];
	if(m_P==NULL)return false;	
	m_L=new double[m_nRPCgcpNum];
	if(m_L==NULL)return false;	
	m_H=new double[m_nRPCgcpNum];
	if(m_H==NULL)return false;	
	m_Y=new double[m_nRPCgcpNum];
	if(m_Y==NULL)return false;	
	m_X=new double[m_nRPCgcpNum];
	if(m_X==NULL)return false;


	m_nRPCgcpNum=0;		
	double h=m_MinElevation;
	double x, y, lat, lon;
	DPoint3D gpo;
	DPoint2D orpo;
	int ii, j;
	
	for(ii=0;ii<m_nz;ii++)
	{
		for(i=0;i<m_ny;i++)
		{	
			y=m_Miny+i*dy;	
			for(j=0;j<m_nx;j++)
			{
				x=m_Minx+j*dx;

				DPoint2D trpo(x,y);
				if ( !m_ModelRange.InPolygon( trpo ) ) 
					continue;
				((ImageMappingModel*)m_pImageTrans)->TransformedToOriginal( trpo.x ,trpo.y, orpo.x, orpo.y);
				((RSGeometryModel*)m_pOrignalModel)->FromImageToGround(orpo.x, orpo.y, h, lon, lat);
				((ObjectMappingModel*)m_pGroundTrans)->OriginalToTransformed( lon, lat, h, gpo.lon, gpo.lat, gpo.h);
				
				m_RPCgcp[ii].ddy[i*m_nx+j]=gpo.lat;
				m_RPCgcp[ii].ddx[i*m_nx+j]=gpo.lon;			
				if (m_bflag == true)
				{
					m_P[m_nRPCgcpNum]=gpo.lat;
					m_L[m_nRPCgcpNum]=gpo.lon;
					m_Y[m_nRPCgcpNum]=y;
					m_X[m_nRPCgcpNum]=x;
				}
				else
				{
					// 反求rpc模型，此时m_L,m_P表示x,y
					// m_Y,m_X表示lat,lon
					// 那么之后的基准参数什么的都是反过来对应的
					m_P[m_nRPCgcpNum]=y;    
					m_L[m_nRPCgcpNum]=x;
					m_Y[m_nRPCgcpNum]=gpo.lat;
					m_X[m_nRPCgcpNum]=gpo.lon;
				}

				m_H[m_nRPCgcpNum]=gpo.h;
				m_nRPCgcpNum++;
			}
		}
		h+=dz; 
	}

	//计算归一化参数
	helper.Compute_avAnddx(m_P,m_nRPCgcpNum,LAT_OFF,LAT_SCALE);
	helper.Compute_avAnddx(m_L,m_nRPCgcpNum,LONG_OFF,LONG_SCALE);
	helper.Compute_avAnddx(m_H,m_nRPCgcpNum,HEIGHT_OFF,HEIGHT_SCALE);
	helper.Compute_avAnddx(m_Y,m_nRPCgcpNum,LINE_OFF,LINE_SCALE);
	helper.Compute_avAnddx(m_X,m_nRPCgcpNum,SAMP_OFF,SAMP_SCALE);

    // 归一化
	// 
	for(i=0;i<m_nRPCgcpNum;i++)
	{
		m_Y[i]=(m_Y[i]-LINE_OFF)/LINE_SCALE;
		m_X[i]=(m_X[i]-SAMP_OFF)/SAMP_SCALE;
		
		m_P[i]=(m_P[i]-LAT_OFF)/LAT_SCALE;
		m_L[i]=(m_L[i]-LONG_OFF)/LONG_SCALE;
		m_H[i]=(m_H[i]-HEIGHT_OFF)/HEIGHT_SCALE;		
	}


	// computer the value of the rpc model
	IDoRPC3();

	m_GlobalAffineModel.createAffineTranslationModel(this, SAMP_OFF, LINE_OFF, 0, SAMP_SCALE, LINE_SCALE);

	return true;
}



//##ModelId=4E4A89B3012D
void RSRPCModel::DoRPC1()//分母不同，阶数为1
{
	double a[7],aa[49],al[7];
	double bb[49],bl[7];
	
	memset(a,0,sizeof(double)*7);
	memset(aa,0,sizeof(double)*49);
	memset(al,0,sizeof(double)*7);

	memset(bb,0,sizeof(double)*49);
	memset(bl,0,sizeof(double)*7);
	
	for(int i=0;i<m_nRPCgcpNum;i++)
	{
		a[0]=1.0;a[1]=m_L[i];a[2]=m_P[i];a[3]=m_H[i];
		a[4]=-m_Y[i]*m_L[i];a[5]=-m_Y[i]*m_P[i];a[6]=-m_Y[i]*m_H[i];
		helper.pNormal(a,7,m_Y[i],aa,al,1.0);		
		
		a[0]=1.0;a[1]=m_L[i];a[2]=m_P[i];a[3]=m_H[i];
		a[4]=-m_X[i]*m_L[i];a[5]=-m_X[i]*m_P[i];a[6]=-m_X[i]*m_H[i];
		helper.pNormal(a,7,m_X[i],bb,bl,1.0);	
	}	
	helper.Gauss(aa,al,7);
	helper.Gauss(bb,bl,7);

	memset(LNUM,0,sizeof(double)*20);
	memset(LDEN,0,sizeof(double)*20);
	memset(SNUM,0,sizeof(double)*20);
	memset(SDEN,0,sizeof(double)*20);
	
	memcpy(LNUM,al,sizeof(double)*4);
	memcpy(LDEN+1,al+4,sizeof(double)*3);
	memcpy(SNUM,bl,sizeof(double)*4);
	memcpy(SDEN+1,bl+4,sizeof(double)*3);
	
	LDEN[0]=1.0; 	
	SDEN[0]=1.0; 
}


//##ModelId=4E4A89B30138
void RSRPCModel::DoRPC2()//分母不同，阶数为2
{	
	DoRPC1();
	double a[19];
	double aa[361],al[19];
	double bb[361],bl[19];	
	double altl=0,bltl=0;

	memset(a,0,sizeof(double)*19);
	memset(aa,0,sizeof(double)*361);
	memset(al,0,sizeof(double)*19);	
	
	memset(bb,0,sizeof(double)*361);
	memset(bl,0,sizeof(double)*19);
    int j;
	for(int i=0;i<m_nRPCgcpNum;i++)
	{	
		a[0]=1.0;              
		a[1]=m_L[i];               
		a[2]=m_P[i];               
		a[3]=m_H[i];               
		a[4]=m_L[i]*m_P[i];        
		a[5]=m_L[i]*m_H[i];        
		a[6]=m_P[i]*m_H[i];        
		a[7]=m_L[i]*m_L[i];        
		a[8]=m_P[i]*m_P[i];        
		a[9]=m_H[i]*m_H[i];        
		
		for(j=10;j<=18;j++)
		{
			a[j]=-a[j-9]*m_Y[i];
		}
		
		helper.pNormal(a,19,m_Y[i],aa,al,1.0);	

		altl+=m_Y[i];
		

		a[0]=1.0;              
		a[1]=m_L[i];               
		a[2]=m_P[i];               
		a[3]=m_H[i];               
		a[4]=m_L[i]*m_P[i];        
		a[5]=m_L[i]*m_H[i];        
		a[6]=m_P[i]*m_H[i];        
		a[7]=m_L[i]*m_L[i];        
		a[8]=m_P[i]*m_P[i];        
		a[9]=m_H[i]*m_H[i];        

		for(j=10;j<=18;j++)
		{
			a[j]=-a[j-9]*m_X[i];
		}
	
		helper.pNormal(a,19,m_X[i],bb,bl,1.0);
		
		bltl+=m_X[i];
	}
	

    double x[39];
	memcpy(x,LNUM,sizeof(double)*10);
	memcpy(x+10,LDEN+1,sizeof(double)*9);
	helper.GaussExt(aa, al, x, 19);
//	GaussL(aa,al,x,19,altl);
	memcpy(al,x,sizeof(double)*19);

	memcpy(x,SNUM,sizeof(double)*10);
	memcpy(x+10,SDEN+1,sizeof(double)*9);
	helper.GaussExt(bb,bl,x,19);
	//GaussL(bb,bl,x,19,bltl);;
	memcpy(bl,x,sizeof(double)*19);


	memset(LNUM,0,sizeof(double)*20);
	memset(LDEN,0,sizeof(double)*20);
	memset(SNUM,0,sizeof(double)*20);
	memset(SDEN,0,sizeof(double)*20);
	
	memcpy(LNUM,al,sizeof(double)*10);
	memcpy(LDEN+1,al+10,sizeof(double)*9);
	memcpy(SNUM,bl,sizeof(double)*10);
	memcpy(SDEN+1,bl+10,sizeof(double)*9);
	
	LDEN[0]=1.0; 	
	SDEN[0]=1.0; 
}


//##ModelId=4E4A89B30177
void RSRPCModel::IDoRPC3()//分母不同，阶数为3
{
	DoRPC3();


	double a[39],aa[1521],al[39];
	double bb[1521],bl[39];	
	 double aext[20],NumL,DenL,NumS,DenS,P,L,H;
	int j;
	double x[39];

	double sigema0=1e10,sigema=1e10;


	do {

		sigema0=sigema;
		memset(a,0,sizeof(double)*39);
		memset(aa,0,sizeof(double)*1521);
		memset(al,0,sizeof(double)*39);			
		memset(bb,0,sizeof(double)*1521);
		memset(bl,0,sizeof(double)*39);
		
		sigema=0;
		for(int i=0;i<m_nRPCgcpNum;i++)
		{
			P=m_P[i],L=m_L[i],H=m_H[i];
			RPCPLH2a(P,L,H,aext);
			NumL=helper.RPCPLH(LNUM,P,L,H);	
			DenL=helper.RPCPLH(LDEN,P,L,H);	
			NumS=helper.RPCPLH(SNUM,P,L,H);	
			DenS= helper.RPCPLH(SDEN,P,L,H);	
			
			for(j=0;j<20;j++)
			{
				a[j]=aext[j]/DenL;
			}
			for(j=20;j<39;j++)
			{
				a[j]=-aext[j-19]*NumL/DenL/DenL;
			}			
			
			helper.pNormal(a,39,m_Y[i]-NumL/DenL,aa,al,1.0);
			
			sigema+=(m_Y[i]-NumL/DenL)*(m_Y[i]-NumL/DenL);
			
			
			for(j=0;j<20;j++)
			{
				a[j]=aext[j]/DenS;
			}
			for(j=20;j<39;j++)
			{
				a[j]=-aext[j-19]*NumS/DenS/DenS;
			}			
			
			helper.pNormal(a,39,m_X[i]-NumS/DenS,bb,bl,1.0);
			sigema+=(m_X[i]-NumS/DenS)*(m_X[i]-NumS/DenS);
		}
		sigema/=m_nRPCgcpNum;
		if(sigema>sigema0) break;


		memset(x,0,sizeof(double)*39);
		helper.GaussExt(aa,al,x,39);

//		for(j=0;j<39;j++)
//		{
//				aa[j*39+j]+=0.002;
//		}	
//		Gauss(aa,al,39);;	
//		memcpy(x,al,sizeof(double)*39);
		
		for(j=0;j<20;j++)
		{
			LNUM[j]+=x[j];
		}
		for(j=20;j<39;j++)
		{
			LDEN[j-19]+=x[j];
		}

		memset(x,0,sizeof(double)*39);
		helper.GaussExt(bb,bl,x,39);;	

//		for(j=0;j<39;j++)
//		{
//				bb[j*39+j]+=0.002;
//		}	
//		Gauss(bb,bl,39);;	
//		memcpy(x,bl,sizeof(double)*39);


		for(j=0;j<20;j++)
		{
			SNUM[j]+=x[j];
		}
		for(j=20;j<39;j++)
		{
			SDEN[j-19]+=x[j];
		}	
		LDEN[0]=1.0; 	
		SDEN[0]=1.0; 

		break;
		
	} while(sigema<sigema0);
}
//##ModelId=4E4A89B30148
void RSRPCModel::DoRPC3()//分母不同，阶数为3
{

	DoRPC2();
//
//	FILE *fp=fopen("d:\\a.txt","w");
//	fprintf(fp,"%d\n",m_nRPCgcpNum);
	double a[39],aa[1521],al[39];
	double bb[1521],bl[39];	
	memset(a,0,sizeof(double)*39);
	memset(aa,0,sizeof(double)*1521);
	memset(al,0,sizeof(double)*39);			
	memset(bb,0,sizeof(double)*1521);
	memset(bl,0,sizeof(double)*39);
	double altl=0,bltl=0;
	int j;
	for(int i=0;i<m_nRPCgcpNum;i++)
	{	
		a[0]=1.0;              
		a[1]=m_L[i];                    
		a[2]=m_P[i];                    
		a[3]=m_H[i];                    
		a[4]=m_L[i]*m_P[i];             
		a[5]=m_L[i]*m_H[i];             
		a[6]=m_P[i]*m_H[i];             
		a[7]=m_L[i]*m_L[i];             
		a[8]=m_P[i]*m_P[i];             
		a[9]=m_H[i]*m_H[i];            	
		a[10]=m_P[i]*m_L[i]*m_H[i];
		a[11]=m_L[i]*m_L[i]*m_L[i];
		a[12]=m_L[i]*m_P[i]*m_P[i];
		a[13]=m_L[i]*m_H[i]*m_H[i];
		a[14]=m_L[i]*m_L[i]*m_P[i];
		a[15]=m_P[i]*m_P[i]*m_P[i];
		a[16]=m_P[i]*m_H[i]*m_H[i];
		a[17]=m_L[i]*m_L[i]*m_H[i];
		a[18]=m_P[i]*m_P[i]*m_H[i];
		a[19]=m_H[i]*m_H[i]*m_H[i];

		for(j=20;j<39;j++)
		{
			a[j]=-a[j-19]*m_Y[i];
		}

//		for(j=0;j<39;j++)
//		{
//			fprintf(fp,"%20.15lf ",a[j]);
//		}
//		
//		fprintf(fp,"%20.15lf \n",m_Y[i]);
		helper.pNormal(a,39,m_Y[i],aa,al,1.0);	

		altl+=m_Y[i];

		a[0]=1.0;              
		a[1]=m_L[i];                    
		a[2]=m_P[i];                    
		a[3]=m_H[i];                    
		a[4]=m_L[i]*m_P[i];             
		a[5]=m_L[i]*m_H[i];             
		a[6]=m_P[i]*m_H[i];             
		a[7]=m_L[i]*m_L[i];             
		a[8]=m_P[i]*m_P[i];             
		a[9]=m_H[i]*m_H[i];            	
		a[10]=m_P[i]*m_L[i]*m_H[i];
		a[11]=m_L[i]*m_L[i]*m_L[i];
		a[12]=m_L[i]*m_P[i]*m_P[i];
		a[13]=m_L[i]*m_H[i]*m_H[i];
		a[14]=m_L[i]*m_L[i]*m_P[i];
		a[15]=m_P[i]*m_P[i]*m_P[i];
		a[16]=m_P[i]*m_H[i]*m_H[i];
		a[17]=m_L[i]*m_L[i]*m_H[i];
		a[18]=m_P[i]*m_P[i]*m_H[i];
		a[19]=m_H[i]*m_H[i]*m_H[i];		
		
		for(j=20;j<39;j++)
		{
			a[j]=-a[j-19]*m_X[i];
		}		
		
		helper.pNormal(a,39,m_X[i],bb,bl,1.0);	

		bltl+=m_X[i];
	}

	
	double x[39];
	memcpy(x,LNUM,sizeof(double)*20);
	memcpy(x+20,LDEN+1,sizeof(double)*19);
//	for(j=0;j<39;j++)
//	{
//		aa[j*39+j]+=0.001;
//	}	
//	Gauss(aa,al,39);;	

	helper.GaussExt(aa,al,x,39);;	
    memcpy(al,x,sizeof(double)*39);

	memcpy(x,SNUM,sizeof(double)*20);
	memcpy(x+20,SDEN+1,sizeof(double)*19);

//	for(j=0;j<39;j++)
//	{
//		bb[j*39+j]+=0.001;
//	}
//	Gauss(bb,bl,39);;	
	helper.GaussExt(bb,bl,x,39);;	
	memcpy(bl,x,sizeof(double)*39);
	
	memset(LNUM,0,sizeof(double)*20);
	memset(LDEN,0,sizeof(double)*20);
	memset(SNUM,0,sizeof(double)*20);
	memset(SDEN,0,sizeof(double)*20);
	
	memcpy(LNUM,al,sizeof(double)*20);
	memcpy(LDEN+1,al+20,sizeof(double)*19);
	memcpy(SNUM,bl,sizeof(double)*20);
	memcpy(SDEN+1,bl+20,sizeof(double)*19);
	
	LDEN[0]=1.0; 	
	SDEN[0]=1.0; 


}

//##ModelId=4E4A89B30149
void RSRPCModel::DoRPC4()//分母相同不为1，阶数为1
{
	
	double a[11],aa[121],al[11];
	
	memset(a,0,sizeof(double)*11);
	memset(aa,0,sizeof(double)*121);
	memset(al,0,sizeof(double)*11);
	
	for(int i=0;i<m_nRPCgcpNum;i++)
	{	
		memset(a,0,sizeof(double)*11);
		
		a[0]=1.0;a[1]=m_L[i];a[2]=m_P[i];a[3]=m_H[i];
		a[4]=-m_Y[i]*m_L[i];a[5]=-m_Y[i]*m_P[i];a[6]=-m_Y[i]*m_H[i];
		helper.pNormal(a,11,m_Y[i],aa,al,1.0);

		memset(a,0,sizeof(double)*11);		
		a[4]=-m_X[i]*m_L[i];a[5]=-m_X[i]*m_P[i];a[6]=-m_X[i]*m_H[i];
		a[7]=1.0;a[8]=m_L[i];a[9]=m_P[i];a[10]=m_H[i];
		helper.pNormal(a,11,m_X[i],aa,al,1.0);	
	}	
	helper.Gauss(aa,al,11);
	
	memset(LNUM,0,sizeof(double)*20);
	memset(LDEN,0,sizeof(double)*20);
	memset(SNUM,0,sizeof(double)*20);
	memset(SDEN,0,sizeof(double)*20);
	
	memcpy(LNUM,al,sizeof(double)*4);
	memcpy(LDEN+1,al+4,sizeof(double)*3);
	memcpy(SNUM,al+7,sizeof(double)*4);
	memcpy(SDEN+1,al+4,sizeof(double)*3);
	
	LDEN[0]=1.0; 	
	SDEN[0]=1.0; 


}
//##ModelId=4E4A89B30157
void  RSRPCModel::DoRPC5()//分母相同不为1，阶数为2
{
	DoRPC4();
	double a[29],aa[841],al[29];
	
	memset(a,0,sizeof(double)*29);
	memset(aa,0,sizeof(double)*841);
	memset(al,0,sizeof(double)*29);
	int j;
	for(int i=0;i<m_nRPCgcpNum;i++)
	{	
		memset(a,0,sizeof(double)*29);

		a[0]=1.0;              
		a[1]=m_L[i];               
		a[2]=m_P[i];               
		a[3]=m_H[i];               
		a[4]=m_L[i]*m_P[i];        
		a[5]=m_L[i]*m_H[i];        
		a[6]=m_P[i]*m_H[i];        
		a[7]=m_L[i]*m_L[i];        
		a[8]=m_P[i]*m_P[i];        
		a[9]=m_H[i]*m_H[i];        
		
		for(j=10;j<=18;j++)
		{
			a[j]=-a[j-9]*m_Y[i];
		}
		
		helper.pNormal(a,29,m_Y[i],aa,al,1.0);
		
		memset(a,0,sizeof(double)*29);
		
		a[19]=1.0;              
		a[20]=m_L[i];               
		a[21]=m_P[i];               
		a[22]=m_H[i];               
		a[23]=m_L[i]*m_P[i];        
		a[24]=m_L[i]*m_H[i];        
		a[25]=m_P[i]*m_H[i];        
		a[26]=m_L[i]*m_L[i];        
		a[27]=m_P[i]*m_P[i];        
		a[28]=m_H[i]*m_H[i];        
		for(j=10;j<=18;j++)
		{
			a[j]=-a[j+10]*m_X[i];
		}		
	
		helper.pNormal(a,29,m_X[i],aa,al,1.0);	
	}

	double x[29];


	memcpy(x,LNUM,sizeof(double)*10);
	memcpy(x+10,LDEN+1,sizeof(double)*9);
	memcpy(x+19,SNUM,sizeof(double)*10);
	memcpy(x+10,SDEN+1,sizeof(double)*9);


	helper.GaussExt(aa,al,x,29);
	memcpy(al,x,sizeof(double)*29);
	
	memcpy(LNUM,al,sizeof(double)*10);
	memcpy(LDEN+1,al+10,sizeof(double)*9);
	memcpy(SNUM,al+19,sizeof(double)*10);
	memcpy(SDEN+1,al+10,sizeof(double)*9);
	
	LDEN[0]=1.0; 	
	SDEN[0]=1.0; 
}



//##ModelId=4E4A89B30158
void RSRPCModel::DoRPC6()//分母相同不为1，阶数为3
{
	DoRPC5();
	double a[59],aa[3481],al[59];
	
	memset(a,0,sizeof(double)*59);
	memset(aa,0,sizeof(double)*3481);
	memset(al,0,sizeof(double)*59);
	int j;
	for(int i=0;i<m_nRPCgcpNum;i++)
	{	
		memset(a,0,sizeof(double)*59);
		
		a[0]=1.0;              
		a[1]=m_L[i];                    
		a[2]=m_P[i];                    
		a[3]=m_H[i];                    
		a[4]=m_L[i]*m_P[i];             
		a[5]=m_L[i]*m_H[i];             
		a[6]=m_P[i]*m_H[i];             
		a[7]=m_L[i]*m_L[i];             
		a[8]=m_P[i]*m_P[i];             
		a[9]=m_H[i]*m_H[i];            	
		a[10]=m_P[i]*m_L[i]*m_H[i];
		a[11]=m_L[i]*m_L[i]*m_L[i];
		a[12]=m_L[i]*m_P[i]*m_P[i];
		a[13]=m_L[i]*m_H[i]*m_H[i];
		a[14]=m_L[i]*m_L[i]*m_P[i];
		a[15]=m_P[i]*m_P[i]*m_P[i];
		a[16]=m_P[i]*m_H[i]*m_H[i];
		a[17]=m_L[i]*m_L[i]*m_H[i];
		a[18]=m_P[i]*m_P[i]*m_H[i];
		a[19]=m_H[i]*m_H[i]*m_H[i];    
		
		for(j=20;j<=38;j++)
		{
			a[j]=-a[j-19]*m_Y[i];
		}
		
		helper.pNormal(a,59,m_Y[i],aa,al,1.0);
		
		memset(a,0,sizeof(double)*59);
		
		a[39]=1.0;              
		a[40]=m_L[i];                    
		a[41]=m_P[i];                    
		a[42]=m_H[i];                    
		a[43]=m_L[i]*m_P[i];             
		a[44]=m_L[i]*m_H[i];             
		a[45]=m_P[i]*m_H[i];             
		a[46]=m_L[i]*m_L[i];             
		a[47]=m_P[i]*m_P[i];             
		a[48]=m_H[i]*m_H[i];            	
		a[49]=m_P[i]*m_L[i]*m_H[i];
		a[50]=m_L[i]*m_L[i]*m_L[i];
		a[51]=m_L[i]*m_P[i]*m_P[i];
		a[52]=m_L[i]*m_H[i]*m_H[i];
		a[53]=m_L[i]*m_L[i]*m_P[i];
		a[54]=m_P[i]*m_P[i]*m_P[i];
		a[55]=m_P[i]*m_H[i]*m_H[i];
		a[56]=m_L[i]*m_L[i]*m_H[i];
		a[57]=m_P[i]*m_P[i]*m_H[i];
		a[58]=m_H[i]*m_H[i]*m_H[i]; 
		
		for(j=20;j<=38;j++)
		{
			a[j]=-a[j+20]*m_X[i];
		}		
		
		helper.pNormal(a,59,m_X[i],aa,al,1.0);	
	}

	double x[59];


	memcpy(x,LNUM,sizeof(double)*20);
	memcpy(x+20,LDEN+1,sizeof(double)*19);
	memcpy(x+39,SNUM,sizeof(double)*20);
	memcpy(x+20,SDEN+1,sizeof(double)*19);	
	
	helper.GaussExt(aa, al, x, 59);
	memcpy(al,x,sizeof(double)*59);	
	
	memcpy(LNUM,al,sizeof(double)*20);
	memcpy(LDEN+1,al+20,sizeof(double)*19);
	memcpy(SNUM,al+39,sizeof(double)*20);
	memcpy(SDEN+1,al+20,sizeof(double)*19);
	
	LDEN[0]=1.0; 	
	SDEN[0]=1.0; 
}


//##ModelId=4E4A89B30167
void RSRPCModel::DoRPC7()//无分母，阶数为1
{
	double a[4],aa[16],al[4];
	
	memset(a,0,sizeof(double)*4);
	memset(aa,0,sizeof(double)*16);
	memset(al,0,sizeof(double)*4);	
	
	double bb[16],bl[4];
	
	
	memset(bb,0,sizeof(double)*16);
	memset(bl,0,sizeof(double)*4);
	
	for(int i=0;i<m_nRPCgcpNum;i++)
	{		
		a[0]=1.0;a[1]=m_L[i];a[2]=m_P[i];a[3]=m_H[i];
		
		helper.pNormal(a,4,m_Y[i],aa,al,1.0);		
		
		a[0]=1.0;a[1]=m_L[i];a[2]=m_P[i];a[3]=m_H[i];
		
		helper.pNormal(a,4,m_X[i],bb,bl,1.0);	
	}	
	helper.Gauss(aa,al,4);
	helper.Gauss(bb,bl,4);
	
	memset(LNUM,0,sizeof(double)*20);
	memset(LDEN,0,sizeof(double)*20);
	memset(SNUM,0,sizeof(double)*20);
	memset(SDEN,0,sizeof(double)*20);
	
	memcpy(LNUM,al,sizeof(double)*4);	
	memcpy(SNUM,bl,sizeof(double)*4);
	
	
	LDEN[0]=1.0; 	
	SDEN[0]=1.0; 

}
//##ModelId=4E4A89B30168
void RSRPCModel::DoRPC8()//无分母，阶数为2
{
	DoRPC7();
	double a[10],aa[100],al[10];
	double bb[100],bl[10];	
	memset(a,0,sizeof(double)*10);
	memset(aa,0,sizeof(double)*100);
	memset(al,0,sizeof(double)*10);			
	memset(bb,0,sizeof(double)*100);
	memset(bl,0,sizeof(double)*10);
	
	for(int i=0;i<m_nRPCgcpNum;i++)
	{		
		a[0]=1.0;              
		a[1]=m_L[i];               
		a[2]=m_P[i];               
		a[3]=m_H[i];               
		a[4]=m_L[i]*m_P[i];        
		a[5]=m_L[i]*m_H[i];        
		a[6]=m_P[i]*m_H[i];        
		a[7]=m_L[i]*m_L[i];        
		a[8]=m_P[i]*m_P[i];        
		a[9]=m_H[i]*m_H[i];        
		
		helper.pNormal(a,10,m_Y[i],aa,al,1.0);	
		
		
		a[0]=1.0;              
		a[1]=m_L[i];               
		a[2]=m_P[i];               
		a[3]=m_H[i];               
		a[4]=m_L[i]*m_P[i];        
		a[5]=m_L[i]*m_H[i];        
		a[6]=m_P[i]*m_H[i];        
		a[7]=m_L[i]*m_L[i];        
		a[8]=m_P[i]*m_P[i];        
		a[9]=m_H[i]*m_H[i];        

		helper.pNormal(a, 10, m_X[i], bb, bl, 1.0);
	}
	
	helper.Gauss(aa,al,10);
	helper.Gauss(bb, bl, 10);
	
	
	memset(LNUM,0,sizeof(double)*20);
	memset(LDEN,0,sizeof(double)*20);
	memset(SNUM,0,sizeof(double)*20);
	memset(SDEN,0,sizeof(double)*20);
	
	memcpy(LNUM,al,sizeof(double)*10);	
	memcpy(SNUM,bl,sizeof(double)*10);
	
	
	LDEN[0]=1.0; 	
	SDEN[0]=1.0; 
}


//##ModelId=4E4A89B30169
void RSRPCModel::DoRPC9()//无分母，阶数为3
{
	DoRPC8();
	double a[20],aa[400],al[20];
	double bb[400],bl[20];	
	memset(a,0,sizeof(double)*20);
	memset(aa,0,sizeof(double)*400);
	memset(al,0,sizeof(double)*20);			
	memset(bb,0,sizeof(double)*400);
	memset(bl,0,sizeof(double)*20);
	
	for(int i=0;i<m_nRPCgcpNum;i++)
	{		
		a[0]=1.0;              
		a[1]=m_L[i];                    
		a[2]=m_P[i];                    
		a[3]=m_H[i];                    
		a[4]=m_L[i]*m_P[i];             
		a[5]=m_L[i]*m_H[i];             
		a[6]=m_P[i]*m_H[i];             
		a[7]=m_L[i]*m_L[i];             
		a[8]=m_P[i]*m_P[i];             
		a[9]=m_H[i]*m_H[i];            	
		a[10]=m_P[i]*m_L[i]*m_H[i];
		a[11]=m_L[i]*m_L[i]*m_L[i];
		a[12]=m_L[i]*m_P[i]*m_P[i];
		a[13]=m_L[i]*m_H[i]*m_H[i];
		a[14]=m_L[i]*m_L[i]*m_P[i];
		a[15]=m_P[i]*m_P[i]*m_P[i];
		a[16]=m_P[i]*m_H[i]*m_H[i];
		a[17]=m_L[i]*m_L[i]*m_H[i];
		a[18]=m_P[i]*m_P[i]*m_H[i];
		a[19]=m_H[i]*m_H[i]*m_H[i];
		
	
		
		helper.pNormal(a, 20, m_Y[i], aa, al, 1.0);
		
		a[0]=1.0;              
		a[1]=m_L[i];                    
		a[2]=m_P[i];                    
		a[3]=m_H[i];                    
		a[4]=m_L[i]*m_P[i];             
		a[5]=m_L[i]*m_H[i];             
		a[6]=m_P[i]*m_H[i];             
		a[7]=m_L[i]*m_L[i];             
		a[8]=m_P[i]*m_P[i];             
		a[9]=m_H[i]*m_H[i];            	
		a[10]=m_P[i]*m_L[i]*m_H[i];
		a[11]=m_L[i]*m_L[i]*m_L[i];
		a[12]=m_L[i]*m_P[i]*m_P[i];
		a[13]=m_L[i]*m_H[i]*m_H[i];
		a[14]=m_L[i]*m_L[i]*m_P[i];
		a[15]=m_P[i]*m_P[i]*m_P[i];
		a[16]=m_P[i]*m_H[i]*m_H[i];
		a[17]=m_L[i]*m_L[i]*m_H[i];
		a[18]=m_P[i]*m_P[i]*m_H[i];
		a[19]=m_H[i]*m_H[i]*m_H[i];				
		
		helper.pNormal(a, 20, m_X[i], bb, bl, 1.0);
	}

    double x[20];
	memcpy(x,LNUM,sizeof(double)*20);
	helper.GaussExt(aa, al, x, 20);
	memcpy(al,x,sizeof(double)*20);	
	memset(LDEN,0,sizeof(double)*20);
	memcpy(LNUM,al,sizeof(double)*20);	

	memcpy(x,SNUM,sizeof(double)*20);
	helper.GaussExt(bb, bl, x, 20);
	memcpy(bl,x,sizeof(double)*20);	
	memcpy(SNUM,bl,sizeof(double)*20);	
	
	LDEN[0]=1.0; 	
	SDEN[0]=1.0; 
}


//##ModelId=4E4A89B30119
void RSRPCModel::RPCPLH2a(double P,double L,double H,double *a)
{	
	a[0]=1;a[1]=L;a[2]=P;a[3]=H;a[4]=L*P;a[5]=L*H;a[6]=P*H;a[7]=L*L;a[8]=P*P;a[9]=H*H;
	a[10]=P*L*H;a[11]=L*L*L;a[12]=L*P*P;a[13]=L*H*H;a[14]=L*L*P;a[15]=P*P*P;a[16]=P*H*H;
	a[17]=L*L*H;a[18]=P*P*H;a[19]=H*H*H;	
}


void RSRPCModel::RPCError( const char* pPath )
{	
	double m_Minx         ,m_Miny;
	double m_ImgHeight    ,m_ImgWidth;
	m_ImgHeight =m_ImgRange.bottom - m_ImgRange.top;
	m_ImgWidth  =m_ImgRange.right  - m_ImgRange.left;
	double dz = (m_MaxElevation-m_MinElevation)/(m_nz-1);
	m_Minx = m_ImgRange.left;
	m_Miny = m_ImgRange.top;

	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];
    _splitpath_s( pPath,drive, dir, fname, ext );
	char timestring[20];
	_snprintf_s(timestring, 200, "%s_%02d_%02d_%02d_%0d.txt", fname,m_nx,m_ny,m_nz, 3);	// 选用何种平差模型，其最后一个参数应当不一样的
	CString path;
	path.Format("%s%s",drive,dir);


//	RSConfigure &Config = RSGetConfigure();
//	CString path=Config.GetTempFilePath();	
    FILE *fp;
	fopen_s(&fp,path+"\\control"+timestring,"w");

	double xmax,xmin,xv;
	double ymax,ymin,yv;
	double xymax,xymin,xyv;

	fprintf(fp,"%10d\n",m_nRPCgcpNum);
	m_nRPCgcpNum = m_nx * m_ny * m_nz;
	double x,y,xx,yy;
	double* xxxx1=new double[m_nRPCgcpNum];
	double* xxxx2=new double[m_nRPCgcpNum];
	double* xxxx3=new double[m_nRPCgcpNum];

	// 代码有问题，因为这里的x,y肯定在影像范围内
	int num=0;
	int ii,i,j;
	for(ii=0;ii<m_nz;ii++)
	{
		for(i=0;i<m_ny;i++)
		{
			for(j=0;j<m_nx;j++)
			{
				// 这时候检查控制点精度
				FromGroundToImage(m_RPCgcp[ii].ddx[i*m_nx+j],m_RPCgcp[ii].ddy[i*m_nx+j],m_RPCgcp[ii].h,xx,yy);
				x=m_RPCgcp[ii].x0+j*m_RPCgcp[ii].dx;
				y=m_RPCgcp[ii].y0+i*m_RPCgcp[ii].dy;

				DPoint2D trpo(x,y);
 				if ( !m_ModelRange.InPolygon( trpo ) ) continue;


				xxxx1[num]=x-xx;
				xxxx2[num]=y-yy;
				xxxx3[num]=sqrt(xxxx1[num]*xxxx1[num]+xxxx2[num]*xxxx2[num]);
#ifdef _DEBUG
				if ( xxxx3[num] > 0.3)
				{
					double x1,y1,lat1,lon1;
					((RSGeometryModel*)m_pOrignalModel)->FromImageToGround(x, y, m_RPCgcp[ii].h, lat1, lon1);
					((RSGeometryModel*)m_pOrignalModel)->FromGroundToImage(m_RPCgcp[ii].ddx[i*m_nx+j],m_RPCgcp[ii].ddy[i*m_nx+j],
						m_RPCgcp[ii].h,x1,y1);
					FromGroundToImage(m_RPCgcp[ii].ddx[i*m_nx+j],m_RPCgcp[ii].ddy[i*m_nx+j],m_RPCgcp[ii].h,xx,yy);
					m_RPCgcp[ii].ddx[i*m_nx+j]=x-xx;
					m_RPCgcp[ii].ddy[i*m_nx+j]=y-yy;
				}
#endif				

				m_RPCgcp[ii].ddx[i*m_nx+j]=x-xx;
				m_RPCgcp[ii].ddy[i*m_nx+j]=y-yy;


				if(x>m_ImgWidth) continue;
				if(x<0) continue;
				if(y>m_ImgHeight) continue;
				if(y<0) continue;

				
				fprintf(fp,"%10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf\n",x,y,m_RPCgcp[ii].h,xxxx1[num],xxxx2[num],xxxx3[num]);
				num++;
			}
		}			
	}	



	helper.GetMaxMinRMS(xxxx1, num, xmax, xmin, xv);
	helper.GetMaxMinRMS(xxxx2, num, ymax, ymin, yv);
	helper.GetMaxMinRMS(xxxx3, num, xymax, xymin, xyv);


	fprintf(fp,"min%15.8lf %15.8lf %15.8lf\n",xmin,ymin,xymin);
	fprintf(fp,"max%15.8lf %15.8lf %15.8lf\n",xmax,ymax,xymax);
	fprintf(fp,"  v%15.8lf %15.8lf %15.8lf\n",xv,yv,xyv);
	fclose(fp);

	fopen_s(&fp,path+"\\"+fname+".dat","wb");
	fwrite(&m_nz,1,sizeof(int),fp);
	for(i=0;i<m_nz;i++)
	{
		fwrite(&m_RPCgcp[i].h,1,sizeof(double),fp);
		fwrite(&m_RPCgcp[i].x0,1,sizeof(double),fp);
		fwrite(&m_RPCgcp[i].y0,1,sizeof(double),fp);
		fwrite(&m_RPCgcp[i].dx,1,sizeof(double),fp);
		fwrite(&m_RPCgcp[i].dy,1,sizeof(double),fp);
		fwrite(&m_RPCgcp[i].nx,1,sizeof(int),fp);
		fwrite(&m_RPCgcp[i].ny,1,sizeof(int),fp);
		fwrite(m_RPCgcp[i].ddx,1,sizeof(double)*m_RPCgcp[i].nx*m_RPCgcp[i].ny,fp);
		fwrite(m_RPCgcp[i].ddy,1,sizeof(double)*m_RPCgcp[i].nx*m_RPCgcp[i].ny,fp);
	}
	fclose(fp);

	//获得检查点
	double checkblockwidth=m_RPCgcp[0].dx;
	double checkblockheight=m_RPCgcp[0].dy;
	double checkdheight=m_RPCgcp[1].h-m_RPCgcp[0].h;

	double DemH=m_MinElevation+checkdheight/2.0;
	long checkpointnum=(m_nx-1)*(m_ny-1)*(m_nz-1);

	fopen_s(&fp,path+"\\check"+timestring,"w");

	fprintf(fp,"%10d\n",checkpointnum);
    
	num=0;
	double lat,lon;
	for(ii=0;ii<m_nz-1;ii++)
	{	
		for(i=0;i<m_ny-1;i++)
		{	
			y=m_Miny+checkblockheight/2.0+i*checkblockheight;		
			for(j=0;j<m_nx-1;j++)
			{
				double ddx,ddy;
				x=m_Minx+checkblockwidth/2.0+j*checkblockwidth;	


				DPoint2D trpo(x,y),orpo;
				DPoint3D gpo;
				if ( !m_ModelRange.InPolygon( trpo ) ) continue;
				((ImageMappingModel*)m_pImageTrans)->TransformedToOriginal( trpo.x, trpo.y, orpo.x, orpo.y );
				((RSGeometryModel*)m_pOrignalModel)->FromImageToGround(orpo.x, orpo.y, DemH, lon, lat);
				((ObjectMappingModel*)m_pGroundTrans)->OriginalToTransformed( lon, lat, DemH, gpo.lon, gpo.lat, gpo.h);
				FromGroundToImage(gpo.lon,gpo.lat,gpo.h,xx,yy);

                int minhIndex,maxhIndex;
				if(DemH<m_MinElevation) 
				{
					minhIndex=maxhIndex=0;
				}else if(DemH>m_MaxElevation) 
				{
					minhIndex=maxhIndex=m_nz-1;
				}else
				{
					for(int iii=0;iii<m_nz-1;iii++)
					{
						if(DemH>=m_RPCgcp[iii].h&&DemH<m_RPCgcp[iii+1].h)
						{
							minhIndex=iii;
							maxhIndex=iii+1;
						}
					}
				}
				if(minhIndex==maxhIndex)
				{
					CGridInterpolate gi;
					gi.Init(m_RPCgcp[minhIndex]);
					gi.GetZ(x,y,ddx,ddy);
				}
				else
				{
					double ddx1,ddx2,ddy1,ddy2;
					CGridInterpolate gi1;
					gi1.Init(m_RPCgcp[minhIndex]);
					gi1.GetZ(x,y,ddx1,ddy1);

					CGridInterpolate gi2;
					gi2.Init(m_RPCgcp[maxhIndex]);
					gi2.GetZ(x,y,ddx2,ddy2);

					ddx=ddx1+(DemH-m_RPCgcp[minhIndex].h)/(m_RPCgcp[maxhIndex].h-m_RPCgcp[minhIndex].h)*(ddx2-ddx1);
					ddy=ddy1+(DemH-m_RPCgcp[minhIndex].h)/(m_RPCgcp[maxhIndex].h-m_RPCgcp[minhIndex].h)*(ddy2-ddy1);

				}			
//			    RPCcompensation(rpc,xx,yy,DemH);
				
			    ddx=0;ddy=0;
				xxxx1[num]=x-(xx+ddx);
				xxxx2[num]=y-(yy+ddy);
				xxxx3[num]=sqrt(xxxx1[num]*xxxx1[num]+xxxx2[num]*xxxx2[num]);
				if(x>m_ImgWidth) continue;
				if(x<0) continue;
				if(y>m_ImgHeight) continue;
				if(y<0) continue;
				
				fprintf(fp,"%10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf\n",x,y,DemH,xxxx1[num],xxxx2[num],xxxx3[num]);
				num++; 					
				
			}
		}
		
		DemH+=checkdheight; 
	}

	helper.GetMaxMinRMS(xxxx1, num, xmax, xmin, xv);
	helper.GetMaxMinRMS(xxxx2, num, ymax, ymin, yv);
	helper.GetMaxMinRMS(xxxx3, num, xymax, xymin, xyv);

	
//	fprintf(LogInfo,"Check\n");
//	fprintf(LogInfo, "%15.8lf %15.8lf %15.8lf %15.8lf %15.8lf %15.8lf\n" ,xmax,xv,ymax,yv,xymax,xyv);
//	fprintf(LogInfo,"\n");

	fprintf(fp,"min%15.8lf %15.8lf %15.8lf\n",xmin,ymin,xymin);
	fprintf(fp,"max%15.8lf %15.8lf %15.8lf\n",xmax,ymax,xymax);
	fprintf(fp,"  v%15.8lf %15.8lf %15.8lf\n",xv,yv,xyv);
	fclose(fp);

	delete []xxxx1;
	xxxx1=NULL;
	delete []xxxx2;
	xxxx2=NULL;
	delete []xxxx3;
	xxxx3=NULL;

//	fclose( LogInfo);
}
/*
double RSRPCModel::GetLatOffset() const
{
	return LAT_OFF;
}

double RSRPCModel::GetLonOffset() const
{
	return LONG_OFF;
}

double RSRPCModel::GetHeightOffset() const
{
	return HEIGHT_OFF;
}
*/
//##ModelId=4E4A89B3001F

bool   RSRPCModel::SetAffineParameter(double AffineTrans[12])
{
	memcpy( px, AffineTrans, sizeof(double)*3);
	memcpy( py, AffineTrans+3, sizeof(double)*3);
	memcpy( invpx, AffineTrans+6, sizeof(double)*3);
	memcpy( invpy, AffineTrans+9, sizeof(double)*3);
	return true;
}

