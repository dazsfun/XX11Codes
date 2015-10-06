  // GeoRSDataStruct.h: interface for the CGeoRSDataStruct class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GEORSDATASTRUCT_H__261DDC12_AA2A_47F6_97AD_DAF7D223488C__INCLUDED_)
#define AFX_GEORSDATASTRUCT_H__261DDC12_AA2A_47F6_97AD_DAF7D223488C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <afxtempl.h>
#include "EnumString.h"

#ifndef PI
#define  PI 3.1415926535897932384626433832795
#endif

#ifndef GM
#define    GM 3.9860047e+14
#endif

typedef struct 
{
	double x,y,z;
	
} POINT3D;

enum ZY3LineTimeState
{
	ZeroTime,
	Time_Unconsistence,
	Time_Lacked,
	LineTime_Stable
};

Begin_Enum_String(ZY3LineTimeState)
{
	Enum_String(ZeroTime);
	Enum_String(Time_Unconsistence);
	Enum_String(Time_Lacked);
	Enum_String(LineTime_Stable);
}
End_Enum_String;

enum ZY3OrbitState
{
	ZeroOrbit,
	Orbit_Time_Unconsistence,
	Orbit_Accuracy_Large,
	Orbit_Accuracy_ExceedLimit,
	Orbit_Lacked,
	Orbit_Stable
};

Begin_Enum_String(ZY3OrbitState)
{
	Enum_String(ZeroOrbit);
	Enum_String(Orbit_Time_Unconsistence);
	Enum_String(Orbit_Accuracy_Large);
	Enum_String(Orbit_Accuracy_ExceedLimit);
	Enum_String(Orbit_Lacked);
	Enum_String(Orbit_Stable);
}
End_Enum_String;

enum ZY3AttState
{
	ZeroAtt,
	Att_Time_Unconsistence,
	Att_Accuracy_Large,
	Att_Accuracy_ExceedLimit,
	Att_Lacked,
	Att_Stable
};

Begin_Enum_String(ZY3AttState)
{
	Enum_String(ZeroAtt);
	Enum_String(Att_Time_Unconsistence);
	Enum_String(Att_Accuracy_Large);
	Enum_String(Att_Accuracy_ExceedLimit);
	Enum_String(Att_Lacked);
	Enum_String(Att_Stable);
}
End_Enum_String;

typedef struct
{
	double _AverageTime;
	double _RealUsedAttNum;
	double _RealUsedGpsNum;
	double _PhiErrorRMS;
	double _PhiErrorMax;
	double _OmegaErrorRMS;
	double _OmegaErrorMax;
	double _KappaErrorRMS;
	double _KappaErrorMax;
	double _GpsXErrorRMS;
	double _GpsXErrorMax;
	double _GpsYErrorRMS;
	double _GpsYErrorMax;
	double _GpsZErrorRMS;
	double _GpsZErrorMax;
	double _GpsXvErrorRMS;
	double _GpsXvErrorMax;
	double _GpsYvErrorRMS;
	double _GpsYvErrorMax;
	double _GpsZvErrorRMS;
	double _GpsZvErrorMax;
	double _GpsAvgIntTime;
	double _AttAvgIntTime;

	//细节见态度
	ZY3AttState _attState;
	ZY3OrbitState _orbitState;
	ZY3LineTimeState _LineTState;

	//星下点，如果这个都错了，说明不是建模的问题
	POINT3D _NadirPoint;

}ZYModelAccuracy;

typedef struct ZY3ModelCompensateAccuracy
{
	double SCCirculateErrorXMax;
	double SCCirculateErrorXRms;
	double SCCirculateErrorYMax;
	double SCCirculateErrorYRms;
	double SCRpcComSateErrorXMax;
	double SCRpcComSateErrorXRms;
	double SCRpcComSateErrorYMax;
	double SCRpcComSateErrorYRms;
	double SCOriginCirculationErrorXMax;
	double SCOriginCirculationErrorXRms;
	double SCOriginCirculationErrorYMax;
	double SCOriginCirculationErrorYRms;
};


typedef struct 
{
	double  ls2xxyypx[3];
	double  ls2xxyypy[3];
	CString Leve1ImgName;
	CString RPCfilename;
	CString layover_shadowIndexName;
	CString layoverIndexName;
	CString shadowIndexName;
	long  SARSimupleftLine,SARSimupleftsample;
	CString matchedpoint;//from layoverIndexName 2 Leve1ImgName 
}SARBLOCKINPUT;

typedef struct 
{
	double X_OFF,X_SCALE;
	double Y_OFF,Y_SCALE;
	double Z_OFF,Z_SCALE;
	double LINE_OFF,LINE_SCALE;
	double XNUM[20],XDEN[20];
	double YNUM[20],YDEN[20];
	double ZNUM[20],ZDEN[20];
	double Xv_OFF,Xv_SCALE;
	double Yv_OFF,Yv_SCALE;
	double Zv_OFF,Zv_SCALE;
	double XvNUM[20],XvDEN[20];
	double YvNUM[20],YvDEN[20];
	double ZvNUM[20],ZvDEN[20];
}ORBITRPCMODEL;

typedef struct 
{
	double GEC_LINE_OFF,GEC_LINE_SCALE;	
	double GEC_SAMP_OFF,GEC_SAMP_SCALE;

	double SLC_LINE_OFF,SLC_LINE_SCALE;	
	double SLC_SAMP_OFF,SLC_SAMP_SCALE;

	double LNUM[20],LDEN[20],SNUM[20],SDEN[20];	
}GEC2SLCPIXELTRANSFORMATION;

//SLC_LINE=f1(GEC_LINE,GEC_SAMP);
//SLC_SAMP=f1(GEC_LINE,GEC_SAMP);
//a0+a1*GEC_LINE                  +a2*GEC_SAMP                  +a3*GEC_LINE*GEC_LINE         +a4*GEC_SAMP*GEC_LINE         +a5*GEC_SAMP*GEC_SAMP+
//   a6*GEC_LINE*GEC_LINE*GEC_LINE+a7*GEC_LINE*GEC_LINE*GEC_SAMP+a8*GEC_LINE*GEC_SAMP*GEC_SAMP+a9*GEC_SAMP*GEC_SAMP*GEC_SAMP;

// HH HV VH HH 这个结构最多出现四次
typedef struct 
{
	long linenum;
	long pixelnum;
	long line[20];//存储行计数 =new long[linenum];	
	double *NEBZ[20];//=new doulbe[linenum*pixelnum]
	//linenum最多二十行
	double refp[20]; //=new double[linenum]
	long polyDegree[20]; //=new long[linenum]
	double *coefficient[20];//=new double[linenum*polyDegree]
}TXNEBZ;

typedef struct  
{
	double h;//该层的高程
	double x0,y0;//左上角点坐标
	double dx,dy;//x方向和y方向间隔
	int nx,ny;//间隔数目
	double *ddx,*ddy;//在计算前是lat，lon坐标，计算后是x,y方向残差
	double *lon,*lat;//add by fei 用于主图像保存经纬度
}RPCGCPSTRUCT;

// typedef struct  
// {
// 	double h;//该层的高程
// 	double x0,y0;//左上角点坐标
// 	double dx,dy;//x方向和y方向间隔
// 	int nx,ny;//间隔数目
// 	double *ddx,*ddy;//在计算前是lat，lon坐标，计算后是x,y方向残差
// 	double *lon,*lat;//add by fei 用于主图像保存经纬度
// }RPCGCPSTRUCT;
//辐射校正参数
#define ABSOCALIBRATION  0
#define RELACALIBRATION  1
#define PRIORSTAT        0
#define SELFADPTIVESTAT  1 
//控制点数据结构
typedef enum 

{
	CtrlPt,
	CheckPt,
	InvalidePt,
	GPSCtrlPt,
	GPSCheckPt,
	GPSInvalidePt
} CtrlPtType;

typedef enum
{ 
  DOM,
  HAND, 
  LIB,
  JB3,
  MDA,
  GPS
} GcpSource;
typedef enum 
{
	Nearest_Neighbor,//临近点
	Piecewise_Linear,//双线性
	Cubic_Convolution_4P,//双三次
	Truncated_Sinc_6P,
	Cubic_Convolution_6P,//6点三次		
	Knab_6P,	
	Raised_Cosine_6P,
	Kaiser_SINC_24P,
	Raised_Cosine_16P
} ReSampleMode ;




typedef struct
{
	//坐标,相应参数
	CString ID;//控制点库中的ID//char*
	CString ListID;//文件名＋流水号//char*
    CString GcpInfo;//char*
	double x,y;//象素坐标 左上角为原点	
	double L,B,H;//WGS84 H 椭球高
	double X,Y,Z;//WGS84 地心坐标
	CtrlPtType GcpType;
	GcpSource  GcpSrc;
	double xresidual,yresidual,xyprecision;//XY残差以及XY精度

	double east,north;//平面坐标	
	CString CoordSysID;//平面坐标系标识 WKT结构//char*
	double mxy;//平面精度//char*

	double h;//水准面高
	CString HeightSysID;//原始高程系的EPSG编码//char*
	double mh;//高程精度
	
	LPBYTE lpBlock;//控制点影像//BYTE*
	CString ImgFileName;//char*
	long   ImgWidth;
	long   ImgHeight;
	double Lmc,Bmc;//控制点影像底边中点的经纬度

	//
	double nRctImgPosX, nRctImgPosY;     // the rct point position in the rct view
	double nRefImgPosX, nRefImgPosY;	 // the ref point position in the ref view
	BOOL m_bHaveRctPt, m_bHaveRefPt;	 //the flag to show whether the point in the view has created
	int tagRct, tagRef;				     // the color of the point in the view
	long nUsedCount;
	double sample,line;// added by zhang guo 
	double _L,__L,_B,__B,_H,__H; // added by zhang guo 
	double ddx,ddy;       // added by zhang guo  
	double slcx,slcy;
	double correlation; //add by feiwenbo
	long nBandCount;
	
	double phase;//add by feiwenbo
	double coherence;//add by feiwenbo
	double heightresidual;//add by feiwenbo
	
	//裁切设置 add by feiwenbo
	double subset_x,subset_y;

}GCPSTRUCT;

//J2000与WGS84转换关系
typedef struct ut_instant {
	double	j_date; /* julian decimal date, 0 = 01 Jan -4712 12 HR UT */
	long year; /* year, valid range [-4,712, +2,147,483,647] */
	int	month;	/*	[1-12]	*/
	int	day;	/*	[1-31]	*/
	int	i_hour;	/*	[0-23]	*/
	int	i_minute;	/*	[0-59]	*/
	double	second;	/*	[0-59.9999]	*/
	double	d_hour;		/* [0.0-23.9999] includes minute and second */
	double	d_minute;	/*	[0.0-59.9999] includes second	*/
	int	weekday;	/*	[0-6]	*/
	int	day_of_year;	/*	[1-366]	*/
} UTinstant, * UTinstantPtr;

typedef struct{
	double JD;
	double tai_UTC;
}JDTAIUTC;

typedef struct
{
	double MJD;
	double x,y;//" 地球极移
	double UT1R_UTC,UT1R_TAI;//s
	double dPsi,dEpsilon;//0.001"
}JDUTCUT1;

typedef struct 
{
	int NARG[6];
	double XSIN,XCOS,YSIN,YCOS,UTSIN,UTCOS;
}PolarTide;

typedef struct 
{
	int NARG[6];
	double XSIN,XCOS,YSIN,YCOS;
}PolarLunisolar;
//轨道积分
typedef struct{
	
	double majorAxis;	/* semi-major axis (m)		*/
	double minorAxis;	/* semi-minor axis (m)		*/
	double e;		/* numeric eccentricity 	*/
	double e2;		/* numeric eccentricity squared */
	
	double xOff;		/* offset relative to WGS_84 (m)*/
	double yOff;		/* offset relative to WGS_84 (m)*/
	double zOff;		/* offset relative to WGS_84 (m)*/
	
} DATUM;

typedef struct {
	int n,m;
	double C,S;
	double dc,ds;
}EarthGravityModel;

typedef struct 
{
	double JD;
	double X,Y,Z;
}CelestialBodyPos;


typedef struct 
{
	double latitude;
	double longitude;
	double altitude;
	double X;
	double Y;
	double Z;
	double Xv;
	double Yv;
	double Zv;
	double UT;
	int year,month,day,hour,minute;
	double second;
} POINT_Ep;

typedef struct TimeLine
{
	double times;
	double timeInterval;
}TimeL,*pTimeL;



typedef struct
{
	double a;                       // 长半轴
	double e;                       // 扁率
	double i;                       // 轨道倾角
	double capitalomiga;            // 升交点赤径
	double omiga;                   // 近地点俯角
	double f;                       // 
	double u;                       //
	double t0;                      //
	double r;                       //
	double t;                       //
	double E;                       //
	double M;                       //
} ORBIT;



typedef struct 
{
	double h;
	double min;
	double max;
}HPdensitymodel;

//星历姿态等数据
typedef struct 
{
	double UT;
	double YAW;
	double PITCH;
	double ROLL;
	double vYAW;
	double vPITCH;
	double vROLL;
	double q1,q2,q3,q0;
	double q4;                       // it should be q0,it is used in simulation, add by pan
	double vq1,vq2,vq3;
	double frameID;
	double t;
	int year,month,day,hour,minute;
	double second; //add by zhang guo 
} Attitude;


typedef struct 
{
	double LINE_PERIOD;
	double SCENE_CENTER_TIME;
	double SCENE_CENTER_LINE;
	double SCENE_CENTER_COL;
} TimeStamp;


typedef struct 
{
	double PSI_X;
	double PSI_Y;
} Look_Angle;


typedef struct 
{
	double Xs[6];
	double Ys[6];
	double Zs[6];
	double Xvs[6];
	double Yvs[6];
	double Zvs[6];
	double roll[4];
	double pitch[4];
	double yaw[4];
} OrbitAttitudeModel;


typedef struct 
{
	double Xs[21];
	double Ys[21];
	double Zs[21];
	double Xvs[21];
	double Yvs[21];
	double Zvs[21];
	double roll[21];
	double pitch[21];
	double yaw[21];
	double x[21];

	///Information about accuracy
	double _XsErrorRMS;
	double _XsMaxError;
	double _YsErrorRMS;
	double _YsMaxError;
	double _ZsErrorRMS;
	double _ZsMaxError;

	double _XsVErrorRMS;
	double _XsVMaxError;
	double _YsVErrorRMS;
	double _YsVMaxError;
	double _ZsVErrorRMS;
	double _ZsVMaxError;

	double _RollErrorRMS;
	double _RollMaxError;
	double _PitchErrorRMS;
	double _PitchMaxError;
	double _YawErrorRMS;
	double _YawMaxError;
	
	double _xErrorRMS;
	double _xMaxError;
} OrbitAttitudeRegressionModel;



typedef struct {
	double  x0,y0;
	double  sinA,cosA;
	double  dx,dy;
	long	nx,ny;
	double	*z;
} demGRID;

//订单数据结构

typedef struct
{	
	double r[4];
	double captainangle[4];
	double i[4];
	double u[4];
	double roll[4];
	double pitch[4];
	double yaw[4];
}OrbitPolyModel;

typedef struct
{	
	double r[4];
	double a[4];
	double e[4];
	double captainangle[4];
	double i[4];
	double u[4];
	double f[4];
	double omiga[4];
	double roll[4];
	double pitch[4];
	double yaw[4];
}OrbitPolyModelEx;

typedef struct 
{
	CString OriLevel;//原始影像级别
	CString Level;//校正级别	CString	Level1, Level2, Level3, Level4
	ReSampleMode ResamMode;//重采样内核	CString	"NN"，"Bilinear"，?"CC"
	CString HeightMode,HeightParaNmae ;//高程模式和高程辅助参数
	CString Leve0ImgName;//待纠正0级影像
	CString Leve1ImgName;//待纠正1级影像
	CString OrbitName, AttName;//待纠正轨道姿态文件
	CString RectImgCatagName;//待纠正影像编目文件
	CString Level0StartScene;//起算景0级影像
	CString Level1StartScene;//起算景1级影像
	CString StartSceneOrbitName, StartSceneAttName;
	
	
	//坐标系统
	CString OutWKT,OutHeightSysID;//输出影像的平面和高程系统
	CString Leve234ImgName,Level234CatagName;//输出影像
	CString BrowserMapName, ThumbMap;//浏览图和拇指图
	/////辐射校正参数
	int	    RefCCD;//参考ccd,
	CString	AbsoParaFile;//均一化辐射校正文件,
	CString	PriorFile;//先验统计文件,
	int     Level1proMode;//一级产品生产模式,
	int     Level1StatMode;//相对辐射校正下的统计模式,
    /////输出的报告文件名    
	CString ProuctReport;
	//
	CString RPCFileName;//add by zhang guo 
	//
	CString DIMFileName;//add by zhang guo

	CString attFile,ephFile,geoFile,imdFile;//add by zhang guo
	double m_ResamRes_x,m_ResamRes_y;

		//add by fei
	CString phasefile;
	CString coherencefile;
	CString interauxfile;
}OrderForm;


struct SimilarMEASURE
{
	int k;
	float xi,yi;
	float coef;
	float h0,h1;
};


typedef struct 
{
	//////从外界得到的信息
	CString mapProjection;
    int zoneNumber;
	CString earthEllipsoid;
	CString resamplingKernel;
    CString elevationCorrection;
	CString baseElevation;
	
	//////从生产中获得的信息
	long productSize;	
	long numLines,numPixels;
	double lineSpacing,pixelSpacing;
	double width,height;
	double centreLocation_latitude,centreLocation_longitude;
	double upperLeft_latitude,upperLeft_longitude;
	double upperRight_latitude,upperRight_longitude;
	double lowerRight_latitude,lowerRight_longitude;
	double lowerLeft_latitude,lowerLeft_longitude;
}METADATASTRUCT;

typedef struct 
{
	double X,Y,Z;
	double Xv,Yv,Zv;
	double UT;
	int year,month,day,hour,minute;
	double second;
	BYTE GPS_STATE;
}CBERS2GPSDATA;

typedef struct 
{  
  double Inclination;
  double Longitude_of_Ascending_Node;
  double Eccentricity;
  double Argument_of_Perigee;
  double Semi_Major_Axis;
  double Ascending_Node_Crossing_Day;
  double Ascending_Node_Crossing_Time;	 
}CBERS2ORBIT;
//on site test
typedef struct 
{
	double x,y;//原始影象坐标
	double east,north;//WGS84UTM0
	double x_,y_;//控制点投影到影象面的坐标	
	double dx,dy,dxy,sitaimage;//绝对误差
	double rdx,rdy,rdxy,rsitaimage;//相对误差	
}GCPSTRUCTEX; //for Absolute Accuracy Test

typedef struct 
{
	GCPSTRUCTEX p1,p2;	
	double A,B;
	double Ax,Bx;
	double Ay,By;
	double Cx,Cy,C;//距离差
	double rCx,rCy,rC;//距离相对误差
	double six,siy,si;//sample interval	
	double Nx,Ny, N;
	double sita,rsita;
}GCPSTRUCTEXEX; //for relative Accuracy Test

typedef struct         
{                      
	double X[512]; 
	double Y[512]; 
	int PointNum;  
	float Angle;   
	int Dis;       
	double a[3];   
                       
}Line; 



struct fPOINT2D	{
	float x,y;
}	; 


typedef struct 
{
	int LineIndex1;
	int LineIndex2;
	double x0,y0;
	double angle;//度
}TwoLine;

typedef struct 
{
	double LNUM[20],LDEN[20],SNUM[20],SDEN[20];
	double LAT_OFF,LAT_SCALE,LONG_OFF,LONG_SCALE,HEIGHT_OFF,HEIGHT_SCALE;
	double LINE_OFF,LINE_SCALE,SAMP_OFF,SAMP_SCALE;
}RPCPARA;


typedef struct VectorPoint
{
	fPOINT2D point;
	VectorPoint *next;
} VectorPoint;


typedef struct VectorLine
{
	int LineNum;
	int LineLength;
	VectorPoint StartPoint;
	VectorPoint EndPoint;
	VectorPoint *MiddlePoint;
	VectorLine  *next;
} VectorLine;

typedef struct
{
	CString ID;
	int mapScale;
	CString WKT;
	double X[4],Y[4];
	double Lat[4],Lon[4];
	double TopLeftX,TopLeftY;
	double I[4],J[4];
}MAPSERIAL;

typedef struct 
{
	double I,J;
	double XY;
	BOOL IsUse;
	
}HV1000POS;


typedef struct 
{
	double I,J;
	double X;
	
}H1000POS;

typedef struct 
{
	double I,J;
	double Y;
	
}V1000POS;

typedef struct 
{
	double I,J;
	double X;
	CString  Xstring;
}H1000POSEX;

typedef struct 
{
	double I,J;
	double Y;
	CString Ystring;
}V1000POSEX;

typedef struct 
{
	double FRAME_LON;
	double FRAME_LAT;
	int    FRAME_ROW;
	int    FRAME_COL;
} Vertex;

typedef struct 
{
	double FRAME_LON;
	double FRAME_LAT;
	int    FRAME_ROW;
	int    FRAME_COL;
} Scene_Center;

typedef struct 
{
	int    p_num;                 //点的序号
	BYTE   property;              //点的属性,0-无效点，1-控制点，2-核查点           
	double CRS_X;                 //控制点X方向
	double CRS_Y;                 //控制点Y方向
	double CRS_Z;                 //控制点Z方向
	double DATA_X;                //图像坐标x
	double DATA_Y;                //图像坐标y
	double mX;                    //x方向的残差
	double mY;                    //y方向的残差
} Tie_Point;
typedef struct 
{
	int nColnum;
	int nRownum;
	int nBandnum;
} RAS_D;
typedef struct 
{
	double Satellite_Altitude;
	double NaDir_Lat;
	double NaDir_Lon;
	CArray <POINT_Ep,POINT_Ep > Ep_Array;
} Ephemeris_Para;

typedef struct st_Quaternion
{
	double x,y,z,w;
}stQuaternion;

typedef struct{
	double LNUM[20],LDEN[20],SNUM[20],SDEN[20];
	double LAT_OFF,LAT_SCALE,LONG_OFF,LONG_SCALE,HEIGHT_OFF,HEIGHT_SCALE;
	double LINE_OFF,LINE_SCALE,SAMP_OFF,SAMP_SCALE;
	double px[3],py[3];
	double invpx[3],invpy[3];
	int nz;
	RPCGCPSTRUCT *rpcgcp;
}RPCMODEL;
#include <Afxtempl.h>
typedef struct DOUBLEPOINT
{
	double x;
	double y;
	
	BOOL operator==(DOUBLEPOINT& dPoint) const
	{
		if(x != dPoint.x || y != dPoint.y)
			return FALSE;
		return TRUE;
		
	}
	BOOL operator!=(DOUBLEPOINT dPoint) const
	{
		return !operator==(dPoint);		
	}
	
}POINT2D;
#define EMPOINT POINT2D 




typedef POINT3D PointStruct;


typedef struct
{
	double x1;
	double y1;
	double x2;
	double y2;
}LINE;

//I0,J0 像主点在扫描坐标系统的坐标
//x,y   象点在像平面坐标系统的坐标
//I,J   像点在扫描坐标系统的坐标
//从像平面坐标系统到扫描坐标系统
//I=I0+(a11*x+a12*y)/pixelsize
//J=J0+(a21*x+a22*y)/pixelsize
//从扫描坐标系统到像平面坐标系统
//x=(b11*(I-I0)+b12*(J-J0))*pixelsize
//y=(b21*(I-I0)+b22*(J-J0))*pixelsize
//pixelsize 象素大小
typedef struct
{
	//I0,J0 像主点在扫描坐标系统的坐标
	double I0,J0;
	//焦距
	double f;
	//从像平面坐标系统到扫描坐标系统
	double a11,a12,a21,a22;
	//从扫描坐标系统到像平面坐标系统
	double b11,b12,b21,b22;
	double pixelsize;
	
}INTERIORORIENTATION;


typedef struct
{
	double Xs;
	double Ys;
	double Zs;
	double fai;
	double omiga;
	double kapa;
}OUTERORIENTATION;


typedef struct  {
	double f;//m
	double pixel0,line0;//m
	double LINE_pixelL,LINE_lineL;//m
	double LINE_angle;	
	double LINE_pixel_size;//m
	double LINE_resample_size;//m
	double lineRange;//m
	double pixelRange;//m	
	int LINE_DYNAMIC_SAMPLE_NUM;
	int LINE_pixel_num;
	int Line_line_num;
	double LINE_PERIOD;

	// these tow are used in camera simulation
	double pixel_size;            // add by pan, it should be LINE_pixel_size
	double pixel_num;             // add by pan, it should be LINE_pixel_num
}CAMERA;
// add by zhangguo 20081004


typedef struct {
	
	double f;//m
	double sample0,line0;//mm
	double LINE_sampleL,LINE_lineL;//m
	double LINE_angle;//无效
    double lineRange;
	double pixel_size;

	double anglesampleerror;
	double anglelineerror;
	
	
	
}CAMERAEx;
// add by zhangguo 20101013

typedef struct
{
	int     pyramid_nlayer;
	CSize   image_size[5];	
	LPBYTE  pyramid_image[5];
}PYRAMID;

typedef struct
{
	double x,y;
	double coef;			
}cof2POINT;

typedef struct 	
{
	int trip;
	int order;
	CString  filename;	
	CString  cmrname;
}IMAGEFILE;

typedef struct  
{
	double dx,dy,dz;
	double fai,omiga,kapa;
	double lanmda;
}SIMILARITYTRANSFORMATION;

typedef struct
{
	int nImgID;
	POINT2D offsetPoint;			
}OFFSETPOINT;
typedef struct
{
	int nPointID;
	OFFSETPOINT topleftPoint;
	OFFSETPOINT toprightPoint;
	OFFSETPOINT bottomleftPoint;
	OFFSETPOINT bottomrightPoint;
}ARRAY_OFFSETPOINT;

typedef struct 
{
	int left;
	int right;
}LR_tri;

typedef struct 
{
	int p_index[3];
	int tri_index[3];
}Triangle;

typedef struct tagTiePoint
{
	CString id;	
	CArray <int,int> overlay;
	CArray <POINT2D,POINT2D> point;
	tagTiePoint(tagTiePoint& tag)
	{
		id=tag.id;
		overlay.Copy(tag.overlay);
		point.Copy(tag.point);
	};
	tagTiePoint()
	{
		id="";
		overlay.RemoveAll();
		point.RemoveAll();
	};
	void operator=(tagTiePoint& tag) 
	{
		id="";
		overlay.RemoveAll();
		point.RemoveAll();
		id=tag.id;
		overlay.Copy(tag.overlay);
		point.Copy(tag.point);		
	}
}TIEPOINTEX;


typedef struct{
	double x;
	double y;	
	double dx,dy;//R-L	
} MATCHSTRUCT;

typedef struct  {
	double x,y;
	double w;
	__int64 id;
}POINT2DINDEX;



typedef struct 
{
	__int64 id;
	double east,north,h;
	double lat,lon;
	int state;//0 检查点 1 控制点 2 待定点 3 结束符号
}GCPINDEX;

typedef struct 
{
	__int64 id;
	int imgnum;
	double *x,*y;
	double *w;
	int *imgid;
}TIEPOINT;


typedef struct 
{
	double EPH_X,EPH_Y,EPH_Z;
	double STAR_ROLL,STAR_PITCH,STAR_YAW;
	double CAMERA_X,CAMERA_Y,CAMERA_Z;
	double CAMERA_ROLL,CAMERA_PITCH,CAMERA_YAW;
}EPH_STAR_CAMERA_CONFIG;


//add by zhangguo 20081002
typedef struct  
{
	double dx,dy,dz; 
	double sita_z,sita_y,sita_x; //绕XYZ轴旋转
}SimilarityTransformation;


typedef struct 
{
	double east_lon;
	double north_lat;
	double h;
	CString  east_lonstr;
	CString  north_latstr;
	CString  hstr;

	double lon;
	double lat;
	double H;


	CString wkt;
	CString datum;
	CString project;
	CString information;
}GPSCONTROLPOINT;

struct matchPOINT	{
	float xl, yl;
	float xr, yr;
};

typedef struct 	
{
	float xl, yl;
	float xr, yr;
	float c;
}matchPOINT_ZG;




typedef struct
{
	float x,y;
	float coef;			
}cof2POINTfloat;

typedef struct 	
{
	float xl, yl;
	cof2POINTfloat xyr[12];
	int num;
	int Isdel;
}matchPOINT1;

// 控制点匹配候选
struct mtCANDIDATE1
{	
	float  xi, yi;	// 影像坐标
	float coef;		// 相关系数
	float pb;		// 匹配概率
	float vx, vy;	// 残差
};

// 松弛节点
struct mtRelaxNODE	
{
	int	numOfCnds;			// 候选点数
	mtCANDIDATE1  mtCnds[5];	// 候选点
	float p0;		// 零匹配概率
};

typedef struct 
{
	double Xs[21];
	double Ys[21];
	double Zs[21];
	double Xvs[21];
	double Yvs[21];
	double Zvs[21];
	double UT_OFF,UT_SCALE;
	double Xs_OFF,Xs_SCALE;
	double Ys_OFF,Ys_SCALE;
	double Zs_OFF,Zs_SCALE;
	double Xvs_OFF,Xvs_SCALE;
	double Yvs_OFF,Yvs_SCALE;
	double Zvs_OFF,Zvs_SCALE;
} OrbitModel;


typedef struct  
{
	int IsBIASANGLE; //1 校正 0 不校正
	double BIASANGLEControlError;//度
    double YAW0;//姿态初始数值
	double PITCH0;//姿态初始数值
	double ROLL0;//姿态初始数值
	int sidelookControlNum;//侧摆控制次数
    double *sidelookControlStartTime;//侧摆控制起始时间
	double *sidelookControlTimeRange;//侧摆控制持续周期
	double *dYAW0;//侧摆控制速度
	double *dPITCH0;//侧摆控制速度
	double *dROLL0;//侧摆控制速度
	double MaxYAW;//最大角度控制
	double MaxPITCH;//最大角度控制
	double MaxROLL;//最大角度控制
	double YAWSta;//姿态稳定度
	double PITCHSta;//姿态稳定度
	double ROLLSta;//姿态稳定度
	double YAWError;//姿态测量误差
	double PITCHError;//姿态测量误差
	double ROLLError;//姿态测量误差
	double ATTTIMESpan;//姿态测量周期
	double ATTTIMEError;//姿态量测延时误差, add by pan
}AttitudeInput;



typedef struct  
{
	int IsBIASANGLE; //1 校正 0 不校正
	double BIASANGLEControlError;//度

	CString OscillationFileName;

	Attitude att0; //成像时候姿态 有效为 ROLL PITCH YAW

	int sidelookControlNum;//侧摆控制次数
    double *sidelookControlStartTime;//侧摆控制起始时间
	double *sidelookControlTimeRange;//侧摆控制持续周期
	double *dYAW0;//侧摆控制速度
	double *dPITCH0;//侧摆控制速度
	double *dROLL0;//侧摆控制速度
	double MaxYAW;//最大角度控制
	double MaxPITCH;//最大角度控制
	double MaxROLL;//最大角度控制



	Attitude attSta;//姿态稳定度 有效为 ROLL PITCH YAW

	int AttMode;//0 星敏陀螺联合处理  1 陀螺单独处理  2  星敏单独处理
    
	Attitude attError;//有效为 ROLL PITCH YAW UT 代表时间同步误差

	//针对姿态模式2 
	Attitude STARInstall;//  有效为 ROLL PITCH YAW 星敏安装矩阵
	Attitude STARInstallError;//  有效为 ROLL PITCH YAW 星敏安装矩阵误差
	CString  STARInstallcyclefile;//星敏安装周期文件，为目录模式，

	double ATTTIMESpan;//姿态测量周期

	//相机安装矩阵部分

	Attitude cameraInstall;// roll pitch yaw 有效
	Attitude cameraInstallError;// roll pitch yaw 有效	
	CString camerainstallcyclefilename;//相机安装周期变化文件

}AttitudeInputEx;






// added by li yang
typedef struct 
{
	BOOL change;//true 变化 false 没有变化
	POINT *boundaryPointList;// x 行方向 y列方向
	long boundaryPointNum;
	POINT topleft,bottomright;//区域范围 
	long w,h;
	long bandnum;
	float *bandinfor;
	long MPC;//区域内部点个数
	
	double *distance;
	
	double cof;
	
}OBJECT;

// added by li yang
typedef struct 
{
	POINT *PointList;
	POINT startP;
	POINT endP;
	int ptNum;
	BOOL bUsed;
}SECTION;


// add by pan
enum ErrorTypeForSIM
{
	// the first number is the error type, and the last number is the modal
	// 总共4个模块，分别表示轨道、姿态、CCD和成像
	// 误差共7种类型，轨道观测、轨道同步、姿态误差、姿态同步、安装误差、CCD倾斜和成像同步误差
	OrbitPosErr = 101,              // 轨道位置误差
	OrbitVelErr = 111,              // 轨道速度误差
	OrbitTimeErr = 201,             // 轨道同步误差
	AttCombineErr = 302,            // 姿态联合误差
	AttGyroErr = 312,               // 陀螺观测误差
	StarSObErr = 322,               // 星敏观测误差
	StarSInstErr = 332,             // 星敏安装误差
	AttTimeErr = 402,               // 姿态同步误差
//	InstPosErr = 502,               // 安装位置误差
	InstAttErr = 512,               // 安装姿态误差
	CCDSlantErr = 603,              // CCD倾斜误差
	ImageTimeErr = 704              // 成像同步误差
};

// add by pan
typedef struct  
{
	CString   ErrorName;     // ?
	ErrorTypeForSIM eErrType;
	double    ErrBegin;
	double    ErrEnd;
	double    ErrStep;

	BOOL      flag;         // 1使用，0无效;
	int       Num;          // 误差个数
} SimulatorError;

typedef struct					//ADD BY JYH
{
	double x;
	double y;
	double lat;
	double lon;
	double ex;
	double ey;
	double exy;
}ErrAnalyData;


typedef signed char         INT8;
typedef signed short        INT16;
typedef signed int          INT32;
typedef signed __int64      INT64;
typedef unsigned char       UINT8;
typedef unsigned short      UINT16;
typedef unsigned int        UINT32;
typedef unsigned __int64    UINT64;
//增加一些定义
#endif // !defined(AFX_GEORSDATASTRUCT_H__261DDC12_AA2A_47F6_97AD_DAF7D223488C__INCLUDED_)
