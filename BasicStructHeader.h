#ifndef _BASIC_ELEMENT_TYPEDEFSTRUCT_H_
#define  _BASIC_ELEMENT_TYPEDEFSTRUCT_H_
#include <iostream>
#include <vector>
using namespace std;

#ifndef PI
#define  PI 3.1415926535897932384626433832795
#endif

typedef struct ImgP
{
       double x;
       double y;
}ImgS,*pImgS;

typedef struct CheckResult
{
	double _maxerror;
	double _minerror;
	double _rmserror;
}ChR,*pCheckResult;

enum ErrorType_1
{
	No_Error,
	Img_ItsNotAnImage,
	Img_NoBaseImage,
	Img_CreateNewImageFailed,
	Cache_AllocateCacheFailed,
	Misclosure_ProjectionTrack_LEFT_RIGHT_LEFT,
	Fitting_ExceedLimation_LST_EPIPOLARLINE,
	Web_API_BingMAP_FAILED
};

enum ClassName
{
	EpiMakerClass
};

static string ClassNames[1] = { "EpiImg" };

typedef struct GroundP
{
         double  X;
         double  Y;
         double  Z;

		 double lat;
		 double lon;
		 double h;
}GP,*pGp;

typedef struct FPoint    //Control Point
{
      char id[2];
      ImgP Pixel;
      GroundP ground;

	  friend ostream &operator<<(ostream &os,const FPoint &c)
	  {
		  os<<"Point:"<<c.id<<endl;
		  os<<"Point type:"<<"Control point"<<endl; 
		  os<<"Image coordinates:"<<endl;
		  os<<c.Pixel.x<<"	"<<c.Pixel.y<<endl;

		  os<<"latitude and longitude coordinates:"<<endl;
		  os<<c.ground.lat<<"	"<<c.ground.lon<<"	"<<c.ground.h<<endl;

		  os<<"Ground coordinates:"<<endl;
		  os<<c.ground.X<<"	"<<c.ground.Y<<"	"<<c.ground.Z<<endl;//Adding
		 
		  os<<"	"<<endl;

		  return os;
	  }
};

typedef struct SPoint  //Check Point
{
      char id[256];
      ImgP Pixel;
      GroundP ground;

	  friend ostream &operator<<(ostream &os,const SPoint &c)///by reload operator "<<",we can print the object in way we feel like to
	  {
		  os<<"Point:"<<c.id<<endl; //Adding, 
		  os<<"Point type:"<<"Check Point"<<endl;
		  os<<"Image coordinates:"<<endl;
		  os<<c.Pixel.x<<"	"<<c.Pixel.y<<endl;
		  os<<"Ground coordinates:"<<endl;
		  os<<c.ground.X<<"	"<<c.ground.Y<<"	"<<c.ground.Z<<endl;
		  os<<"	"<<endl;

		  return os;
	  }
};
#endif