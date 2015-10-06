#ifndef RCPROCESS_H
#define RCPROCESS_H
/*---------------------------------------------------------------------------------------------------------/
*函数名称：bool X11RCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[])
*参数：const char *str_FileName  //影像文件绝对路径名称
*      const char *str_RCType    //处理方式，可选LibCoe,Stat，Histo，MeanStat 其中Stat 和 MeanStat 速度较快，Histo 很慢
*                                //对X11建议Stat && MeanStat
*      int OverLapPixelNum[]     //各片CCD重叠像素个数，这里处理X11，默认8片CCD，需要提供7个重叠区像素个数
*返回值：true -- success, false -- failed
*功能说明：对X11数据进行初步辐射处理，其余类似
*---------------------------------------------------------------------------------------------------------*/
bool __declspec(dllexport) X11RCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[]);
bool __declspec(dllexport) X10RCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[]);
bool __declspec(dllexport) O2CRCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[]);
bool __declspec(dllexport) X6RCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[]);
bool __declspec(dllexport) SJ9RCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[]);
#endif
