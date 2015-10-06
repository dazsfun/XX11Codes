#ifndef RCPROCESS_H
#define RCPROCESS_H
/*---------------------------------------------------------------------------------------------------------/
*�������ƣ�bool X11RCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[])
*������const char *str_FileName  //Ӱ���ļ�����·������
*      const char *str_RCType    //����ʽ����ѡLibCoe,Stat��Histo��MeanStat ����Stat �� MeanStat �ٶȽϿ죬Histo ����
*                                //��X11����Stat && MeanStat
*      int OverLapPixelNum[]     //��ƬCCD�ص����ظ��������ﴦ��X11��Ĭ��8ƬCCD����Ҫ�ṩ7���ص������ظ���
*����ֵ��true -- success, false -- failed
*����˵������X11���ݽ��г������䴦����������
*---------------------------------------------------------------------------------------------------------*/
bool __declspec(dllexport) X11RCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[]);
bool __declspec(dllexport) X10RCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[]);
bool __declspec(dllexport) O2CRCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[]);
bool __declspec(dllexport) X6RCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[]);
bool __declspec(dllexport) SJ9RCProcessed(const char *str_FileName,const char *str_RCType,int OverLapPixelNum[]);
#endif
