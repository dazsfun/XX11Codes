// ConsoleApplication1.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include "safestring.h"
#include "VirtualMachine.h"
#include "jsonSerializer.h"
#include "AttitudeAPI.h"
#include "OrbitAPI.h"
#include "LineTimeAPI.h"
#include <string>
#include <fstream>
#include "ImgFacotryAPI.h"
#include <gdal_priv.h>
#include <windows.h>
#include "RandomObject.h"
using namespace std;


#define   _MC    1     // MLCmalloc()�ı�־
#define   _FE    0     // MLCfree()�ı�־
#define  _Check  -1    //����ڴ�ı�־




typedef struct _Memory         //��������ṹ
{
	unsigned int address;
	struct _Memory *next;
}Mem, *LinkMem;

int MemoryContrl(size_t addr, char flag)
{
	static Mem head;        //��������ͷ     
	LinkMem tmp = head.next, node = &head;
	if (_MC == flag)         //��������ڵ�
	{
		while (tmp)
		{
			node = tmp;
			tmp = tmp->next;
		}
		node->next = (LinkMem)malloc(sizeof(Mem));
		if (node->next == NULL)
			return 1;   //����ʧ��!
		node = node->next;
		node->address = addr;
		node->next = NULL;
	}
	else if (_FE == flag)     //ɾ���ڵ�
	{
		while (tmp && (tmp->address != addr))
		{
			node = tmp;
			tmp = tmp->next;
		}
		if (tmp == NULL)
			return -1;  //free��ַ����!
		else
		{
			node->next = tmp->next;
			free(tmp);
		}
	}
	else                           //��й©���ص�ַ�����򷵻�null
		return (int)(head.next ? head.next->address : 0);
	return 0;
}


///////////////�ڴ����//////////////////////


void *MLCmalloc(size_t size)
{
	void *pMalc = malloc(size);
	if (pMalc == NULL)
		return NULL;
	else if (MemoryContrl((size_t)pMalc, _MC) == 1)
	{
		free(pMalc);   //��̬������ʧ��ʱ�ͷŴ˴η���Ŀռ�
		return NULL;
	}
	else
		return pMalc;
}

///////////////�ڴ��ͷ�//////////////////////

void MLCfree(void *ptr)
{
	if (MemoryContrl((size_t)ptr, _FE) == -1)
	{
		//printf("%s/n", "���ͷŵ��ڴ��ַ����!");
		return;
	}
	else
		free(ptr);
}


struct Memery_Node {
	//char *pcFileName;
	unsigned int ulSize;
	//unsigned int ulLine;
	void *memAddr;
};
struct Memry_List {
	struct Memery_Node *Mem_Node;
	struct Memry_List *Next;
};

struct Memry_List gHead_Node = { NULL, NULL };
#define PID_TEST_MLC 100
#define DYNAMIC_MEM_PT 1
void *MLC_MemAlloc( unsigned int ulSize)
{
	struct Memery_Node *Mem_Node = NULL;
	struct Memry_List *List_Node = NULL;
	struct Memry_List *Curr_Node = NULL;
	struct Memry_List *Pre_node = NULL;
	void *pAddr = NULL;
	unsigned int ulPidSid = PID_TEST_MLC;
	unsigned char ucPtNo = DYNAMIC_MEM_PT;


	Mem_Node = (struct Memery_Node *)malloc(sizeof(struct Memery_Node));
	if (NULL == Mem_Node)
		return NULL;
	memset(Mem_Node, 0, sizeof(struct Memery_Node));

	pAddr = (void *)malloc(ulSize*sizeof(char));
	if (NULL == pAddr) {
		free(Mem_Node);
		Mem_Node = NULL;
		return NULL;
	}
	memset(pAddr, 0, ulSize);

	Curr_Node = (struct Memry_List *)malloc(sizeof(struct Memry_List));
	if (NULL == Curr_Node) {
		cout << "nonono" << endl;
		free(pAddr);
		pAddr = NULL;
		free(Mem_Node);
		Mem_Node = NULL;
		return NULL;
	}
	memset(Curr_Node, 0, sizeof(struct Memry_List));

	Mem_Node->ulSize = ulSize;
	//Mem_Node->pcFileName = pcFileName;
	//Mem_Node->ulLine = ulLine;
	Mem_Node->memAddr = pAddr;

	Curr_Node->Mem_Node = Mem_Node;
	Curr_Node->Next = NULL;

	List_Node = &gHead_Node;
	while (NULL != List_Node) {
		Pre_node = List_Node;
		List_Node = List_Node->Next;
	}

	Pre_node->Next = Curr_Node;
	return pAddr;

}

void MLC_MemFree( void *pAddr)
{
	/*��������ɹ���*/
	struct Memery_Node *Mem_Node = NULL;
	struct Memry_List *List_Node = NULL;
	struct Memry_List *Curr_Node = NULL;
	struct Memry_List *Pre_node = NULL;

	if (NULL == pAddr || (NULL == (gHead_Node.Next)->Mem_Node) || gHead_Node.Next == NULL)
		return;

	Curr_Node = &gHead_Node;
	while (NULL != Curr_Node)
	{
		Pre_node = Curr_Node;
		Curr_Node = Curr_Node->Next;
		if (NULL == Curr_Node) {
			return;
		}
		Mem_Node = Curr_Node->Mem_Node;
		if (Mem_Node->memAddr == pAddr)
			break;
	}
	Pre_node->Next = Curr_Node->Next;

	free(Curr_Node->Mem_Node);
	Curr_Node->Mem_Node = NULL;
	free(Mem_Node->memAddr);
	Mem_Node->memAddr = NULL;
	free(Curr_Node);
	Curr_Node = NULL;
	return;
}
//////////���й©//////////
unsigned int MLC_CheckIfMemLeak(unsigned int *pulAddrNum, void* pAddrArray[])
{
	/*��������ɹ���*/
	struct Memery_Node *Mem_Node = NULL;
	struct Memry_List *List_Node = NULL;
	struct Memry_List *Curr_Node = NULL;
	struct Memry_List *Pre_node = NULL;
	unsigned int ulAddrNum = 0;



	Curr_Node = gHead_Node.Next;
	while (NULL != Curr_Node) {
		pAddrArray[ulAddrNum++] = Curr_Node->Mem_Node->memAddr;
		Curr_Node = Curr_Node->Next;
	}

	if (ulAddrNum > 0)
		return -1;

	return 0;
}


void MLC_CollectMemory(void)
{
	/*��������ɹ���*/
	struct Memery_Node *Mem_Node = NULL;
	struct Memry_List *List_Node = NULL;
	struct Memry_List *Curr_Node = NULL;
	struct Memry_List *Pre_node = NULL;
	if (NULL == gHead_Node.Next || NULL == (gHead_Node.Next)->Mem_Node)
		return;

	Curr_Node = gHead_Node.Next;
	while (NULL != Curr_Node) {
		memcpy(Pre_node, Curr_Node, sizeof(struct Memry_List));
		free(Curr_Node->Mem_Node->memAddr);
		Curr_Node->Mem_Node->memAddr = NULL;
		free(Curr_Node->Mem_Node);
		Curr_Node->Mem_Node = NULL;
		free(Curr_Node);
		Curr_Node = Pre_node->Next;
	}
	gHead_Node.Next = NULL;
}
size_t CheckLeak()
{
	return (size_t)(MemoryContrl((size_t)0, _Check));
}

#ifndef nullptr
#define nullptr 0

#include <memory.h>
void* g_paMemLeakAddrArray[2048] = { 0 };
int main()
{
	char woringmode[1024];
	cout << "��ѡ����ģʽ" << endl;
	cout << "ȫ���̲���:" << "SysTest000" << endl;
	cout << "��ͶӰ��Ԫ����:" << "UniTest001" << endl;
	cout << "����У����Ԫ����" << "UniTest002" << endl;
	cin >> woringmode;
	if (safestring::compare(woringmode, "UniTest001", 10))
	{
		char unitestconfigfilepath[1024];
		cout << "��������ͶӰ��Ԫ���������ļ�·��:" << endl;
		//cin >> unitestconfigfilepath;
		sprintf_s(unitestconfigfilepath, "%s", "D:\\2015-3-31-14-58-56\\imagemodel");
		VirtualMachine vir("UnitTest", false);
		double testx, testy;
		testx = 5070; testy = 240;
		double height = 100;
		double lontest, latest;
		vir.FromXY2LonLatUnitTest(testx, testy, height, lontest, latest, 0, static_cast<string>(unitestconfigfilepath));
	}

	else if (safestring::compare(woringmode, "SysTest000", 10))
	{
		string hei = "D:\\IRS\\JB11-1_IRS_000116063_003_002_L0\\JB11-1_IRS_000116063_003.FRED";
		VirtualMachine envAPI("C93INFRAED201502170516PMsDateStr.log", hei);
	}

	else if (safestring::compare(woringmode, "UniTest002", 10))
	{
		RandomObject robject(480, 10786);

		string hei = "Whole.tif";
		ImgFacotryAPI imgapitemp(480, 10786, 12, "hei");
		GDALAllRegister();

		int currentImageID = 62;
		imgapitemp.Begin_Match_WorkFlow("C:\\Temp\\match\\left_.tif", "c:\\temp\\match\\right_", nullptr);
		//imgapitemp.Begin_Match_WorkFlow("D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\Whole.tif21_.tif",
		//	"D:\\IRS\\JB11-1_IRS_000116063_003_002_L0\\JB11-1_IRS_000116063_003_004_002_L3C_TOABT.tiff",
		//	nullptr);
		//imgapitemp.Mosaic(imagenames, imgframewidth * imagenames.size(), imgframeheight, false);
		//imgapitemp.DetectorHistoMatchCalibration<unsigned short>(hei, 480, 10876, 12, GDT_UInt16);
		/*
		vector<string> heiyou;
		heiyou.push_back("c:\\temp\\testRC\\Whole.tif0_.tif");
		heiyou.push_back("c:\\temp\\testRC\\Whole.tif1_.tif");
		heiyou.push_back("c:\\temp\\testRC\\Whole.tif2_.tif");
		imgapitemp.StitchFile(heiyou, "whatever.tif", 480, 10876);
		*/
		//imgapitemp.ImageScale("D:\\OpenCVWithCuda\\RPCModel\\ConsoleApplication1\\ConsoleApplication1\\whateverits.tif", 0, 10500, 12080, 11905, 9.9 / 6.5, 1);
		
	}
	else
	{
		cout << "����ģ��δ��ʶ��" << endl;
		cout << "��������˳�" << endl;
		int good;
		cin >> good;
		return 0;
	}

}
#endif