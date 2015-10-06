#pragma once
#include "AuxFactory.h"
#include "GeoBase.h"
#include "jsonSerializer.h"
#include "Safeauxilary.h"
/*
��������Ҫ��̳�jsonSerializer���ԭ��
��ÿ�γ���ģ�����ݣ���Ϣ���������̼�������ǿɼ�¼����׷�ݣ��ɱȽϵ�
�������ܼ�С���ֲ����Ŀ���
*/
class ImageModel : public JSONOBJECTSerializer
{
public:
	///���Ĭ�Ϲ��캯������ȡ��
	CSatOrbit helper;
	ImageModel(OrbitAPI* orbitComp,AttitudeAPI* attComp,LineTimeAPI* linetimeapi);
	ImageModel(ImageModel* p_model) :currentOrbit(OrbitAPI("currentOrbit")), currentAttitude(AttitudeAPI("currentAttitude"))
	{
		*this = *p_model;
	}
	OrbitAPI currentOrbit;
	AttitudeAPI currentAttitude;
public:
	//��������
	ImageModel(const ImageModel& _model):currentOrbit(OrbitAPI("currentOrbit")), currentAttitude(AttitudeAPI("currentAttitude"))
	{
		*this = _model;
	}
	~ImageModel();

public:
	virtual void SetPropertys() {};
	virtual void SetUTCAPIByTimestring() {};


	////�μǵĿ��������ƣ��Լ���������
	/*
	std::vector<type> source;  // �ٶ�source���Ѿ���������������
std::vcctor<type> dst;     // �����source�е����ݿ����� dst��

// ��ʼ����dstû�з����ڴ棬 �����Ƿ����ͬ��С���ڴ� ��Ҳ���Ը�����Ҫ����
dst.resize(source.size());

// ��ס vector�з������һ�����ڴ棬 Ȼ���ٴλ���������̬���ݶ���
// ��Ա������ȡ��������ĸ����� ���Եó��� sizeof(type)
memcpy(&dst[0], &source[0], source.size() * sizeof(type)); //�㶨*/ ///ע:�˴�������ʹ��������������

	void operator=(const ImageModel& _info)
	{
		//_orbitComponents.resize(_info._orbitComponents.size());
		//_attComponents.resize(_info._attComponents.size());
		//_timeComponents.resize(_info._timeComponents.size());
		
		for (int i = 0; i < _info._orbitComponents.size(); i++)
		{
			_orbitComponents.push_back(_info._orbitComponents[i]);
		}

		for (int i = 0; i < _info._attComponents.size(); i++)
		{
			_attComponents.push_back(_info._attComponents[i]);
		}

		for (int i = 0; i < _info._timeComponents.size(); i++)
		{
			_timeComponents.push_back(_info._timeComponents[i]);
		}

	}

private:
	vector<OrbitAPI> _orbitComponents;
	vector<AttitudeAPI> _attComponents;
	vector<LineTimeAPI> _timeComponents;
	/*
	OrbitAPI* _orbitComponents;
	AttitudeAPI* _attComponents;
	LineTimeAPI* _timeComponents;
	*/
	int attType;
	GeoBase m_base;

public:
	//��ȡ�ٶ�
	double VelX(double time)
	{
		OrbitAPI currentOrbit("Point");
		m_base.LagrangianInterpolation(_orbitComponents, time, _orbitComponents.size(), currentOrbit);
		return currentOrbit.GPSXV;
	}

	///��ȡ�ٶ�ʸ������
	double VelAngle(double time)
	{
		OrbitAPI currentOrbit("currentPoint");
		m_base.LagrangianInterpolation(_orbitComponents, time, _orbitComponents.size(), currentOrbit);

		return atan(currentOrbit.GPSYV / currentOrbit.GPSXV);
	}

	bool SelfCheck();
	void IBody2WGS84(double time, double * body2WGS84Rot,double* position);
	void IBody2WGS84(string _configfile, double* body2WGS84, double* positionXYZ);
	void IBody2WGS84Test(double time, double* body2WGS84Rot, double* position);
protected:
	template<typename T> void ClearVectorSame(vector<T>& _content)
	{
		vector<T>::iterator itr = _content.begin();
		do
		{
			if (*itr == *(itr + 1))
			{
				itr = _content.erase(itr);
			}
			else
				itr++;//�������������³���
		} while (itr != _content.end());
	}
};

