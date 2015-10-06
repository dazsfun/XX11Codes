#pragma once
#include "AuxFactory.h"
#include "GeoBase.h"
#include "jsonSerializer.h"
#include "Safeauxilary.h"
/*
让所有主要类继承jsonSerializer类的原因：
让每次成像模型数据，信息，建立过程及结果都是可记录，可追溯，可比较的
并尽可能减小这种操作的开销
*/
class ImageModel : public JSONOBJECTSerializer
{
public:
	///这个默认构造函数将被取消
	CSatOrbit helper;
	ImageModel(OrbitAPI* orbitComp,AttitudeAPI* attComp,LineTimeAPI* linetimeapi);
	ImageModel(ImageModel* p_model) :currentOrbit(OrbitAPI("currentOrbit")), currentAttitude(AttitudeAPI("currentAttitude"))
	{
		*this = *p_model;
	}
	OrbitAPI currentOrbit;
	AttitudeAPI currentAttitude;
public:
	//不允许拷贝
	ImageModel(const ImageModel& _model):currentOrbit(OrbitAPI("currentOrbit")), currentAttitude(AttitudeAPI("currentAttitude"))
	{
		*this = _model;
	}
	~ImageModel();

public:
	virtual void SetPropertys() {};
	virtual void SetUTCAPIByTimestring() {};


	////牢记的拷贝，复制，以及拷贝构造
	/*
	std::vector<type> source;  // 假定source中已经保存了您的数据
std::vcctor<type> dst;     // 你想把source中的数据拷贝到 dst中

// 初始化的dst没有分配内存， 这里是分配等同大小的内存 你也可以根据需要分配
dst.resize(source.size());

// 记住 vector中分配的是一个块内存， 然后再次基础上作动态数据定义
// 成员函数获取的是数组的个数， 所以得乘以 sizeof(type)
memcpy(&dst[0], &source[0], source.size() * sizeof(type)); //搞定*/ ///注:此处并不可使用这种做法拷贝

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
	//获取速度
	double VelX(double time)
	{
		OrbitAPI currentOrbit("Point");
		m_base.LagrangianInterpolation(_orbitComponents, time, _orbitComponents.size(), currentOrbit);
		return currentOrbit.GPSXV;
	}

	///获取速度矢量方向
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
				itr++;//这里迭代器会更新出错
		} while (itr != _content.end());
	}
};

