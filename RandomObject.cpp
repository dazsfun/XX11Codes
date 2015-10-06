#include "stdafx.h"
#include "RandomObject.h"





RandomObject::~RandomObject()
{
	delete[] _rightx;
	delete[] _righty;
}
