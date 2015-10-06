#include "stdafx.h"
#include "Safeauxilary.h"

#include "log.h"

double Safeauxilary::InterPolation(const double& seed, const vector<double>& _seeds)
{
	int seedend1 = static_cast<int>(seed);
	if (seedend1 <0 || seedend1 > _seeds.size() - 1)
	{
		cout << "超出插值有效范围" << endl;
		LogAPI logs;
		cout << seedend1 << endl;
		logs.IAddLog("超出插值有效范围");
		return -999999.0;
	}

	if (seedend1 == _seeds.size() - 1)
		return _seeds[_seeds.size() - 1];

	double edge = seed - seedend1;
	return _seeds[seedend1] * (1 - edge) + _seeds[seedend1 + 1] * edge;
}



bool Safeauxilary::GetJ2000ToWGS84(int year, int month, int day, int hour, int  minute, double second, double* rotMatrix, double* rotMatrixDot)
{
	int jd, i, j;
	double sod;
	char file_table[256] = "file_table";
	DIR_TABLES(file_table, 256);



	MJD_SOD(&year, &month, &day, &hour, &minute, &second, &jd, &sod);

	double rotmat[3][3] = { 0 };
	double matrixdot[3][3] = { 0 };
	ECIF2ECEF(&jd, &sod, &rotmat[0][0], &matrixdot[0][0]);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			rotMatrix[j * 3 + i] = rotmat[i][j];
			rotMatrixDot[j * 3 + i] = matrixdot[i][j];
		}
	}
	return true;
}