// ConsoleApplication1.cpp : Defines the entry point for the console application.
//
#include <intrin.h>
#include "stdafx.h"
#include <iostream>
#include <omp.h>

using namespace std;

#pragma intrinsic(__rdtsc);

unsigned __int64 t0_init, t_init;
const int N = 2000;
//================================================
int A_int[N][N];
int B_int[N][N];
int C_int[N][N];

float A_float[N][N];
float B_float[N][N];
float C_float[N][N];

double A_double[N][N];
double B_double[N][N];
double C_double[N][N];
//================================================

void init()
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A_int[i][j] = 1;
			B_int[i][j] = 1;
			C_int[i][j] = 0;

			A_float[i][j] = 1.0;
			B_float[i][j] = 1.0;
			C_float[i][j] = 0.0;

			A_double[i][j] = 1.0;
			B_double[i][j] = 1.0;
			C_double[i][j] = 0.0;
		}
	}
}

double time_count(int N, char type, int parallelism_enabled = 0, double fr = 4000000000)
{
	omp_set_num_threads(6);
	//extern int parallelism_enabled;
	unsigned __int64 t0 = __rdtsc();
	switch (type)
	{
	case 'i':
		#pragma omp parallel for if(parallelism_enabled)
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					C_int[i][j] += A_int[i][k] * B_int[k][j];
				}
			}
		}
		break;
	case 'f':
		#pragma omp parallel for if(parallelism_enabled)
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					C_float[i][j] += A_float[i][k] * B_float[k][j];
				}
			}
		}
		break;
	case 'd':
		#pragma omp parallel for if(parallelism_enabled)
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					C_double[i][j] += A_double[i][k] * B_double[k][j];
				}
			}
		}
		break;
	}
	unsigned __int64 t = __rdtsc();
	return (t-t0) / fr;
}

int _tmain(int argc, _TCHAR* argv[])
{
	t0_init = __rdtsc();
	init();
	t_init = __rdtsc();
	cout << "time for init: " << (t_init-t0_init)/ 4000000000 << endl;

	//Outing table
	cout << "posled" << endl;
	cout << "type\t10\t100\t200\t400\t800\t1000\t2000" << endl;
	cout << "--------------------------------------------------------------" << endl;
	char type = 'i';
	cout << "i\t" << time_count(10,type) << "\t" << time_count(100, type) << "\t" << time_count(200, type) << "\t" << time_count(400, type) << "\t" << time_count(800, type) << "\t" << time_count(1000, type) << "\t" << time_count(2000, type) << endl;
	//-------------------------
	type = 'f';
	cout << "f\t" << time_count(10, type) << "\t" << time_count(100, type) << "\t" << time_count(200, type) << "\t" << time_count(400, type) << "\t" << time_count(800, type) << "\t" << time_count(1000, type) << "\t" << time_count(2000, type) << endl;
	//-------------------------
	type = 'd';
	cout << "d\t" << time_count(10, type) << "\t" << time_count(100, type) << "\t" << time_count(200, type) << "\t" << time_count(400, type) << "\t" << time_count(800, type) << "\t" << time_count(1000, type) << "\t" << time_count(2000, type) << endl;
	//-------------------------

	cout << "parallel" << endl;
	cout << "type\t10\t100\t200\t400\t800\t1000\t2000" << endl;
	cout << "--------------------------------------------------------------" << endl;
	type = 'i';
	cout << "i\t" << time_count(10, type, 1) << "\t" << time_count(100, type, 1) << "\t" << time_count(200, type, 1) << "\t" << time_count(400, type, 1) << "\t" << time_count(800, type, 1) << "\t" << time_count(1000, type, 1) << "\t" << time_count(2000, type, 1) << endl;
	//-------------------------
	type = 'f';
	cout << "f\t" << time_count(10, type, 1) << "\t" << time_count(100, type, 1) << "\t" << time_count(200, type, 1) << "\t" << time_count(400, type, 1) << "\t" << time_count(800, type, 1) << "\t" << time_count(1000, type, 1) << "\t" << time_count(2000, type, 1) << endl;
	//-------------------------
	type = 'd';
	cout << "d\t" << time_count(10, type, 1) << "\t" << time_count(100, type, 1) << "\t" << time_count(200, type, 1) << "\t" << time_count(400, type, 1) << "\t" << time_count(800, type, 1) << "\t" << time_count(1000, type, 1) << "\t" << time_count(2000, type, 1) << endl;
	//-------------------------
	
	cout << "accsel" << endl;
	cout << "type\t10\t100\t200\t400\t800\t1000\t2000" << endl;
	cout << "--------------------------------------------------------------" << endl;
	type = 'i';
	cout << "i\t" << time_count(10, type, 0)/time_count(10, type, 1) << "\t" << time_count(100, type, 0)/time_count(100, type, 1) << "\t" << time_count(200, type, 0) / time_count(200, type, 1) << "\t" << time_count(400, type, 0) / time_count(400, type, 1) << "\t" << time_count(800, type, 0) / time_count(800, type, 1) << "\t" << time_count(1000, type, 0) / time_count(1000, type, 1) << "\t" << time_count(2000, type, 0) / time_count(2000, type, 1) << endl;
	//-------------------------
	type = 'f';
	cout << "f\t" << time_count(10, type, 0)/time_count(10, type, 1) << "\t" << time_count(100, type, 0)/time_count(100, type, 1) << "\t" << time_count(200, type, 0) / time_count(200, type, 1) << "\t" << time_count(400, type, 0) / time_count(400, type, 1) << "\t" << time_count(800, type, 0) / time_count(800, type, 1) << "\t" << time_count(1000, type, 0) / time_count(1000, type, 1) << "\t" << time_count(2000, type, 0) / time_count(2000, type, 1) << endl;
	//-------------------------
	type = 'd';
	cout << "d\t" << time_count(10, type, 0)/time_count(10, type, 1) << "\t" << time_count(100, type, 0)/time_count(100, type, 1) << "\t" << time_count(200, type, 0) / time_count(200, type, 1) << "\t" << time_count(400, type, 0) / time_count(400, type, 1) << "\t" << time_count(800, type, 0) / time_count(800, type, 1) << "\t" << time_count(1000, type, 0) / time_count(1000, type, 1) << "\t" << time_count(2000, type, 0) / time_count(2000, type, 1) << endl;
	//-------------------------
	int buf;
	cin >> buf;
	return 0;
}
