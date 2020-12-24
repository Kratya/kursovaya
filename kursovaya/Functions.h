#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<iomanip>
#include<cmath>
#include<string>
#include<algorithm>
#include<iterator>
#include<conio.h>
#define test1
using namespace std;

int n, m;
vector <vector<double>> nodes;
vector<vector<int>> elems, conds1, conds2, conds3;

vector<vector<double>> local_matrix(3), G_matrix(3), M_matrix(3);
vector<double> local_b(3);

vector<vector<double>> global_A; 
vector<double> global_b;

vector<int> ig, jg;
vector<double> ggl, di;

#ifdef test1
double function(double x, double y)
{
	return 0;
}

double Diffusion_coef()
{
	return 2;
}

double Gamma_coef()
{
	return 3;
}

double Kraev_us1(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return 5;
	case 1:
		return 0;
	default: 0;
	}
}

double Kraev_us2(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return 2 * Diffusion_coef() * y;
	case 1:
		return -3 * Diffusion_coef();
	}
}
#endif

#ifdef test2
double function(double x, double y)
{
	return 4*(x + y);
}

double Diffusion_coef()
{
	return 2;
}

double Gamma_coef()
{
	return 4;
}

double Kraev_us1(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return x + y;
	case 1:
		return 0;
	default: 0;
	}
}

double Kraev_us2(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return 1;
	case 1:
		return 1;
	default:
		return 0;
	}
}
#endif 

#ifdef test3
double function(double x, double y)
{
	return x * x;
}

double Diffusion_coef()
{
	return 2;
}

double Gamma_coef()
{
	return 3;
}

double Kraev_us1(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return x * x;
	case 1:
		return 0;
	default: 0;
	}
}

double Kraev_us2(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return 0;
	case 1:
		return 0;
	default:
		return 0;
	}
}
#endif 

#ifdef test4
double function(double x, double y)
{
	return 5;
}

double Diffusion_coef()
{
	return 2;
}

double Gamma_coef()
{
	return 3;
}

double Kraev_us1(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return 5;
	case 1:
		return 0;
	default: 0;
	}
}

double Kraev_us2(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return 0;
	case 1:
		return 0;
	default:
		return 0;
	}
}
#endif 

#ifdef test5
double function(double x, double y)
{
	return x * y;
}

double Diffusion_coef()
{
	return 1;
}

double Gamma_coef()
{
	return 1;
}

double Kraev_us1(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return x * y;
	case 1:
		return 0;
	default: 0;
	}
}

double Kraev_us2(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return 0;
	case 1:
		return 0;
	default:
		return 0;
	}
}
#endif 

double mes_G(double x1, double x2, double y1, double y2)
{
	return (sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)));
}

void uc_kraev1()
{
	ifstream fin("FirstCondition.txt");
	int p, v1, v2, num_f;		//кол-во ребер, на которых задано 1 кр.условие
	long B = 1e+10;
	fin >> p;
	for (int i = 0; i < p; i++)
	{
		fin >> v1 >> v2 >> num_f;
		double x1, y1, x2, y2;
		x1 = nodes[v1][0];
		y1 = nodes[v1][1];
		x2 = nodes[v2][0];
		y2 = nodes[v2][1];
		double res1 = Kraev_us1(num_f, x1, y1);
		double res2 = Kraev_us1(num_f, x2, y2);
		di[v1] = B;
		di[v2] = B;
		global_b[v1] = B * res1;
		global_b[v2] = B * res2;
	}
	fin.close();
}



void uc_kraev2()
{
	ifstream fin("SecondCondition.txt");
	int p, v1, v2, num_f;
	fin >> p;
	for (int i = 0; i < p; i++)
	{
		fin >> v1 >> v2 >> num_f;
		double x1, y1, x2, y2;
		x1 = nodes[v1][0];
		y1 = nodes[v1][1];
		x2 = nodes[v2][0];
		y2 = nodes[v2][1];
		double mesG = mes_G(x1, x2, y1, y2) / 6.;
		double tetta1 = Kraev_us2(num_f, x1, y1);
		double tetta2 = Kraev_us2(num_f, x2, y2);
		global_b[v1] += mesG * (2 * tetta1 + tetta2);
		global_b[v2] += mesG * (tetta1 + 2 * tetta2);
	}
	fin.close();
}