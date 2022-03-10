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
#define test5
using namespace std;

int n, m;
vector <vector<double>> nodes;
vector<vector<int>> elems, conds1, conds2, conds3;

vector<vector<double>> local_matrix(3), G_matrix(3), M_matrix(3);
vector<double> local_b(3);

vector<vector<double>> global_A; 
vector<double> global_b;

vector<int> ig, jg;
vector<double> ggl, di, ggu;

#ifdef test1
double function(double x, double y)
{
	return 0;
}

double Diffusion_coef()
{
	return 1;
}

double Gamma_coef()
{
	return 0;
}

double Kraev_us1(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return x;
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
		return  Diffusion_coef() * 0;
	case 1:
		return -1*Diffusion_coef();
	case 2:
		return 1 * Diffusion_coef();
	}
}
#endif

#ifdef test2
double function(double x, double y)
{
	return 0;
}

double Diffusion_coef()
{
	return 1;
}

double Gamma_coef()
{
	return 0;
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
		return  Diffusion_coef() * 0;
	case 1:
		return -1 * Diffusion_coef();
	case 2:
		return 1 * Diffusion_coef();
	}
}
#endif 

#ifdef test3
double function(double x, double y)
{
	return 0;
}

double Diffusion_coef()
{
	return 1;
}

double Gamma_coef()
{
	return 0;
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
		return  Diffusion_coef() * 0;
	case 1:
		return -2*x * Diffusion_coef();
	case 2:
		return 2*x * Diffusion_coef();
	}
}
#endif 


#ifdef test4
double function(double x, double y)
{
	return 0;
}

double Diffusion_coef()
{
	return 1;
}

double Gamma_coef()
{
	return 0;
}

double Kraev_us1(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return x*x + y*y;
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
		return  Diffusion_coef() * 0;
	case 1:
		return -2 * x * Diffusion_coef();
	case 2:
		return 2 * x * Diffusion_coef();
	case 3:
		return -2 * y * Diffusion_coef();
	case 4:
		return 2 * y * Diffusion_coef();
	}
}
#endif 

#ifdef test5
double function(double x, double y)
{
	return 0;
}

double Diffusion_coef()
{
	return 1;
}

double Gamma_coef()
{
	return 0;
}

double Kraev_us1(int num_f, double x, double y)
{
	switch (num_f)
	{
	case 0:
		return x * x * x;
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
		return  Diffusion_coef() * 0;
	case 1:
		return -3 * x*x * Diffusion_coef();
	case 2:
		return 3 * x*x * Diffusion_coef();
	}
}
#endif 


double mes_G(double x1, double x2, double y1, double y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
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
		/*
		for (int k = ig[v1]; k < ig[v1 + 1]; k++)
		{
			ggl[k] = 0;
		
		}

		for (int k = ig[v2]; k < ig[v2 + 1]; k++)
		{
			ggl[k] = 0;

		}

		for (int k = v1+1; k < n; k++)
		{
			for (int l = ig[v1]; l < ig[v1 + 1]; l++)
			{
				if (jg[l] == v1)
				{
					ggu[l] = 0;
					break;
				}

				if (jg[l] > v1)
					break;

			}

		}

		for (int k = v2 + 1; k < n; k++)
		{
			for (int l = ig[v2]; l < ig[v2 + 1]; l++)
			{
				if (jg[l] == v2)
				{
					ggu[l] = 0;
					break;
				}

				if (jg[l] > v2)
					break;

			}

		}*/



	
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