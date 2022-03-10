#pragma once
#include "Functions.h"

vector<double> ggl_new, di_new, ggu_new, Ma;

void Mult(vector<double>& x, vector<double>& Ma)
{
	for (int i = 0; i < n; ++i) {
		int gi = ig[i], gi_1 = ig[i + 1];
		Ma[i] = di[i] * x[i];
		for (int j = gi; j < gi_1; ++j)
		{
			int column = jg[j];
			Ma[i] += ggl[j] * x[column];
			Ma[column] += ggl[j] * x[i];
		}
	}
}

double Mult_scal(vector<double> vec1, vector<double> vec2)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += vec1[i] * vec2[i];
	return sum;
}

void LOS()
{
	int k = 0;
	double norm_eps;
	vector <double>x0(n, 0);
	vector <double>r0(n, 0);
	vector <double>z0(n, 0);
	vector <double>p0(n, 0);
	vector <double>rk(n, 0);
	vector <double>zk(n, 0);
	vector <double>pk(n, 0);
	vector <double>xk(n, 0);
	Ma.resize(n);
	double maxiter = 1000;
	double eps = 1e-12;

	for (int i = 0; i < n; ++i)
	{
		x0[i] = 1;
	}
	int count = 0;
	Mult(x0, Ma);
	for (int i = 0; i < n; ++i)
	{
		r0[i] = global_b[i] - Ma[i];
		z0[i] = r0[i];
	}
	Mult(z0, p0);
	double sr = Mult_scal(r0, r0);
	while (sr > eps && count <= maxiter)
	{
		double pp = Mult_scal(p0, p0);
		double ak = Mult_scal(p0, r0) / pp;
		for (int i = 0; i < n; ++i)
		{
			x0[i] = x0[i] + ak * z0[i];
			r0[i] = r0[i] - ak * p0[i];
		}
		Mult(r0, Ma);
		double bk = -Mult_scal(p0, Ma) / pp;
		for (int i = 0; i < n; ++i)
		{
			z0[i] = r0[i] + bk * z0[i];
			p0[i] = Ma[i] + bk * p0[i];
		}
		sr = sqrt(Mult_scal(r0, r0));
		++count;

	}
	for (int i = 0; i < n; i++)
		cout << x0[i] << '\n';
}

//----------метод сопряженных градиентов---------------

void SoprGrad()
{
	cout << "Reshenie metodom sopryjennyh gradientov:" << endl;
	int k = 0;
	vector <double>y(n, 0);
	vector <double>x(n, 0);
	vector <double>ax(n, 0);
	vector <double>ap(n, 0);
	vector <double>z(n, 0);
	vector <double>z1(n, 0);
	vector <double>p(n, 0);

	double a, b, nz;
	
	double e = 1e-8;

	//  cout<<"Nevyazka na4al'nogo priblijeniya ravna:"<<endl;
	Mult(x, ax);
	for (int i = 0; i < n; i++) 
		z[i] = global_b[i] - ax[i];

	//  PrintV(z,n);
	if (Mult_scal(z, z) != 0)
	{
		for (int i = 0; i < n; i++)  
			p[i]= z[i];
		nz = 1000.;
		while (nz > e) {
			Mult(p, ap);
			a = Mult_scal(z, p) / Mult_scal(z, ap);
			for (int i = 0; i < n; i++) 
			{
				y[i] = x[i] + a * p[i];
				z1[i] = z[i] - a * ap[i];
			}
			nz = sqrt(Mult_scal(z1, z1));
			b = Mult_scal(z1, ap) / Mult_scal(p, ap);
			for (int i = 0; i < n; i++) 
			{
				p[i] = z1[i] - b * p[i];
				z[i] = z1[i];
				x[i] = y[i];
			}
			k++;
		}	
	}

	for (int i = 0; i < n; i++)
		cout << y[i] << '\n';
}


