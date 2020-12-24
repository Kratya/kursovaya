#pragma once
#include "Functions.h"
#include "LOS.h"
#include <Set>

void Input()
{
	ifstream fnodes, felems, fcond1;
	string str;

	fnodes.open("Nodes.txt");
	fnodes >> n;					
	nodes.resize(n);
	for (int i = 0; i < n; i++)
	{
		nodes[i].resize(2);
		fnodes >> nodes[i][0] >> nodes[i][1];
	}
	fnodes.close();

	felems.open("Elems.txt");
	felems >> m;					
	elems.resize(m);
	for (int i = 0; i < m; i++)
	{
		elems[i].resize(4);
		felems >> elems[i][0] >> elems[i][1] >> elems[i][2] >> elems[i][3];
	}
	felems.close();
}

void Portrait()
{
	vector<set<int>> list(n);

	for (int s = 0; s < elems.size(); s++)
		for (int p = 0; p < 3; p++)
			for (int j = p + 1; j < 3; j++)
			{
				int ind1 = elems[s][p];
				int ind2 = elems[s][j];
				if (ind1 < ind2) swap(ind1, ind2);
				list[ind1].insert(ind2);
			}

	ig.resize(n + 1);
	ig[0] = ig[1] = 0;

	for (int i = 2; i < n + 1; i++) 
	{
		int col = ig[i - 1];
		ig[i] = col + list[i - 1].size();
	}
	jg.resize(ig[n]);
	for (int i = 1, k = 0; i < n; i++) 
	{
		for (int j : list[i]) 
		{
			jg[k] = j;
			k++;
		}
	}
}

void G(vector<double>& x, vector<double>& y)
{
	double det = abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]));
	double mult = Diffusion_coef() * det / 2;
	vector<vector<double>> a(3);

	for (int i = 0; i < 3; i++)
		a[i].resize(3);

	a[0][0] = (x[1] * y[2] - y[1] * x[2]);
	a[0][1] = (-(y[2] - y[1]));
	a[0][2] = (x[2] - x[1]);

	a[1][0] = (-x[0] * y[2] + x[2] * y[0]);
	a[1][1] = (y[2] - y[0]);
	a[1][2] = (-x[2] + x[0]);

	a[2][0] = (x[0] * y[1] - x[1] * y[0]);
	a[2][1] = (-y[1] + y[0]);
	a[2][2] = (x[1] - x[0]);

	for (int i = 0; i < 3; i++)
	{
		G_matrix[i].resize(3);
		for (int j = 0; j < 3; j++)
			G_matrix[i][j] = mult * ((a[i][1]) * (a[j][1]) + (a[i][2]) * (a[j][2]));
	}
}

void M(vector<double>& x, vector<double>& y)
{
	double det = abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]));
	double* f = new double[3], mult = Gamma_coef() * det / 24;
	for (int i = 0; i < 3; i++)
	{
		M_matrix[i].resize(3);
		for (int j = 0; j < 3; j++)
			(i == j) ? M_matrix[i][j] = 2 * mult : M_matrix[i][j] = mult;
	}

	f[0] = mult * function(x[0], y[0]);
	f[1] = mult * function(x[1], y[1]);
	f[2] = mult * function(x[2], y[2]);

	local_b[0] = 2 * f[0] + f[1] + f[2];
	local_b[1] = f[0] + 2 * f[1] + f[2];
	local_b[2] = f[0] + f[1] + 2 * f[2];
}

void Local_Matrix(vector<double> x, vector<double> y, vector<double>& local_f)
{
	G(x, y);
	M(x, y);
	for (int i = 0; i < 3; i++)
	{
		local_matrix[i].resize(3);
		for (int j = 0; j < 3; j++)
			local_matrix[i][j] = M_matrix[i][j] + G_matrix[i][j];
	}

}

void AddLocalToGlobal(vector<int>elems)
{
	for (int i = 0; i < 3; i++)
	{
		di[elems[i]] += local_matrix[i][i];
		global_b[elems[i]] += local_b[i];
		for (int j = 0; j < i; j++)
		{
			auto a = elems[i];
			auto b = elems[j];
			if (a < b) swap(a, b);

			auto begin = jg.begin() + ig[a];
			if (ig[a + 1] > ig[a])
			{
				auto end = jg.begin() + ig[a + 1] - 1;
				auto iter = lower_bound(begin, end, b); 
				auto index = iter - jg.begin();
				ggl[index] += local_matrix[i][j];
			}
		}
	}

}