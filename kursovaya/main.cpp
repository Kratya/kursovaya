#include "Galerkin.h"
#define USE_SECOND_CONDITION
int main()
{
	vector<double>x(3), y(3);
	Input();
	Portrait();

	di.resize(n);
	global_b.resize(n);
	ggl.resize(ig[n] - ig[0]);
	ggu.resize(ig[n] - ig[0]);

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int point = elems[i][j];
			x[j] = nodes[point][0];
			y[j] = nodes[point][1];
		}
		Local_Matrix(x, y, local_b);
		AddLocalToGlobal(elems[i]);
	}
#ifdef USE_SECOND_CONDITION
	uc_kraev2();
#endif
	uc_kraev1();
	SoprGrad();
	//LOS();
}
