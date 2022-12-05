#include"opt_alg.h"
double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		//Tu wpisz kod funkcji

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	
	solution Xopt;
        
	solution XB(x0), XB_old, X;
	XB.fit_fun(ff, ud1, ud2);
	//output_HJ << t(XB.x(0)) << t(XB.x(1)) << endl;    //zapis iteracje
	while (true)
	{
		X = HJ_trial(ff, XB, s, ud1, ud2); //odpalenie etapu probnego
		if (X.y < XB.y) //sprawdzamy czy etap probny przyniosl poprawe
			{
			while (true) //etap roboczy wykonywany co chwile
				{
				XB_old = XB;
				XB = X;
				X.x = 2 * XB.x - XB_old.x;
				X.fit_fun(ff, ud1, ud2);
				X = HJ_trial(ff, X, s, ud2, ud2);
				if (X.y >= XB.y)
					break; //przerwanie etapu roboczego
				if (solution::f_calls > Nmax)
					return XB;
				}
			}
		else //zmniejszamy dlugosc kroku
			s *= alpha;
		if (s<epsilon || solution::f_calls>Nmax) { //warunki stopu
			XB.flag = 0;
			return XB;
		}
		//output_HJ << t(XB.x(0)) << t(XB.x(1)) << endl;    //zapis iteracje
	}
	Xopt.flag = 1;
	return Xopt;
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	int n = get_dim(XB); //dlugosc wektora X
	matrix D = ident_mat(n); //kierunki macierz jednostkowa
	//Etap probny konczy sie porazka gdy zostajemy na podstawowym punkcie, nastyepnie zmniejszamy kroki az krok bedzie mniejszy
	// od epsilon
	//Po zakonczeniu etapu probnego, jesli znalezlismy nowy punkt wykonujemy etap roboczy, odbicie lustrzane starej bazy wzgledem nowej bazy
	// z punkty X (odbicia) odpalamy etap probny, to co zwroci  porownujemy z punktem symetrii, jesli jest lepszy wykonujemy kolejny raz etap roboczy
	// jesli jest gorszy anulujemy etap roboczy , wracamy do bazy i rozpoczynamy iteracje 
	solution X;
	for (int i = 0; i < n; ++i)
	{
		X.x = XB.x + s * D[i];
		X.fit_fun(ff, ud1, ud2);
		if (X.y < XB.y)
			XB = X;
		else
		{
			X.x = XB.x - s * D[i];
			X.fit_fun(ff, ud1, ud2);
			if (X.y < XB.y)
				XB = X;
		}
	}
	return XB;
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
