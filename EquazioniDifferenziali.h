#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;


void EuleroEsplicito(double (*f)(double, double, double, double*), double (*g)(double, double, double, double*), double* args, double* init, int res)
{

	double t = *(init), x = *(init + 2), y = *(init + 3);
	double h = (*(init + 1) - *(init)) / (double)res;
	double xprec, yprec;
	for (int i = 0; i < res; i++)
	{
		yprec = y;
		xprec = x;
		x += h * f(t, x, y, args);
		y += h * g(t, xprec, y, args);
		t += h;
		if (isnan(y))
		{
			*init = t - h;
			*(init + 2) = xprec;
			*(init + 3) = yprec;
			break;
		}
	}
	if (!isnan(y))
	{
		*init = t;
		*(init + 2) = x;
		*(init + 3) = y;
	}
}
//------------------------
void RK4(double (*f)(double, double, double, double*), double (*g)(double, double, double, double*), double* args, double* init, int res)
{

	double t = *(init), x = *(init + 2), y = *(init + 3);
	double h = (*(init + 1) - *(init)) / (double)res;
	double k[4], l[4];
	double xprec, yprec;
	for (int i = 0; i < res; i++)
	{

		
		xprec = x;
		yprec = y;
		k[0] = f(t, x, y, args);
		l[0] = g(t, x, y, args);
		k[1] = f(t + h / 2, x + (k[0] * h) / 2, y + (l[0] * h) / 2, args);
		l[1] = g(t + h / 2, x + (k[0] * h) / 2, y + (l[0] * h) / 2, args);

		k[2] = f(t + h / 2, x + (k[1] * h) / 2, y + (l[1] * h) / 2, args);
		l[2] = g(t + h / 2, x + (k[1] * h) / 2, y + (l[1] * h) / 2, args);

		k[3] = f(t + h, x + (k[2] * h), y + (l[2] * h), args);
		l[3] = g(t + h, x + (k[2] * h), y + (l[2] * h), args);
		x += (h / 6) * (k[0] + 2 * k[1] + 2 * k[2] + k[3]);
		y += (h / 6) * (l[0] + 2 * l[1] + 2 * l[2] + l[3]);
		t += h;
		if (isnan(y))
		{
			*init = t - h;
			*(init + 2) = xprec;
			*(init + 3) = yprec;
			break;
		}
	}
	if (!isnan(y))
	{
		*init = t;
		*(init + 2) = x;
		*(init + 3) = y;
	}

}



//------------------------
/*
void AdaptiveEulero(double (*f)(double, double, double, double*), double (*g)(double, double, double, double*), double* args, double* init,double tol)
{

	double t = *(init), x = *(init + 2), y = *(init + 3);
	double h = 1;
	double xprec, yprec;
     while (t < *(init+1))
	{
		yprec = y;
		xprec = x;
		double x1 = x + h * f(t,x, y,args);
		double x2 = x + h * f(t+h, x1, y, args);
		double y1 = y + h * g(t, xprec, y, args);
		double y2 = y + h * g(t + h, xprec, y1, args);
		if (isnan(x1) || isnan(x2) || isnan(y1) || isnan(y2))
		{
			h /= 2;
		}
		else
		{
			double err = abs((x2 - x1) / (1 - pow(2, -1)));

			if (err < tol) {
				x = x2;
				y = y2;
				t += h;
			}

			h *= pow(tol / err, 0.5);

			if (isnan(y))
			{
				*init = t - h;
				*(init + 2) = xprec;
				*(init + 3) = yprec;
				break;
			}
		}
		
	}
	if (!isnan(y))
	{
		*init = t;
		*(init + 2) = x;
		*(init + 3) = y;
	}
}

//------------------------
void RK4Step(double (*f)(double, double, double, double*), double (*g)(double, double, double, double*), double* args,double* xy,double t,double x,double y,double h) 
{
	double k[4], l[4];
	k[0] = f(t, x, y, args);
	l[0] = g(t, x, y, args);
	k[1] = f(t + h / 2, x + (k[0] * h) / 2, y + (l[0] * h) / 2, args);
	l[1] = g(t + h / 2, x + (k[0] * h) / 2, y + (l[0] * h) / 2, args);
	k[2] = f(t + h / 2, x + (k[1] * h) / 2, y + (l[1] * h) / 2, args);
	l[2] = g(t + h / 2, x + (k[1] * h) / 2, y + (l[1] * h) / 2, args);
	k[3] = f(t + h, x + (k[2] * h), y + (l[2] * h), args);
	l[3] = g(t + h, x + (k[2] * h), y + (l[2] * h), args);

	*(xy) =x+(h / 6) * (k[0] + 2 * k[1] + 2 * k[2] + k[3]);
	*(xy+1) =y+(h / 6) * (l[0] + 2 * l[1] + 2 * l[2] + l[3]);
}


void AdaptiveRK4(double (*f)(double, double, double, double*), double (*g)(double, double, double, double*), double* args, double* init, double tol)
{
		double t = *(init), x = *(init + 2), y = *(init + 3);
		double h = 1;
		double xprec, yprec,err;
		double xy1[2], xy2[2], xy2h[2];
		while (t < *(init + 1))
		{
			yprec = y;
			xprec = x;

			RK4Step(f, g, args, xy1, t, x, y, h);
			RK4Step(f, g, args, xy2, t, x, y, h / 2);
			RK4Step(f, g, args, xy2h, t+h/2, xy2[0], xy2[1], h / 2);
		
			if (isnan(xy1[0]) || isnan(xy1[1]) || isnan(xy2[0]) || isnan(xy2[1]) || isnan(xy2h[0]) || isnan(xy2h[1]))
			{
				h /= 2;
			}
			else
			{
				err = abs((xy2h[0] - xy1[0]) / (1 - pow(2, -1)));
				if (err < tol) {

					t += h;
					x = xy2h[0];
					y = xy2h[1];
				}
				h *= pow(tol / err, 0.25);

			}
			if (err==0)
			{
				*init = t - h;
				*(init + 2) = xprec;
				*(init + 3) = yprec;
				break;
			}
		}
		if (!isnan(y))
		{
			*init = t;
			*(init + 2) = x;
			*(init + 3) = y;
		}


}
*/