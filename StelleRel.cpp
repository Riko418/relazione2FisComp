#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "EquazioniDifferenziali.h"
using namespace std;


double rho(double t, double x, double y,double* args)
{
	double g = *args;
	double k = *(args + 1);
	return 	pow(y/(k*(g - 1.)),1./g);
}
double epsilon(double t, double x, double y, double* args)
{
	double g = *args;
	double k = *(args + 1);
	return 	rho(t,x,y,args)+k*pow(rho(t, x, y, args),g);
}
double frel(double t, double x, double y, double* fargs)
{
	return epsilon(t,x,y,fargs)*t*t;
}
double grel(double t, double x, double y, double* gargs)
{
	if (t == 0.)
	{
		return 0;
	}
	else
	{
		return  -((y+ epsilon(t, x, y, gargs))*(x+t*t*t*y))/(t*t-2.*x*t);
	}
}

int main()
{
	ofstream file("./GammaRel4-3.csv");
	double pi = 3.141592653589793;
	double hc = 197.327;
	double Ms = 1.1157467e60;
	double Mnc = 939.565;
	double rhoS = 0.16 * Mnc;
	double G = hc * 6.67259e-45;
	double Rn = 10.;
	double P0 = rhoS;
	double Gamma = 4./3.;
	double K = 0.1;
	double parametri[2] = {Gamma,K};
	
    double R0 = sqrt(P0 / (4. * pi * G *rhoS * rhoS ));
	double R0_m = R0 * 1e-15; //R0 in metri
	double M0 = P0 * R0 / (G * rhoS);
	
	double Ps = 1e-6; //Pressione centrale
	double init[4] = {0,2e5 / R0_m,0,Ps};
	
	double r = 5e4;
	//cout<< 1e5 / (R0_m * 100000)<<endl;
	
	while (r >= 3e3)
	{

		RK4(frel, grel, parametri, init, 10000);
		if (init[0] * R0_m >= 3e3 && init[0] * R0_m <= 1.9e5)
		{
			cout<< init[0] * R0_m << " " << init[2] * (M0 / Ms) <<" " <<Ps<< endl;
			file << init[0] * R0_m << " " << init[2] * (M0 / Ms) <<" " <<Ps<< endl;
		}
		r = init[0] * R0_m;
		init[0] = 0, init[2] = 0, init[3] = Ps;
		Ps *= 1.001;
		cout << Ps << endl;

	}
	file.close();
	
	return 0;
}

