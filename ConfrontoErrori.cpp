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
	return 	pow(y/(k*(g - 1)),1./g);
}

double f(double t, double x, double y, double* fargs)
{
	return rho(t,x,y,fargs)*t*t;
}
double g(double t, double x, double y, double* gargs)
{
	if (t == 0)
	{
		return 0;
	}
	else
	{
		return  (-x * rho(t, x, y, gargs)) / t / t;
	}
}
int main()
{
	//Costanti
	double pi = 3.141592653589793;
	double hc = 197.327;
	double Ms = 1.1157467e60;
	double Mnc = 939.565;
	double rhoS = 0.16 * Mnc;
	double G = hc * 6.67259e-45;
	double Rn = 10.;
	double P0 = rhoS;
	double Gamma = 2.54;
	double K = 0.01;
	double parametri[2] = {Gamma,K};
	
	//Calcolo R0 e M0 adimensionali
    double R0 = sqrt(P0 / (4 * pi * G * rhoS * rhoS));
	double R0_m = R0 * 1e-15; //R0 in metri
	double M0 = P0 * R0 / (G * rhoS);
	
	double Ps = 10; //Pressione centrale
	double init[4] = {0,1e5 / R0_m,0,Ps};//Vettore dei parametri iniziali (r0,r1,m iniziale,pressione iniziale)

	RK4(f, g, parametri, init, 1e7); //Calcolo la soluzione con un h molto piccolo
	double Raggio = init[0];
	double Massa = init[2];
	double h = 1e5 / (R0_m * 1e7);
	double N = 1.;
	//cout << Raggio << " " << Massa << " " << h << endl;
	init[0] = 0, init[2] = 0, init[3] = Ps;

	ofstream file("./Errori.csv");
	for (int i = 0; i < 1500; i++)
	{
		N *= 1.01;
		h = 1e5 / (R0_m * int(N));
		RK4(f, g, parametri, init, int(N));
		file << abs(init[0] - Raggio) / Raggio << "," << abs(init[2] - Massa) / Massa << ",";
		init[0] = 0., init[2] = 0., init[3] = Ps;
		EuleroEsplicito(f, g, parametri, init, int(N));
		file << abs(init[0] - Raggio) / Raggio << "," << abs(init[2] - Massa) / Massa << "," << h << endl;
		init[0] = 0., init[2] = 0., init[3] = Ps;

	}
	cout << "Fine" << endl;

	return 0;
}

