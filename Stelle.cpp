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
	return 	pow(y/(k*(g - 1)),1/g);
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
	
	double Ps = 1e-6; //Pressione centrale
	double init[4] = {0.,1e5 / R0_m,0.,Ps};//Vettore dei parametri iniziali (r0,r1,m iniziale,pressione iniziale)
	double r = 4e3;
	ofstream file("./DatiGamma5-3.csv");//File salvataggio dati
	while (r >= 3e3)
	{
		RK4(f, g, parametri, init, 10000);
		if (init[0] * R0_m >= 3e3 && init[0] * R0_m <= 5e4) 
		{
			file << init[0] * R0_m << "," << init[2] * (M0 / Ms) << endl;
		}
		r = init[0] * R0_m;
		init[0] = 0., init[2] = 0., init[3] = Ps;
		Ps *= 1.01;
		cout << Ps << endl;

	}
	file.close();
	//Gas Fermioni relativistico
	ofstream file2("./DatiGamma4-3.csv");
	parametri[0] = 4. / 3.;
	parametri[1] = 0.1;
	Ps = 1e-6;
	init[0] = 0., init[2] = 0., init[3] = Ps;
	r = 5e4;
	cout << "4/3" << endl;
	while (r >= 3e3 && r <= 5e4 )
	{

		RK4(f, g, parametri, init, 10000);
		if (init[0] * R0_m >= 2.9e3 && init[0] * R0_m <= 5e4)
		{
			file2 << init[0] * R0_m << "," << init[2] * (M0 / Ms) << endl;
			r = init[0] * R0_m;
		}
		init[0] = 0., init[2] = 0., init[3] = Ps;
		Ps *= 1.01;
		cout << Ps << endl;

	}
	file2.close();
	//Stella di Neutroni
	ofstream file3("./DatiGamma2,5.csv");
	parametri[0] = 2.54;
	parametri[1] = 0.01;
	Ps = 1e-6;
	init[0] = 0., init[2] = 0., init[3] = Ps;
	r = 1e4;
	cout << "Neutroni" << endl;
	double stabilita;
	while (r <= 5e4)
	{

		RK4(f, g, parametri, init, 10000);
		if (init[0] * R0_m >= 3e3 && init[0] * R0_m <= 5.1e4)
		{	
		    stabilita= pow(init[0] * R0_m, 3. * Gamma - 4.) * pow(init[2] * (M0 / Ms), 2. - Gamma);
			file3 << init[0] * R0_m << " " << init[2] * (M0 / Ms) <<" " << stabilita << endl;
			cout << init[0] * R0_m << " " << init[2] * (M0 / Ms) <<" " << stabilita << endl;
			
			r = init[0] * R0_m;
		}

		init[0] = 0., init[2] = 0., init[3] = Ps;
		Ps *= 1.01;
		cout << Ps << endl;

	}
	file3.close();


	return 0;
}

