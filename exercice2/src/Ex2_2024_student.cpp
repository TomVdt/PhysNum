#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <valarray>
#include <cmath> // Se usi fun
using namespace std;

class Exercice2
{

private:
	double t, dt, tFin;
	double m, g, L, alpha;
	double d, Omega, kappa;
	double theta, thetadot;
	int N_excit, nsteps_per;
	int sampling;
	int last;
	ofstream* outputFile;

	void printOut(bool force)
	{
		if ((!force && last >= sampling) || (force && last != 1))
		{
			// TODO: vérifier
			double emec = Emec(theta, thetadot, t);
			double pnc = Pnonc(theta, thetadot, t);

			*outputFile << t << " "
						<< theta << " "
						<< thetadot << " "
						<< emec << " "
						<< pnc << " "
						<< length(t) << endl;
			last = 1;
		}
		else
		{
			last++;
		}
	}

	// TODO: vérifier
	valarray<double> acceleration(double x, double v, double t_)
	{
		valarray<double> acc = valarray<double>(2);
		const double l = length(t_);

		acc[0] = -g * sin(x) / l; // accélération dépendant uniquement de x
		acc[1] = -2 * lendot(t_) * v / l; // accélération dépendant uniquement de v

		return acc;
	}

	// TODO: vérifier
	double Emec(double x, double v, double t_)
	{
		const double l = length(t_);
		return  1.0 / 2.0 * m * pow(l * v, 2) + m * g * l * cos(x);
	}

	// TODO: FAUX
	double Pnonc(double x, double v, double t_)
	{
		return (-m * length(t_) * v*v - m * g * cos(x) * lendot(t_));
	}

	double length(double t_)
	{
		return L + alpha * t + d * sin(Omega * t_);
	}

	double lendot(double t_)
	{
		return alpha + Omega * d * cos(Omega * t_);
	}

	double lendotdot(double t_)
	{
		return -Omega * Omega * d * sin(Omega * t_);
	}

	void step()
	{
		// TODO: vérfier
		valarray<double> a = acceleration(theta, thetadot, t);
		const double dt2 = dt*dt;
		const double a_tot = a.sum();

		theta = theta + thetadot * dt + 1.0 / 2.0 * a_tot * dt2;

		const double v_plus_dt_2 = thetadot + 1.0 / 2.0 * a_tot;
		// v doesn't matter for this, explicitly mark it to avoid ~confusion~
		const double new_a1 = acceleration(theta, 0.0, t + dt)[0];
		// x doesn't matter for this, explicitly mark it to avoid ~confusion~
		const double new_a2 = acceleration(0.0, v_plus_dt_2, t + dt / 2.0)[1];

		thetadot = thetadot + 1.0 / 2.0 * (a[0] + new_a1) * dt + new_a2 * dt;
	}


public:
	Exercice2(int argc, char* argv[])
	{
		const double pi = 3.1415926535897932384626433832795028841971e0;
		string inputPath("configuration.in"); // Fichier d'input par defaut
		if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
			inputPath = argv[1];

		ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

		for (int i(2); i < argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
			configFile.process(argv[i]);

		tFin = configFile.get<double>("tFin");      // t final (réécrit si N >0)
		d = configFile.get<double>("d");         // amplitude du changement periodic de longueur
		Omega = configFile.get<double>("Omega");     // fréquence du changement periodic de longueur 
		kappa = configFile.get<double>("kappa");     // coefficient de frottement
		m = configFile.get<double>("m");         // mass
		g = configFile.get<double>("g");         // accélération de la gravité
		L = configFile.get<double>("L");         // longueur
		alpha = configFile.get<double>("alpha");     // Coefficient linéaire du changement de longueur
		theta = configFile.get<double>("theta0");    // condition initiale en theta
		thetadot = configFile.get<double>("thetadot0"); // condition initiale en thetadot
		sampling = configFile.get<int>("sampling");     // fréquence d'écriture des données
		N_excit = configFile.get<int>("N");            // Nombre de périodes d'excitation simulées
		nsteps_per = configFile.get<int>("nsteps");       // nombre de pas de temps par période (si N>0), nombre total de pas de temps si N=0

		// Ouverture du fichier de sortie
		outputFile = new ofstream(configFile.get<string>("output").c_str());
		outputFile->precision(15);

		if (N_excit > 0) {
			// TODO Définir le pas temporel et le temp final si N_excit est défini.
			tFin = 10.0;
			dt = 1.0;
		}
		else {
			// TODO: vérifier
			dt = tFin / nsteps_per;
		}
	}

	~Exercice2()
	{
		outputFile->close();
		delete outputFile;
	};

	void run()
	{
		t = 0.;
		last = 0;
		printOut(true);

		while (t < tFin - 0.5 * dt)
		{
			step();

			t += dt;
			printOut(false);
		}
		printOut(true);
	};

};

int main(int argc, char* argv[])
{
	Exercice2 engine(argc, argv);
	engine.run();

	return 0;
}
