#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <valarray>
#include <cmath> 
#include <numeric>
using namespace std;

class Exercice3
{

private:
	long double t, dt, tFin;
	long double r0, r1, v0, a;
	long double xs, xt;
	long double G = 6.674e-11l;
	long double mtot = 0e0l;
	long double L2x, L2y;
	valarray<long double> m = std::valarray<long double>(0.e0l, 2);
	// int N_excit;
	int nsteps;
	int sampling;
	int last;
	int  nsel_physics;
	bool adapt;
	long double alpha = 0e0l;
	long double beta = 0e0l;
	long double tol = 0e0l;
	valarray<long double> x0 = valarray<long double>(0.e0l, 4); // Correctly initialized
	valarray<long double> x = valarray<long double>(0.e0l, 4); // Correctly initialized
	ofstream* outputFile;

	long double omega;

	void printOut(bool write) {
		if ((!write && last >= sampling) || (write && last != 1)) {
			long double Energy = compute_energy(x);
			*outputFile
				<< t << " "
				<< x[0] << " "
				<< x[1] << " "
				<< x[2] << " "
				<< x[3] << " "
				<< Energy << " "
				<< nsteps << " "
				<< dt << "\n";
			last = 1;
		}
		else {
			last++;
		}
	}

	valarray<long double> get_f(const valarray<long double>& x, long double t) {
		valarray<long double> xdot(0.0, 4);
		if (nsel_physics == 1) {

			const long double prefact = -G * m[1] / pow(norm(x), 3);

			xdot[2] = x[0] * prefact;
			xdot[3] = x[1] * prefact;
		}
		else if (nsel_physics == 2) {
			const long double Rs = sqrt(x[1] * x[1] + (x[0] - xs) * (x[0] - xs));
			const long double Rt = sqrt(x[1] * x[1] + (x[0] - xt) * (x[0] - xt));

			const long double c1 = G * m[0] / (Rs * Rs * Rs);
			const long double c2 = G * m[1] / (Rt * Rt * Rt);

			xdot[2] = -c1 * (x[0] + alpha * a)
				- c2 * (x[0] - beta * a)
				+ omega * omega * x[0]
				+ 2.0l * omega * x[3];
			xdot[3] = -c1 * x[1]
				- c2 * x[1]
				+ omega * omega * x[1]
				- 2.0l * omega * x[2];
		}
		else {
			cerr << "No dynamics corresponds to this index" << endl;
			return xdot;
		}

		// Velocities
		xdot[0] = x[2];
		xdot[1] = x[3];

		return xdot;
	}


	// Function to compute potential energy per mass in R (nsel_physics=1) or in R'(nsel_physics=2)
	long double get_Epot(const valarray<long double>& x) {
		long double energy_pot = 0;
		if (nsel_physics == 1) {
			energy_pot = -G * m[1] / norm(x);
		}
		else if (nsel_physics == 2) {
			energy_pot = 
				- 1.0 / 2.0 * omega * omega * (x[0] * x[0] + x[1] * x[1])
				- G * m[0] / sqrt(pow(x[0] - xs, 2) + x[1] * x[1])
				- G * m[1] / sqrt(pow(x[0] - xt, 2) + x[1] * x[1]);
		}
		return energy_pot;
	}

	// Function to compute mechanical energy per mass in R'
	long double compute_energy(const valarray<long double>& x) {
		return 1.0 / 2.0 * (x[2] * x[2] + x[3] * x[3])
			// + omega * (x[3] * x[0] - x[2] * x[1])
			+ get_Epot(x);
	}

	// Norm of the position, does not take the velocity
	long double norm(const valarray<long double>& vect) {
		return sqrt(vect[0] * vect[0] + vect[1] * vect[1]);
	}


	void initial_condition(void) {
		if (nsel_physics == 1) {			
			x0[0] = -r0;
			x0[1] = 0.0l;

			// To calculate the kinetic energy
			// omega = 0;

			v0 = r1 * sqrt(2.0l * G * m[1] * (1.0l / r0 - 1.0l / r1) / (r1 * r1 - r0 * r0));
			x0[2] = 0.0;
			x0[3] = v0;
		}
		else if (nsel_physics == 2) {
			// TODO go back through this to check correct initialisation
			xs = -m[1] * a / mtot;
			xt = m[0] * a / mtot;
			omega = sqrt(G * m[0] / (a * a * xt));

			x0[0] = L2x;
			x0[1] = L2y;

			x0[2] = 0.0l;
			x0[3] = -0.1l;

		}
		else {
			cerr << "No dynamics corresponds to this index" << endl;
		}
	}

	valarray<long double> RK4_do_onestep(const valarray<long double>& y_old, long double t, long double dt) {
		valarray<long double> k1(0.0, 4), k2(0.0, 4), k3(0.0, 4), k4(0.0, 4), ynew(0.0, 4);

		k1 = dt * get_f(y_old, t);
		k2 = dt * get_f(y_old + k1 / 2.0, t + dt / 2.0);
		k3 = dt * get_f(y_old + k2 / 2.0, t + dt / 2.0);
		k4 = dt * get_f(y_old + k3, t + dt);

		ynew = y_old + 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		return ynew;
	}


public:
	Exercice3(int argc, char* argv[]) {
		const long double pi = 3.1415926535897932384626433832795028841971e0l;
		string inputPath("configuration.in"); // Fichier d'input par defaut
		if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
			inputPath = argv[1];

		ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

		for (int i(2); i < argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
			configFile.process(argv[i]);

		tFin = configFile.get<long double>("tFin");            // t final (overwritten if N_excit >0)
		//  dt       = configFile.get<long double>("dt");            // time step (overwritten if nsteps_per >0)
		m[0] = configFile.get<long double>("m1");              // mass of the sun
		m[1] = configFile.get<long double>("m2");              // mass of the earth
		r0 = configFile.get<long double>("r0");              // r0
		r1 = configFile.get<long double>("r1");              // r1
		L2x = configFile.get<long double>("L2x");              // L2x
		L2y = configFile.get<long double>("L2y");              // L2y
		a = configFile.get<long double>("a");               // demi grand-axe (solei-terre hyp MCU)
		nsel_physics = configFile.get<int>("nsel_physics");       //1) one body problem around mass#2 or 2) one body in rotating reference frame of {1,2}
		adapt = configFile.get<bool>("adapt");             //if 1=true -> adaptive dt, if 0=false -> fixed dt
		tol = configFile.get<long double>("tol");             //tolerance of the adaptive scheme
		sampling = configFile.get<int>("sampling");     // number of time steps between two writings on file
		nsteps = configFile.get<int>("nsteps");        // number of time step per period
		mtot = m[0] + m[1];
		alpha = m[1] / mtot;
		beta = m[0] / mtot;


		// Ouverture du fichier de sortie
		outputFile = new ofstream(configFile.get<string>("output").c_str());
		outputFile->precision(20);
		//TO DO initialize tFin for nsel_physics=1 and initialize dt for both nsel_physics
		if (nsel_physics == 1) {
			a = (r0 + r1)/2.0l;
			tFin = 2.0l * pi * sqrt(pow(a, 3.0l) / (G * m[1])) / 2.0l;
		}

		dt = tFin / nsteps;
	}

	~Exercice3() {
		outputFile->close();
		delete outputFile;
	};

	void run() {
		t = 0.0l;
		initial_condition();
		x = x0;
		last = 0;
		printOut(true);
		valarray<long double> y1(0.0l, 4);
		valarray<long double> y2(0.0l, 4);
		valarray<long double> y_tilde(0.0l, 4);

		if (!adapt) {
			while (t < tFin - dt * 0.5l) {
				x = RK4_do_onestep(x, t, dt);
				t += dt;
				printOut(false);
			}
		}
		else {
			nsteps = 0;
			long double d = 0.0;
			const long double f(0.99);

			while (t < tFin) {
				dt = min(dt, tFin - t);
				++nsteps;

				y1 = RK4_do_onestep(x, t, dt);
				y_tilde = RK4_do_onestep(x, t, dt / 2.0l);
				y2 = RK4_do_onestep(y_tilde, t, dt / 2.0l);

				d = norm(y1 - y2);

				if (d > tol) {
					do {
						dt = f * dt * pow(tol / d, 1.0l / (4.0l + 1.0l));			// TODO: RK4 converge ordre 4

						y1 = RK4_do_onestep(x, t, dt);
						y_tilde = RK4_do_onestep(x, t, dt / 2.0l);
						y2 = RK4_do_onestep(y_tilde, t, dt / 2.0l);

						d = norm(y1 - y2);
					} while (d > tol);
					x = y2;
					t += dt;
				}
				else {
					x = y2;
					t += dt;
					dt = dt * pow(tol / d, 1.0l / (4.0l + 1.0l));
					// Avoid division by 0
					// if (d != 0.0l) {
					// }
				}
				printOut(false);
			}
		}

		printOut(true); // ecrire le dernier pas de temps
	}

};


int main(int argc, char* argv[]) {
	Exercice3 engine(argc, argv);
	engine.run();

	return 0;
}
