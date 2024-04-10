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
	double t, dt, tFin;
	double r0, r1, v0, a;
	double G = 6.674e-11;
	double mtot = 0e0;
	double L2x,L2y;
	valarray<double> m = std::valarray<double>(0.e0, 2);
	// int N_excit;
	int nsteps;
	int sampling;
	int last;
	int  nsel_physics;
	bool adapt;
	double alpha = 0e0;
	double beta = 0e0;
	double tol = 0e0;
	valarray<double> x0 = valarray<double>(0.e0, 4); // Correctly initialized
	valarray<double> x  = valarray<double>(0.e0, 4); // Correctly initialized
	ofstream *outputFile;

	// TODO this or calculate it with r0 and r1, and omega?
	double d = 149.598023e9;
	double omega;

	void printOut(bool write) {
		if((!write && last>=sampling) || (write && last!=1)) {
			double Energy = compute_energy(x[0],x[1],x[2],x[3]);
			*outputFile << t << " " << x[0] << " " << x[1] << " "<< x[2] << " " << x[3] << " " \
			<< Energy << " " << nsteps << endl; // write output on file
			last = 1;
		}
		else {
			last++;
		}
	}

	valarray<double> get_f(const valarray<double>& x, double t) {
		valarray<double> xdot(0.0, 4);
		//TO DO cheeeeeeck formulae
		if (nsel_physics == 1) {

			const double prefact = -G * m[1] / pow(x[0]*x[0] + x[1]*x[1], 3/2);

			xdot[2] = x[0] * prefact;
			xdot[3] = x[1] * prefact;
		}
		else if (nsel_physics == 2) {
		// TODO same cheeeeeck formulae

			const valarray<double> rs_vect({x[0] + r0, x[1]});
			const double rs = norm(rs_vect);

			const valarray<double> rt_vect({x[0] - r1, x[1]});
			const double rt = norm(rt_vect);

			const double prefact_s = G * m[0]/(pow(rs,3));
			const double prefact_t = G * m[1]/(pow(rt,3));


			xdot[2] = -prefact_s*(x[0] + alpha*d) - prefact_t*(x[0] - beta*d) + omega*omega*x[0] + 2*omega*x[3];
			xdot[3] = -prefact_s*x[1] - prefact_t*x[1] + omega*omega*x[1] - 2*omega*x[2];
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
	double get_Epot(double xx, double yy) {
		//TO DO
		return xx + yy;
	}

	// Function to compute mechanical energy per mass in R'
	double compute_energy(double xx, double yy, double vx, double vy) {
		//TO DO
		return vx + vy + get_Epot(xx, yy);
	}

	double norm(const valarray<double>& vect) {
		return sqrt(vect[0]*vect[0] + vect[1]*vect[1]);
	}

	// double norm4(const valarray<double>& vect) {
	// 	return sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2] + vect[3]*vect[3]);
	// }

	void initial_condition(void){
		if (nsel_physics==1) {
			x0[0] = -r0;
			x0[1] = 0.0;

			x0[2] = 0.0;		// TODO expression of v0
			x0[3] = v0;
		}
		else if (nsel_physics==2) {
			// TODO go back through this to check correct initialisation
			omega = sqrt(G*m[0] / (d*d * d*beta));


			x0[0] = L2x;
			x0[1] = L2y;

			x0[2] = 0.0;
			x0[3] = -0.1;

		}
		else {
			cerr << "No dynamics corresponds to this index" << endl;
		}
	}

	valarray<double> RK4_do_onestep(const valarray<double>& y_old, double t, double dt) {
		valarray<double> k1, k2, k3, k4, ynew;

		k1 = dt * get_f(y_old, t);
		k2 = dt * get_f(y_old + k1/2, t + dt/2);
		k3 = dt * get_f(y_old + k2/2, t + dt/2);
		k4 = dt * get_f(y_old + k3, t + dt);

		ynew = y_old + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4);
		// ynew = y_old + k1;
		return ynew;
	}


public:
	Exercice3(int argc, char* argv[]) {
		// const double pi=3.1415926535897932384626433832795028841971e0;
		string inputPath("configuration.in"); // Fichier d'input par defaut
		if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
			inputPath = argv[1];

		ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

		for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
			configFile.process(argv[i]);

		tFin         = configFile.get<double>("tFin");            // t final (overwritten if N_excit >0)
		//  dt       = configFile.get<double>("dt");            // time step (overwritten if nsteps_per >0)
		m[0]         = configFile.get<double>("m1");              // mass of the sun
		m[1]         = configFile.get<double>("m2");              // mass of the earth
		r0           = configFile.get<double>("r0");              // r0
		r1           = configFile.get<double>("r1");              // r1
		L2x          = configFile.get<double>("L2x");              // L2x
		L2y          = configFile.get<double>("L2y");              // L2y
		a            = configFile.get<double>("a");               // demi grand-axe (solei-terre hyp MCU)
		nsel_physics = configFile.get<int>("nsel_physics");       //1) one body problem around mass#2 or 2) one body in rotating reference frame of {1,2}
		adapt        = configFile.get<bool>("adapt");             //if 1=true -> adaptive dt, if 0=false -> fixed dt
		tol          = configFile.get<double>("tol");             //tolerance of the adaptive scheme
		sampling     = configFile.get<int>("sampling");     // number of time steps between two writings on file
		nsteps       = configFile.get<int>("nsteps");        // number of time step per period
		mtot = m[0]+m[1];
		alpha = m[1] / mtot;
		beta = m[0] / mtot;

		//TO DO	cheeeeeeck
		v0 = r1 * sqrt(2*G*m[1] * (1/r0 - 1/r1)/(r1*r1 - r0*r0));


		// Ouverture du fichier de sortie
		outputFile = new ofstream(configFile.get<string>("output").c_str());
		outputFile->precision(15);
		//TO DO initialize tFin for nsel_physics=1 and initialize dt for both nsel_physics
		dt=tFin/nsteps;
	}

	~Exercice3() {
		outputFile->close();
		delete outputFile;
	};

	void run() {
		t = 0.;
		initial_condition();
		x = x0;
		last = 0;
		printOut(true);
		valarray<double> y1;
		valarray<double> y2;
		valarray<double> y_tilde;

		if (adapt==false) {
			//TODO fixed dt scheme
			while (t < tFin - 0.5 * dt) {
				x = RK4_do_onestep(x, t, dt);
				t += dt;

				printOut(false);

			}
		}
		else {
			//TODO adaptive case, verify algorithm from polycopie
			nsteps = 0;
			const double f(0.99);

			// TODO verifier cette condition mais je crois ça marche quand même
			while (t < tFin - 0.5 * dt) {
				dt = min(dt, tFin - t);
				++nsteps;

				y1 = RK4_do_onestep(x, t, dt);
				y_tilde = RK4_do_onestep(x, t, dt/2);
				y2 = RK4_do_onestep(y_tilde, t, dt/2);

				d = norm(y1 - y2);

				if (d > tol) {
					do {
						dt = f * dt * pow(tol/d, 1/(4 + 1));			// TODO: RK4 converge ordre 4

						y1 = RK4_do_onestep(x, t, dt);
						y_tilde = RK4_do_onestep(x, t, dt/2);
						y2 = RK4_do_onestep(y_tilde, t, dt/2);

						d = norm(y1 - y2);
					} while (d > tol);
					x = y2;
					t += dt;
				} else {
					x = y2;
					t += dt;
					dt = dt * pow(tol/d, 1/(4 + 1));
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
