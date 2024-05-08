#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.tpp"
#include <algorithm>

using namespace std;
constexpr double PI = 3.1415926535897932384626433832795028841971e0;
constexpr double g = 9.81;

enum Boundary {
	FIXED, FREE, EXIT
};
enum Direction {
	LEFT, RIGHT, NONE
};
enum Equation {
	EQ1, EQ2
};
enum Initialisation {
	MODE, DEFAULT
};

void boundary_condition(vector<double>& fnext, vector<double>& fnow,
	double t, double dt,
	const vector<double>& beta2,
	double A, Boundary bc_l, Boundary bc_r, size_t N)
{
	(void)t; (void)dt; (void)A;
	// TODO: Insert boundary conditions
	// TODO: verify
	if (bc_l == FIXED) {
		fnext.at(0) = fnow.at(0);
	}
	if (bc_r == FIXED) {
		fnext.at(N - 1) = fnow.at(N - 1);
	}

	if (bc_l == FREE) {
		fnext.at(0) = fnext.at(1);
	}
	if (bc_r == FREE) {
		fnext.at(N - 1) = fnext.at(N - 2);
	}

	if (bc_l == EXIT) {
		fnext.at(0) = (
			fnow.at(0)
			+ sqrt(beta2.at(0)) * (fnow.at(1) - fnow.at(0))
		);
	}
	if (bc_r == EXIT) {
		fnext.at(N - 1) = (
			fnow.at(N - 1)
			- sqrt(beta2.at(N - 1)) * (fnow.at(N - 1) - fnow.at(N - 2))
		);
	}
}

double finit(double x, double A, double x1, double x2, double xL, double n_init, double xR, Initialisation init) {
	double finit_(0.0);
	if (init == DEFAULT) {
		// TODO verify
		if (x <= x1) {
			finit_ = 0.0;
		} else if (x >= x2) {
			finit_ = 0.0;
		} else {
			finit_ = A / 2.0 * (1.0 - cos(2.0 * PI * (x - x1) / (x2 - x1)));
		}
	} else if (init == MODE) {
		finit_ = A * sin(PI * (2.0 * n_init + 1.0) / (2.0 * (xR - xL)) * (x - xL));
	}
	return finit_;
}

//
// Surcharge de l'operateur pour ecrire les elements d'un tableau
//
template<class T> ostream& operator<<(ostream& o, vector<T> const& v) {
	size_t len(v.size());
	for (size_t i(0); i < (len - 1); ++i) {
		o << v[i] << " ";
	}
	if (len > 0) {
		o << v[len - 1];
	}
	return o;
}

Boundary string_to_boundary(const string& str) {
	if (str == "fixed") {
		return FIXED;
	} else if (str == "free") {
		return FREE;
	} else if (str == "exit") {
		return EXIT;
	} else {
		throw runtime_error("Unknown boundary condition type!");
	}
}

Direction string_to_direction(const string& str) {
	if (str == "left") {
		return LEFT;
	} else if (str == "right") {
		return RIGHT;
	} else if (str == "static") {
		return NONE;
	} else {
		throw runtime_error("Unknown initial direction type!");
	}
}

Initialisation string_to_initialisation(const string& str) {
	if (str == "mode") {
		return MODE;
	} else if (str == "default") {
		return DEFAULT;
	} else {
		throw runtime_error("Unknown initialisation type!");
	}
}

Equation string_to_equation(const string& str) {
	if (str == "Eq1") {
		return EQ1;
	} else if (str == "Eq2") {
		return EQ2;
	} else {
		throw runtime_error("Unknown equation type!");
	}
}

//
// Main
//
int main(int argc, char* argv[]) {
	double dx;
	double dt;
	double t;
	double Nsteps;
	int stride(0);

	string inputPath("configuration.in"); // Fichier d'input par defaut
	if (argc > 1) { // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
		inputPath = argv[1];
	}

	ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

	for (int i(2); i < argc; ++i) {// Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
		configFile.process(argv[i]);
	}

	// Parametres de simulation :
	const double tfin = configFile.get<double>("tfin");
	const int nx = configFile.get<int>("nx"); // nb d'intervalles
	const size_t N = nx + 1;                                // nb de pts de maillage
	double CFL = configFile.get<double>("CFL");
	const double nsteps = configFile.get<double>("nsteps");
	const double A = configFile.get<double>("A");
	const double n_init = configFile.get<double>("n_init");
	const double hL = configFile.get<double>("hL");
	const double hR = configFile.get<double>("hR");
	const double hC = configFile.get<double>("hC");
	const double h00 = configFile.get<double>("h00"); // profondeur, cas uniforme
	const double x1 = configFile.get<double>("x1");
	const double x2 = configFile.get<double>("x2");
	const double xa = configFile.get<double>("xa");
	const double xb = configFile.get<double>("xb");
	const double xc = configFile.get<double>("xc");
	const double xd = configFile.get<double>("xd");
	const double xL = configFile.get<double>("xL");
	const double xR = configFile.get<double>("xR");
	const int n_stride = configFile.get<int>("n_stride");

	// Conditions aux bords:
	const Boundary bc_l = string_to_boundary(configFile.get<string>("cb_gauche"));
	const Boundary bc_r = string_to_boundary(configFile.get<string>("cb_droite"));

	// Type de forme initiale de la vague: selon donnée Eq.(4) ou mode propre
	// (par exemple 'mode' pour mode propre, autrement Eq.(4))
	const Initialisation initialisation = string_to_initialisation(configFile.get<string>("initialisation"));

	// Onde partant initialement vers la gauche ou vers la droite ou statique
	// (par exemple 'left', 'right', 'static')
	const Direction initial_state = string_to_direction(configFile.get<string>("initial_state"));

	// Selecteur pour le cas h0 uniforme:
	const bool v_uniform = configFile.get<bool>("v_uniform");

	// Selecteur pour choisir le pas de temps:
	// true --> dt=tfin/nsteps; t final est exactement tfin
	// false --> dt tel que beta_CFL=1; attention, t final n'est pas exactement tfin
	const bool impose_nsteps = configFile.get<bool>("impose_nsteps");

	vector<double> h0(N);   // profondeur aux points de maillage
	vector<double> vel2(N); // u^2 = g*h_0 aux points de maillage
	vector<double> x(N);    // positions des points de maillage
	vector<double> fpast(N), fnow(N), fnext(N), beta2(N);

	dx = (xR - xL) / (N - 1);
	bool ecrire_f = configFile.get<bool>("ecrire_f"); // Exporter f(x,t) ou non

	// Eq.(1) ou Eq.(2) [ou Eq.(6) (facultatif)]: Eq1, Eq2 ou Eq6
	const Equation equation_type = string_to_equation(configFile.get<string>("equation_type"));

	for (size_t i(0); i < N; ++i) {
		const double x_ = i * dx;
		x.at(i) = xL + x_;

		// TODO initialize the depth h0 and the velocity^2 vel2
		// TODO: verify
		if (v_uniform) {
			h0.at(i) = h00;
			vel2.at(i) = g * h00;
		} else {
			double h0_curr = 0.0;
			if (xL <= x_ && x_ <= xa) {
				h0_curr = hL;
			} else if (xa < x_ && x_ < xb) {
				h0_curr = (1.0 / 2.0) * (hL + hC) + (1.0 / 2.0) * (hL - hC) * cos(PI * (x_ - xa) / (xb - xa));
			} else if (xb <= x_ && x_ <= xc) {
				h0_curr = hC;
			} else if (xc < x_ && x_ < xd) {
				h0_curr = (1.0 / 2.0) * (hR + hC) - (1.0 / 2.0) * (hR - hC) * cos(PI * (x_ - xc) / (xd - xc));
			} else if (xd <= x_ && x_ <= xR) {
				h0_curr = hR;
			}
			h0.at(i) = h0_curr;
			vel2.at(i) = g * h0_curr;
		}
	}

	auto max_vel2 = std::max_element(vel2.begin(), vel2.end());

	// TODO
	// define the dt according to CLF input 
	dt = CFL * dx / sqrt(*max_vel2);
	if (impose_nsteps) {
		// define the dt and CLF when you want to fix nsteps
		dt = tfin / nsteps;
		CFL = sqrt(*max_vel2) * dx / dt;
	}

	// Fichiers de sortie :
	string output = configFile.get<string>("output");

	ofstream fichier_x((output + "_x.out").c_str());
	fichier_x.precision(15);

	ofstream fichier_v((output + "_v.out").c_str());
	fichier_v.precision(15);

	ofstream fichier_f((output + "_f.out").c_str());
	fichier_f.precision(15);

	ofstream fichier_h0((output + "_h0.out").c_str());
	fichier_h0.precision(15);

	// Initialisation des tableaux du schema numerique :

	//TODO initialize f and beta^2
	for (size_t i(0); i < N; ++i)
	{
		// fpast.at(i) = 0.;
		fnow.at(i) = finit(x.at(i), A, x1, x2, xL, n_init, xR, initialisation);
		beta2.at(i) = vel2.at(i) * dt*dt / (dx*dx);

		// TODO initialize beta2, fnow and fpast according to the requests
		// TODO verify
		switch (initial_state) {
			case LEFT: fpast.at(i) = finit(x.at(i) - sqrt(vel2.at(i)) * dt, A, x1, x2, xL, n_init, xR, initialisation); break;
			case RIGHT: fpast.at(i) = finit(x.at(i) + sqrt(vel2.at(i)) * dt, A, x1, x2, xL, n_init, xR, initialisation); break;
			case NONE: fpast.at(i) = fnow.at(i); break;
		}
	}


	cout << "beta2[0] is " << beta2[0] << endl;
	cout << "dt is " << dt << endl;


	// Boucle temporelle :
	for (t = 0.0; t < tfin - 0.5 * dt; t += dt)
	{
		// Ecriture :
		if (stride % n_stride == 0)
		{
			if (ecrire_f) fichier_f << t << " " << fnow << endl;
		}
		++stride;

		// Evolution :
		for (size_t i(1); i < N - 1; ++i)
		{
			// TODO: write the expressions for fnext
			// TODO: verify
			if (equation_type == EQ1) {
				fnext.at(i) = (
					1.0/4.0 * (beta2.at(i+1) - beta2.at(i-1)) * (fnow.at(i+1) - fnow.at(i-1))
					+ beta2.at(i) * (fnow.at(i+1) - 2 * fnow.at(i) + fnow.at(i-1))
					+ 2 * fnow.at(i)
					- fpast.at(i)
				);
			} else if (equation_type == EQ2) {
				fnext.at(i) = (
					2.0 * (1.0 - beta2.at(i)) * fnow.at(i)
					- fpast.at(i)
					+ beta2.at(i) * (fnow.at(i+1) + fnow.at(i-1))
				);
			}
		}

		// TODO add boundary conditions
		// TODO: verify
		boundary_condition(fnext, fnow, t, dt, beta2, A, bc_l, bc_r, N);

		// TODO: faire la mise a jour
		// TODO: verify
		fpast = fnow;
		fnow = fnext;
	}

	if (ecrire_f) fichier_f << t << " " << fnow << endl;
	fichier_x << x << endl;
	fichier_v << vel2 << endl;
	fichier_h0 << h0 << endl;

	fichier_f.close();
	fichier_x.close();
	fichier_v.close();
	fichier_h0.close();

	return 0;
}
