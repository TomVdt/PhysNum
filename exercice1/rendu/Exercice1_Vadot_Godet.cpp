#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // classical math library
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
// Fichier .tpp car inclut fonctions template
#include <numeric>

using namespace std; // opening namespace with c++ basic library

/* La class Engine est le moteur principale de ce code. Il contient
   les methodes de base pour lire / initialiser les inputs,
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{
private:
	// Existing private members of Engine...
	const double pi = 3.1415926535897932384626433832795028841971e0;

	// EngineEuler specific members
	unsigned int maxit; // maximum number of iterations
	double tol;        // iterative method tolerance
	double alpha;        // parameter for choosing which Euler method

	// vVariables definitions
	double tfin;         // final time
	unsigned int nsteps; // number of steps
	double mu; 	       // Magnus parameter
	double mass;         // masse of the ball
	double R;            // radius of the ball
	double omega;        // rotational speed of the ball
	double rho;          // fluid density
	double g;            // gravitational acceleration
	double Ct;           // air drag coefficient 
	double S;            // section surface of the ball (pi * R^2) 

	valarray<double> y0 = std::valarray<double>(0.e0, 4); // Correctly initialized
	valarray<double> y = std::valarray<double>(0.e0, 4); // Correctly initialized

	double t, dt;  // Current time and time step

	unsigned int sampling;  // Number of steps between each recorded diagnostics
	unsigned int last;       // Number of steps since last recorded diagnostic
	ofstream* outputFile;    // Pointer towards output file

	/* Calculating and writing diagnostics in an output file:
	   write: (bool) writing every sample if false
	*/
	void printOut(bool write)
	{
		// y = (x, y, vx, vy)
		double Energy = (
			1.0 / 2.0 * mass * (y[2]*y[2] + y[3]*y[3])
			+ 1.0 / 5.0 * mass * R*R * omega*omega
			+ mass * g * y[1]
		);

		// Writing at every [sampling] steps, except if write is true
		if ((!write && last >= sampling) || (write && last != 1))
		{
			*outputFile << t << " " << y[0] << " " << y[1] << " " \
				<< y[2] << " " << y[3] << " " << Energy << endl; // write output on file
			last = 1;
		}
		else
		{
			last++;
		}
	}

	// Computes value of f from ODE
	void compute_f(valarray<double>& f) const
	{
		const double Frict = 1.0 / 2.0 * Ct * rho * S * sqrt(y[2]*y[2] + y[3]*y[3]) / mass;
		const double R3 = pow(R, 3);
		f[0] = y[2];
		f[1] = y[3];
		f[2] = -mu * R3 * rho * omega * y[3] / mass - Frict * y[2];
		f[3] = mu * R3 * rho * omega * y[2] / mass - g - Frict * y[3];
	}


	// New step method from EngineEuler
	void step()
	{
		unsigned int iteration = 0;
		double error = 999e0;
		valarray<double> f = valarray<double>(0.e0, 4);
		valarray<double> y_old = valarray<double>(y);
		valarray<double> y_control = valarray<double>(y);
		valarray<double> delta_y_EE = valarray<double>(y);

		if (alpha >= 0. && alpha <= 1.0) {
			// Once
			// Using y = y_n, because y *is* y_n initially
			compute_f(f);
			delta_y_EE = alpha * f * dt;

			// Lööps
			while (error > tol && iteration < maxit) {
				// Using y = y_n^k (before updating y)
				compute_f(f);
				y = y_old + delta_y_EE + (1 - alpha) * f * dt;
				
				// Using y = y_n^{k+1} (after updating y)
				compute_f(f);
				y_control = y - y_old - delta_y_EE - (1 - alpha) * f * dt;
				// Taking the norm to get the error
				error = sqrt(y_control[0]*y_control[0] + y_control[1]*y_control[1] + y_control[2]*y_control[2] + y_control[3]*y_control[3]);
				
				// Don't forget to increment
				iteration++;
			}

			// Next position has been reached, increment everything else
			t += dt;
		}
		else
		{
			cerr << "alpha not valid" << endl;
		}

	}

public:
	// Modified constructor
	Engine(ConfigFile configFile)
	{
		// Stocking simulation parameter in class attributes
		tfin = configFile.get<double>("tfin", tfin);	        // reading final time
		nsteps = configFile.get<unsigned int>("nsteps", nsteps); // reading number of steps
		y0[0] = configFile.get<double>("x0", y0[0]);  // initial x position	    
		y0[1] = configFile.get<double>("y0", y0[1]);  // initial y position      
		y0[2] = configFile.get<double>("vx0", y0[2]); // initial speed along x	       
		y0[3] = configFile.get<double>("vy0", y0[3]); // initial speed along y	    
		mass = configFile.get<double>("mass", mass);
		g = configFile.get<double>("g", g);
		omega = configFile.get<double>("omega", omega);
		mu = configFile.get<double>("mu", mu);
		R = configFile.get<double>("R", R);
		rho = configFile.get<double>("rho", rho);
		Ct = configFile.get<double>("Ct", Ct);
		sampling = configFile.get<unsigned int>("sampling", sampling);
		tol = configFile.get<double>("tol", tol);
		maxit = configFile.get<unsigned int>("maxit", maxit);
		alpha = configFile.get<double>("alpha", alpha);
		dt = tfin / nsteps; // calculating the time step

		// Opening output file
		outputFile = new ofstream(configFile.get<string>("output", "output.out").c_str());
		outputFile->precision(15); // Numbers will be written with 15 decimals
	};


	// Virtual destructor
	virtual ~Engine()
	{
		outputFile->close();
		delete outputFile;
	};
	// Entire simulation
	void run()
	{
		S = pi * R * R;
		t = 0.e0; // initialising time
		y = y0;   // initialising y vector
		last = 0; // initialising writing parameter
		printOut(true); // writing initial conditions

		for (unsigned int i(0); i < nsteps; ++i) // loop on steps
		{
			step();  // doing one time step
			printOut(false); // writing current time step
		}
		printOut(true); // writing last time step
	};
};

// main program
int main(int argc, char* argv[])
{
	// Existing main function implementation
	// ...
	string inputPath("configuration.in.example"); // Default input file
	if (argc > 1) // Specific input file from user ("./Exercice2 config_perso.in")
		inputPath = argv[1];

	ConfigFile configFile(inputPath); // Parameters are read and stocked in a "map" of strings.

	for (int i(2); i < argc; ++i) // Additional inputs ("./Exercice2 config_perso.in input_scan=[valeur]")
		configFile.process(argv[i]);

	Engine* engine;

	// Create an instance of Engine instead of EngineEuler
	engine = new Engine(configFile);

	engine->run(); // run simulation

	delete engine; // erasing simulation class
	cout << "Fin de la simulation." << endl;
	return 0;
}


