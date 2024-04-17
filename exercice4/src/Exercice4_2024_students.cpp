#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ConfigFile.tpp"


using namespace std;

// Résolution d'un système d'équations linéaires par élimination de
// Gauss-Jordan:
template<class T>
vector<T> solve(
    const vector<T>& diag,
    const vector<T>& lower,
    const vector<T>& upper,
    const vector<T>& rhs
) {
    vector<T> solution(diag.size());
    vector<T> new_diag(diag);
    vector<T> new_rhs(rhs);

    for (int i = 1; i < diag.size(); ++i) {
        double pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution[diag.size() - 1] =
      new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    for (int i(diag.size() - 2); i >= 0; --i)
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];

    return solution;
}

//@TODO code kappa
double kappa() {
    return 0.0;
}

//@TODO code Source
double Source() {
    return 0.0;
}

int main(int argc, char* argv[]) {
    // Read the default input
    string inputPath = "configuration.in.example";
    // Optionally override configuration file.
    if (argc > 1) {
        inputPath = argv[1];
    }

    ConfigFile configFile(inputPath);
    // Override settings
    for (int i = 2; i < argc; i++) {
        configFile.process(argv[i]);
    }

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose   = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    const double R      = configFile.get<double>("R");      // radius of cylinder
    const double S0     = configFile.get<double>("S0");     // source parameter
    const double r0     = configFile.get<double>("r0");     // source parameter
    const double kappa0 = configFile.get<double>("kappa0"); // conductivity parameter
    const double kappaR = configFile.get<double>("kappaR"); // conductivity parameter
    const double sigma  = configFile.get<double>("sigma");  // source parameter
    const double alpha  = configFile.get<double>("alpha");  // parameter that allows to switch from equidistant in r to equidistant in r^2
    const double TR     = configFile.get<double>("TR");     // Temperature boundary condition
    const double N      = configFile.get<int>("N");         // Number of finite element intervals
    string fichier_T    = configFile.get<string>("output");
    string fichier_heat = fichier_T+"_heat.out";

    // Create our finite elements
    const int pointCount =  N + 1; // Number of grid points

    // Position of elements @TODO code r[i]
    vector<double> r(pointCount);
    for (int i = 0; i < N + 1; ++i) {
        r[i] = 1.0;
    }

    // Distance between elements @TODO code h[i]
    vector<double> h(pointCount - 1);
    for (int i = 0; i < h.size(); ++i) {
        h[i] = 1.0;
    }

    // Construct the matrices
    vector<double> diagonal(pointCount, 0.0);  // Diagonal
    vector<double> lower(pointCount - 1, 0.0); // Lower diagonal
    vector<double> upper(pointCount - 1, 0.0); // Upper diagonal
    vector<double> rhs(pointCount, 0.0);       // right-hand-side
    
    for (int k = 0; k < N; ++k) { 
        // Matrix  and right-hand-side 
        // @TODO insert contributions from interval k 

        upper[k]        += 0.0;
        lower[k]        += 0.0;
        diagonal[k]     += 1.0; 
        diagonal[k + 1] += 1.0;

        rhs[k]     += 1.0; 
        rhs[k + 1] += 1.0;
    }

    // Boundary conditions @TODO insert boundary conditions


    // Solve the system of equations (do not change the following line!)
    vector<double> temperature = solve(diagonal, lower, upper, rhs);

    // Calculate heat flux
    vector<double> heatFlux(temperature.size() - 1, 0);
    for (int i = 0; i < heatFlux.size(); ++i) {
        //@TODO compute heat flux at mid intervals, use finite element representation
        heatFlux[i]         = 0.0;
    }

    // Export data
    {
        // Temperature
        ofstream ofs(fichier_T);
        ofs.precision(15);

        if (r.size() != temperature.size())
            throw std::runtime_error("error when writing temperature: r and "
                                     "temperature does not have size");
        for (int i = 0; i < temperature.size(); ++i) {
            ofs << r[i] << " " << temperature[i] << endl;
        }
    }

    {
        // Heat flux
        ofstream ofs(fichier_heat);
        ofs.precision(15);

        if (r.size() != (heatFlux.size() + 1))
            throw std::runtime_error("error when writing heatFlux: size of "
                                     "heatFlux should be 1 less than r");
        for (int i = 0; i < heatFlux.size(); ++i) {
            const double midPoint = 0.5 * (r[i + 1] + r[i]);
            ofs << midPoint << " " << heatFlux[i] << endl;
        }
    }

    return 0;
}

