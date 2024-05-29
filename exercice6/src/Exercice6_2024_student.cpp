#include "ConfigFile.tpp"
#include <chrono>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
typedef vector<complex<double>> vec_cmplx;

// Global constants
constexpr complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i
constexpr double hbar = 1.0;
constexpr double m = 1.0;

enum PotentialType {
    DOUBLE_WELL,
    STEP
};

// Fonction resolvant le systeme d'equations A * solution = rhs
// où A est une matrice tridiagonale
template<class T>
void triangular_solve(
    vector<T> const& diag, vector<T> const& lower, vector<T> const& upper,
    vector<T> const& rhs, vector<T>& solution
) {
    vector<T> new_diag = diag;
    vector<T> new_rhs = rhs;

    // forward elimination
    for (size_t i(1); i < diag.size(); ++i) {
        T pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution.resize(diag.size());

    // solve last equation
    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    // backward substitution
    for (int i = diag.size() - 2; i >= 0; --i) {
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];
    }
}

// Potential V(x)
double V_calculate(double x, double V0, double n_v, double xL, double xR, PotentialType type) {
    if (type == DOUBLE_WELL) {
        return 1.0/2.0 * V0 * (1.0 + cos(2.0 * M_PI * n_v * (x - xL)/(xR - xL)));
    } else if (type == STEP) {
        return ((x - xL)/(xR - xL) > 0.5) * V0;
    }
    return 0.0;
}

// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule entre les points from.dx et to.dx,
//  - E:    calcule son energie moyenne,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.


// Calculate proba between two points given their index
double prob(const vec_cmplx& psi, double dx, size_t from, size_t to) {
    complex<double> cum(0.0, 0.0);

    for (size_t i(from); i < to; i++) {
        cum += (
            conj(psi.at(i)) * psi.at(i)
            + conj(psi.at(i+1)) * psi.at(i+1)
        ) / 2.0;
    }
    cum *= dx;

    return cum.real();
}

// Calculate energy of the system
double E(const vec_cmplx& psi, const vec_cmplx& H_psi, double dx) {
    complex<double> cum(0.0, 0.0);
    for (size_t i(0); i < psi.size() - 1; i++) {
        cum += (
            conj(psi.at(i)) * H_psi.at(i)
            + conj(psi.at(i+1)) * H_psi.at(i+1)
        ) / 2.0;
    }
    cum *= dx;

    return cum.real();
}

// Caculate average position
double xmoy(const vector<double>& x, const vec_cmplx& psi, double dx) {
    complex<double> cum(0.0, 0.0);

    for (size_t i(0); i < x.size() - 1; i++) {
        cum += (
            conj(psi.at(i)) * x.at(i) * psi.at(i)
            + conj(psi.at(i+1)) * x.at(i+1) * psi.at(i+1)
        ) / 2.0;
    }
    cum *= dx;

    return cum.real();
}

// Calculate average squared position
double x2moy(const vector<double>& x, const vec_cmplx& psi, double dx) {
    complex<double> cum(0.0, 0.0);

    for (size_t i(0); i < x.size() - 1; i++) {
        cum += (
            conj(psi.at(i)) * x.at(i) * x.at(i) * psi.at(i)
            + conj(psi.at(i+1)) * x.at(i+1) * x.at(i+1) * psi.at(i+1)
        ) / 2.0;
    }
    cum *= dx;

    return cum.real();
}

// Calculate average momentum
double pmoy(const vec_cmplx& psi, double dx) {
    complex<double> cum(0.0, 0.0);
    size_t end = psi.size() - 1;

    // Calculate first interval
    cum += (
        -conj(psi.at(0)) * complex_i * hbar * ((psi.at(1) - psi.at(0)) / (dx))
        -conj(psi.at(1)) * complex_i * hbar * ((psi.at(2) - psi.at(0)) / (2.0 * dx))
    ) / 2.0;

    // Calculate middle intervals
    for (size_t i(1); i < end - 1; i++) {
        cum += (
            -conj(psi.at(i)) * complex_i * hbar * ((psi.at(i+1) - psi.at(i-1)) / (2.0 * dx))
            -conj(psi.at(i+1)) * complex_i * hbar * ((psi.at(i+2) - psi.at(i)) / (2.0 * dx))
        ) / 2.0;
    }

    // Calculate last interval
    cum += (
        -conj(psi.at(end - 1)) * complex_i * hbar * ((psi.at(end) - psi.at(end - 2)) / (2.0 * dx))
        -conj(psi.at(end)) * complex_i * hbar * ((psi.at(end) - psi.at(end - 1)) / (dx))
    ) / 2.0;
    
    cum *= dx;

    return cum.real();
}

// Calculate average squared momentum
double p2moy(const vec_cmplx& psi, double dx) {
    complex<double> cum(0.0, 0.0);
    size_t end = psi.size() - 1;

    // Calculate first interval
    cum += (
        0.0
        -conj(psi.at(1)) * hbar * hbar * ((psi.at(2) - 2.0 * psi.at(1) + psi.at(0)) / (dx * dx))
    ) / 2.0;

    // Calculate middle intervals
    for (size_t i(1); i < end - 1; i++) {
        cum += (
            -conj(psi.at(i)) * hbar * hbar * ((psi.at(i+1) - 2.0 * psi.at(i) + psi.at(i-1)) / (dx * dx))
            -conj(psi.at(i+1)) * hbar * hbar * ((psi.at(i+2) - 2.0 * psi.at(i+1) + psi.at(i)) / (dx * dx))
        ) / 2.0;
    }
    
    // Calculate last interval
    cum += (
        -conj(psi.at(end - 1)) * hbar * hbar * ((psi.at(end) - 2.0 * psi.at(end - 1) + psi.at(end - 2)) / (dx * dx))
        + 0.0
    ) / 2.0;

    cum *= dx;

    return cum.real();
}

// Function to normalise wave function
void normalize(vec_cmplx& psi, double dx) {
    double norm = prob(psi, dx, 0, psi.size() - 1);

    // Modifies given psi
    for (auto& z: psi) {
        z /= sqrt(norm);
    }
}

// Write the observables of the system on the output file
void write_observables(std::ofstream& fichier_observables, 
                        double t, const vector<double>& x, const vec_cmplx& psi, 
                        const vec_cmplx& H_psi, double dx, size_t Nx0, bool convergence = false) {
    if (!convergence) {
        fichier_observables << t << " "
            << prob(psi, dx, 0, Nx0) << " "
            << prob(psi, dx, Nx0, x.size() - 1) << " "
            << E(psi, H_psi, dx) << " "
            << xmoy(x, psi, dx) << " "
            << x2moy(x, psi, dx) << " "
            << pmoy(psi, dx) << " "
            << p2moy(psi, dx) << endl;
    } else if (convergence) {
        fichier_observables << t << " "
            << prob(psi, dx, 0, Nx0) << " "
            << prob(psi, dx, Nx0, x.size() - 1) << " "
            << E(psi, H_psi, dx) << endl;
    }
}


// Calculate a matrix * vector multiplication with a tri-diagonal representation of the matrix
vec_cmplx diag_matrix_vector(const vec_cmplx& diag, const vec_cmplx& lower, const vec_cmplx& upper, const vec_cmplx& psi) {
    size_t Npoints = diag.size();
    vec_cmplx result(Npoints);

    result.at(0) = diag.at(0)*psi.at(0) + upper.at(0)*psi.at(1);
    result.back() = lower.back()*psi.at(Npoints - 2) + diag.back()*psi.back();

    for (size_t i(1); i < Npoints - 1; ++i) {
        result.at(i) = lower.at(i - 1)*psi.at(i - 1) + diag.at(i)*psi.at(i) + upper.at(i)*psi.at(i + 1);
    }
    return result;
}


// SIMULATION
int main(int argc, char** argv) {
    string inputPath("configuration.in.example");
    if (argc > 1) {
        inputPath = argv[1];
    }

    ConfigFile configFile(inputPath);

    for (int i(2); i < argc; ++i)
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Physical parameters :
    const double tfin = configFile.get<double>("tfin");
    const double xL = configFile.get<double>("xL");
    const double xR = configFile.get<double>("xR");
    const double L = xR - xL;
    const double V0 = configFile.get<double>("V0");
    const double n_v = configFile.get<double>("n_v");
    const double n = configFile.get<int>("n"); // Read mode number as integer, convert to double

    // Numerical parameters :
    double dt = configFile.get<double>("dt");
    const int Nintervals = configFile.get<int>("Nintervals");
    const int Npoints = Nintervals + 1;
    const double dx = L / Nintervals;
    // Is a convergence analysis being done
    bool convergence = configFile.get<bool>("convergence");
    PotentialType pot_type;
    {
        string tmp = configFile.get<string>("potential");
        if (tmp == "well") {
            pot_type = DOUBLE_WELL;
        } else if (tmp == "step") {
            pot_type = STEP;
        } else {
            throw runtime_error("Invalid potential type given");
        }
    }
    
    // Parameters for initialisation
    const double x0 = configFile.get<double>("x0");
    const double k0 = 2 * M_PI * n / L;
    const double sigma0 = configFile.get<double>("sigma_norm") * L;

    // Initialise time and index to check Probability
    double t = 0.0;
    size_t Nx0 = round((Npoints - 1.0) / 2.0);

    // Build the x mesh
    vector<double> x(Npoints);
    for (int i(0); i < Npoints; i++) {
        x.at(i) = xL + i * dx;
    }

    // Initialisation of wave function, equation (4.116) from notes
    vec_cmplx psi(Npoints);
    for (int i(0); i < Npoints; i++) {
        psi.at(i) = exp(complex_i * k0 * x[i]) * exp(-pow((x[i] - x0) / sigma0, 2) / 2.0);
    }

    // Modify boundary values :
    psi.front() *= 0.0;
    psi.back() *= 0.0;

    // Normalisation :
    normalize(psi, dx);

    // TODO check matrices
    // Matrices (d: diagonal, a: lower diagonal, c: upper diagonal) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // Hamiltonian
    vec_cmplx dA(Npoints), aA(Nintervals), cA(Nintervals); // Matrix of left hand side term of equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals), cB(Nintervals); // Matrix of right hand side term of equation (4.100)

    complex<double> a = complex_i * hbar * dt / (4.0 * m * dx * dx); // Complex coefficient a of equation (4.100)

    // Vector of potential at each position
    vector<double> V(Npoints);
    for (size_t i(0); i < V.size(); ++i) {
        V.at(i) = V_calculate(x.at(i), V0, n_v, xL, xR, pot_type);
    }

    // Vector of the values of b of equation (4.100)
    vector<complex<double>> b(Npoints);
    for (size_t i(0); i < b.size(); ++i) {
        b.at(i) = complex_i * dt * V.at(i) / (hbar * 2);
    }

    // Compute elements of A, B and H
    // These matrices are in tridiagonal form, d:diagonal, c and a: upper and lower diagonals
    double const_hamiltonian(hbar * hbar / (2.0 * m * dx * dx));
    for (size_t i(0); i < dH.size(); ++i){
        dH.at(i) = 2.0 * const_hamiltonian + V.at(i);

        dA.at(i) = 1.0 + 2.0*a + b.at(i);
        dB.at(i) = 1.0 - 2.0*a - b.at(i);
    }

    for (size_t i(0); i < aH.size(); ++i){
        aH.at(i) = -const_hamiltonian;
        cH.at(i) = -const_hamiltonian;
        
        aA.at(i) = -a;
        cA.at(i) = -a;

        aB.at(i) = a;
        cB.at(i) = a;
    }

    // Modify A and B matrices to satisfy boundary conditions
    // H is not changed as it does not play a role in the computation of the evolution
    dA.front() = 1.0; dA.back() = 1.0;
    aA.front() = 0.0; aA.back() = 0.0; cA.front() = 0.0; cA.back() = 0.0;

    dB.front() = 1.0; dB.back() = 1.0;
    aB.front() = 0.0; aB.back() = 0.0; cB.front() = 0.0; cB.back() = 0.0;


    // Output files
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel((output + "_pot.out").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i) {
        fichier_potentiel << x[i] << " " << V.at(i) << endl;
    }
    fichier_potentiel.close();

    ofstream fichier_psi((output + "_psi2.out").c_str());
    fichier_psi.precision(15);

    ofstream fichier_observables((output + "_obs.out").c_str());
    fichier_observables.precision(15);

    // t0 writing
    for (int i(0); i < Npoints; ++i) {
        fichier_psi << pow(abs(psi[i]), 2) << " " << real(psi[i]) << " " << imag(psi[i]) << " ";
    }
    fichier_psi << endl;

    // Writing observables at t0
    vec_cmplx H_psi = diag_matrix_vector(dH, aH, cH, psi);
    write_observables(fichier_observables, t, x, psi, H_psi, dx, Nx0, convergence);

    // Time loop:    
    while (t < tfin - 0.5 * dt) {

        // Compute right hand side term
        vec_cmplx psi_tmp(Npoints, 0.);
        psi_tmp = diag_matrix_vector(dB, aB, cB, psi);

        // Solving A * psi = psi_tmp :
        triangular_solve(dA, aA, cA, psi_tmp, psi);
        t += dt;

        // Write psi at t
        if (!convergence) {
            for (int i(0); i < Npoints; ++i) {
                fichier_psi << pow(abs(psi[i]), 2) << " " << real(psi[i]) << " " << imag(psi[i]) << " ";
            }
            fichier_psi << endl;
        }

        // Writing observables
        H_psi = diag_matrix_vector(dH, aH, cH, psi);
        write_observables(fichier_observables, t, x, psi, H_psi, dx, Nx0, convergence);
    }


    fichier_observables.close();
    fichier_psi.close();

}
