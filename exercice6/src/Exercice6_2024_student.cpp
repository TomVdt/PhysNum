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

// Potentiel V(x) : @TODO write potential
double V_calculate(double x, double V0, double n_v, double xL, double xR) {
    return 1.0/2.0 * V0 * (1.0 + cos(2.0 * M_PI * n_v * (x - xL)/(xR - xL)));
}

// @TODO compute the folliwing quantities
// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule entre les points nL.dx et nR.dx,
//  - E:    calcule son energie moyenne,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.


// Fonction pour normaliser une fonction d'onde :
void normalize(const vector<double>& x, vec_cmplx& psi);

// Les definitions de ces fonctions sont en dessous du main.
// TODO: PROBA IS FUCKED
double prob(const vector<double>& x, const vec_cmplx& psi, double dx, size_t from, size_t to) {
    complex<double> cum(0.0, 0.0);

    for (size_t i(from); i < to - 1; i++) {
        cum += (
            conj(psi.at(i)) * x.at(i) * psi.at(i)
            + conj(psi.at(i+1)) * x.at(i+1) * psi.at(i+1)
        ) / 2.0;
    }
    cum *= dx;

    return cum.real();
}

double E(const vector<double>& x, const vec_cmplx& psi, const vec_cmplx& H_psi, double dx) {
    complex<double> cum(0.0, 0.0);
    for (size_t i(0); i < x.size() - 1; i++) {
        cum += (
            conj(psi.at(i)) * H_psi.at(i)
            + conj(psi.at(i+1)) * H_psi.at(i+1)
        ) / 2.0;
    }
    cum *= dx;

    return cum.real();
}

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

double pmoy(const vector<double>& x, const vec_cmplx& psi, double dx) {
    complex<double> cum(0.0, 0.0);

    for (size_t i(0); i < x.size() - 1; i++) {
        // TODO: simplify expression
        if (i == 0uz) {
            cum += (
                -conj(psi.at(i)) * complex_i * hbar * ((psi.at(i+1) - psi.at(i)) / (dx))
                -conj(psi.at(i+1)) * complex_i * hbar * ((psi.at(i+2) - psi.at(i)) / (2.0 * dx))
        ) / 2.0;
        } else if (i == x.size() - 2) {
            cum += (
                -conj(psi.at(i)) * complex_i * hbar * ((psi.at(i+1) - psi.at(i-1)) / (2.0 * dx))
                -conj(psi.at(i+1)) * complex_i * hbar * ((psi.at(i+1) - psi.at(i)) / (dx))
        ) / 2.0;
        } else {
            cum += (
                -conj(psi.at(i)) * complex_i * hbar * ((psi.at(i+1) - psi.at(i-1)) / (2.0 * dx))
                -conj(psi.at(i+1)) * complex_i * hbar * ((psi.at(i+2) - psi.at(i)) / (2.0 * dx))
        ) / 2.0;
        }
    }
    cum *= dx;

    return cum.real();
}

double p2moy(const vector<double>& x, const vec_cmplx& psi, double dx) {
    complex<double> cum(0.0, 0.0);

    for (size_t i(0); i < x.size() - 1; i++) {
        // TODO: simplify expression
        if (i == 0uz) {
            cum += (
                0.0
                -conj(psi.at(i+1)) * hbar * hbar * ((psi.at(i+2) - 2.0 * psi.at(i+1) + psi.at(i)) / (dx * dx))
            ) / 2.0;
        } else if (i == x.size() - 2) {
            cum += (
                -conj(psi.at(i)) * hbar * hbar * ((psi.at(i+1) - 2.0 * psi.at(i) + psi.at(i-1)) / (dx * dx))
                - 0.0
            ) / 2.0;
        } else {
            cum += (
                -conj(psi.at(i)) * hbar * hbar * ((psi.at(i+1) - 2.0 * psi.at(i) + psi.at(i-1)) / (dx * dx))
                -conj(psi.at(i+1)) * hbar * hbar * ((psi.at(i+2) - 2.0 * psi.at(i+1) + psi.at(i)) / (dx * dx))
            ) / 2.0;
        }
    }
    cum *= dx;

    return cum.real();
}

// @TODO write a function to normalize psi
// TODO: verify
void normalize(const vector<double>& x, vec_cmplx& psi) {
    double norm = 0.0;
    const double h = x.at(1) - x.at(0);

    for (size_t i(0); i < x.size() - 1; i++) {
        norm += (pow(abs(psi.at(i)), 2) + pow(abs(psi.at(i+1)), 2)) / 2.0;
    }
    norm *= h;
    
    // Modifies given psi
    for (auto& z: psi) {
        z /= sqrt(norm);
    }
}

void write_observables(std::ofstream& fichier_observables, double t, const vector<double>& x, const vec_cmplx& psi, const vec_cmplx& H_psi, double dx) {
    fichier_observables << t << " "
        << prob(x, psi, dx, 0, x.size() / 2) << " "
        << prob(x, psi, dx, x.size() / 2, x.size()) << " "
        << E(x, psi, H_psi, dx) << " "
        << xmoy(x, psi, dx) << " "
        << x2moy(x, psi, dx) << " "
        << pmoy(x, psi, dx) << " "
        << p2moy(x, psi, dx) << endl;
}


// Calculate a matrix * vector multiplication with a tri-diagonal representation of the matrix
vec_cmplx diag_matrix_vector(const vec_cmplx& dH, const vec_cmplx& aH, const vec_cmplx& cH, const vec_cmplx& psi) {
    size_t Npoints = dH.size();
    size_t Nintervals = aH.size();
    vec_cmplx H_psi(Npoints);

    H_psi.at(0) = dH.at(0)*psi.at(0) + cH.at(0)*psi.at(1);
    H_psi.at(Npoints - 1) = aH.at(Nintervals - 1)*psi.at(Npoints - 2) + dH.at(Npoints - 1)*psi.at(Npoints - 1);

    for (size_t i(1); i < Npoints - 1; ++i) {
        H_psi.at(i) = aH.at(i - 1)*psi.at(i - 1) + dH.at(i)*psi.at(i) + cH.at(i)*psi.at(i + 1);
    }
    return H_psi;
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

    // Parametres physiques :
    const double tfin = configFile.get<double>("tfin");
    const double xL = configFile.get<double>("xL");
    const double xR = configFile.get<double>("xR");
    const double L = xR - xL;
    const double V0 = configFile.get<double>("V0");
    const double n_v = configFile.get<double>("n_v");
    const double n = configFile.get<int>("n"); // Read mode number as integer, convert to double

    // Parametres numeriques :
    double dt = configFile.get<double>("dt");
    const int Nintervals = configFile.get<int>("Nintervals");
    const int Npoints = Nintervals + 1;
    const double dx = L / Nintervals;
    
    const double x0 = configFile.get<double>("x0");
    const double k0 = 2 * M_PI * n / L;
    const double sigma0 = configFile.get<double>("sigma_norm") * L;

    // initialization time and position to check Probability
    double t = 0;
    unsigned int Nx0 = floor((xR * 0.5 - xL) / (xR - xL) * Npoints); //chosen xR*0.5 since top of potential is at half x domain

    // Not useful, python already times this shit
    // const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    // @TODO build the x mesh
    // TODO: verifier
    for (int i(0); i < Npoints; i++) {
        x.at(i) = xL + i * dx;
    }

    // Initialisation de la fonction d'onde :
    // TODO: initialiser le paquet d'onde, equation (4.116) du cours
    vec_cmplx psi(Npoints);
    for (int i(0); i < Npoints; i++) {
        psi.at(i) = exp(complex_i * k0 * x[i]) * exp(-pow((x[i] - x0) / sigma0, 2) / 2.0);
    }

    // TODO: Modifications des valeurs aux bords :
    psi.front() *= 0.0;
    psi.back() *= 0.0;

    // TODO Normalisation :
    normalize(x, psi);

    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals), cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals), cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

    complex<double> a = complex_i * hbar * dt / (4.0 * m * dx * dx); // Coefficient complexe a de l'equation (4.100)

    // Vector of potential at each position
    vector<double> V(Npoints);
    for (size_t i(0); i < V.size(); ++i) {
        V.at(i) = V_calculate(x.at(i), V0, n_v, xL, xR);
    }

    // Vector of the values of b of equation (4.100)
    vector<complex<double>> b(Npoints);
    for (size_t i(0); i < b.size(); ++i) {
        b.at(i) = complex_i * dt * V.at(i) / (hbar * 2);
    }

    // TODO: calculer les éléments des matrices A, B et H.
    // Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales supérieures et inférieures
    double const_hamiltonian(hbar * hbar / (2.0 * m * dx * dx));
    for (size_t i(0); i < dH.size(); ++i){
        dH.at(i) = 2.0 * const_hamiltonian + V.at(i);

        dA.at(i) = 1.0 + 2.0*a + b.at(i);
        // TODO: CHECK SIGNES, j'ai rajouté le - au b
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

    // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites
    // H is not changed as it does not play a role in the evolution
    dA.front() = 1.0; dA.back() = 1.0;
    aA.front() = 0.0; aA.back() = 0.0; cA.front() = 0.0; cA.back() = 0.0;

    dB.front() = 1.0; dB.back() = 1.0;
    aB.front() = 0.0; aB.back() = 0.0; cB.front() = 0.0; cB.back() = 0.0;

    // Fichiers de sortie :
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

    // Ecriture des observables :
    vec_cmplx H_psi = diag_matrix_vector(dH, aH, cH, psi);
    write_observables(fichier_observables, t, x, psi, H_psi, dx);

    // Boucle temporelle :    
    while (t < tfin) {

        // TODO Calcul du membre de droite :
        // TODO: verify
        vec_cmplx psi_tmp(Npoints, 0.);
        psi_tmp.front() = dB.at(0) * psi.at(0) + cB.at(0) * psi.at(1);
        for (size_t i(1); i < psi.size() - 1; i++) {
            psi_tmp.at(i) = aB.at(i-1) * psi.at(i-1) + dB.at(i) * psi.at(i) + cB.at(i) * psi.at(i+1);
        }
        psi_tmp.back() = aB.back() * psi.at(psi.size() - 2) + dB.back() * psi.back();

        // Resolution de A * psi = psi_tmp :
        triangular_solve(dA, aA, cA, psi_tmp, psi);
        t += dt;

        // t0 writing
        for (int i(0); i < Npoints; ++i) {
            fichier_psi << pow(abs(psi[i]), 2) << " " << real(psi[i]) << " " << imag(psi[i]) << " ";
        }
        fichier_psi << endl;

        // Ecriture des observables :
        H_psi = diag_matrix_vector(dH, aH, cH, psi);
        write_observables(fichier_observables, t, x, psi, H_psi, dx);
    }


    fichier_observables.close();
    fichier_psi.close();

    // NO SPAMMING COUT PLEASE
    // const auto simulationEnd = std::chrono::steady_clock::now();
    // const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    // std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
    //     << " seconds" << std::endl;
}
