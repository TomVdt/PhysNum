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
    return 1/2 * V0 * (1 + cos(2*M_PI*n_v * (x - xL)/(xR - xL)));
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
vec_cmplx normalize();

// Les definitions de ces fonctions sont en dessous du main.
double prob() {
    return 0.0;
}

double E() {
    return 0.0;
}

double xmoy() {
    return 0.0;
}

double x2moy() {
    return 0.0;
}

double pmoy() {
    return 0.0;
}

double p2moy() {
    return 0.0;
}

//@TODO write a function to normalize psi
vec_cmplx normalize() {
    vec_cmplx psi_norm(1.0, 1.0);
    return psi_norm;
}

void write_observables(std::ofstream& fichier_observables, double t) {
    fichier_observables << t << " "
        << prob() << " "
        << prob() << " " // TODO: huh why twice?
        << E() << " "
        << xmoy() << " "
        << x2moy() << " "
        << pmoy() << " "
        << p2moy() << endl;
}

// SIMULATION
int main(int argc, char** argv) {
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i

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
    const double hbar = 1.0;
    const double m = 1.0;
    const double tfin = configFile.get<double>("tfin");
    const double xL = configFile.get<double>("xL");
    const double xR = configFile.get<double>("xR");
    const double V0 = configFile.get<double>("V0");
    const double n_v = configFile.get<double>("n_v");
    const double n = configFile.get<int>("n"); // Read mode number as integer, convert to double

    // Parametres numeriques :
    double dt = configFile.get<double>("dt");
    const int Nintervals = configFile.get<int>("Nintervals");
    const int Npoints = Nintervals + 1;
    const double dx = (xR - xL) / Nintervals;

    // Not useful, python already times this shit
    // const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    //@TODO build the x mesh

    // Initialisation de la fonction d'onde :
    vec_cmplx psi(Npoints);

    // initialization time and position to check Probability
    double t = 0;
    unsigned int Nx0 = floor((xR * 0.5 - xL) / (xR - xL) * Npoints); //chosen xR*0.5 since top of potential is at half x domain

    const double x0 = configFile.get<double>("x0");
    const double k0 = 2 * M_PI * n / (xR - xL);
    const double sigma0 = configFile.get<double>("sigma_norm") * (xR - xL);
    // TODO: initialiser le paquet d'onde, equation (4.116) du cours

    // TODO: Modifications des valeurs aux bords :

    // TODO Normalisation :
    psi = normalize();

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
        dB.at(i) = 1.0 - 2.0*a + b.at(i);
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
    dA.at(0) = 1.0, dA.at(Npoints - 1) = 1.0;
    aA.at(0) = 0.0, aA.at(Nintervals - 1) = 0.0, cA.at(0) = 0.0, cA.at(Nintervals - 1) = 0.0;

    dB.at(0) = 1.0, dB.at(Npoints - 1) = 1.0;
    aB.at(0) = 0.0, aB.at(Nintervals - 1) = 0.0, cB.at(0) = 0.0, cB.at(Nintervals - 1) = 0.0;

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
    write_observables(fichier_observables, t);

    // Boucle temporelle :    
    while (t < tfin) {

        // TODO Calcul du membre de droite :
        vec_cmplx psi_tmp(Npoints, 0.);


        // Resolution de A * psi = psi_tmp :
        triangular_solve(dA, aA, cA, psi_tmp, psi);
        t += dt;

        // t0 writing
        for (int i(0); i < Npoints; ++i) {
            fichier_psi << pow(abs(psi[i]), 2) << " " << real(psi[i]) << " " << imag(psi[i]) << " ";
        }
        fichier_psi << endl;

        // Ecriture des observables :
        write_observables(fichier_observables, t);
    }


    fichier_observables.close();
    fichier_psi.close();

    // NO SPAMMING COUT PLEASE
    // const auto simulationEnd = std::chrono::steady_clock::now();
    // const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    // std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
    //     << " seconds" << std::endl;
}
