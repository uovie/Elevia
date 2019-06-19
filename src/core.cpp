// standard C++ headers
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>

// Elevia headers
#include "misc.h"
#include "Mol_Int.h"

// Eigen matrix algebra library
#include <Eigen/Dense>

// Libint Gaussian integrals library
#include <libint2.hpp>

// function declarations
std::vector<libint2::Shell> intro_basis(const std::vector<Elevia::atom>& atoms);
size_t Elevia::nbasis(const std::vector<libint2::Shell>& shells);
size_t Elevia::max_nprim(const std::vector<libint2::Shell>& shells);
int Elevia::max_l(const std::vector<libint2::Shell>& shells);
std::vector<size_t> Elevia::map_shell_to_basis_function(const std::vector<libint2::Shell>& shells);
Matrix Elevia::SAD(const std::vector<Elevia::atom>& atoms);
Matrix Elevia::one_body_ints(const std::vector<libint2::Shell>& shells, libint2::Operator obtype,
    const std::vector<Elevia::atom>& atoms = std::vector<Elevia::atom>());
Matrix Elevia::two_body_fock(const std::vector<libint2::Shell>& shells, const Matrix& P);

extern Elevia::fio vie;

void core(void)
{
    using libint2::Shell;
    using libint2::Engine;
    using libint2::Operator;

    // specify precision
    std::cout << std::fixed << std::setprecision(8);
    vie.out   << std::fixed << std::setprecision(8);

    // system infomation
    std::cout << "\nGeneral Infomation\n"
        << "Method: " << vie.sys.method << ", Basis Set: " << vie.sys.basis << ", Total Num: "
        << vie.sys.number << ", Charge: " << vie.sys.charge << ", Spin Multi: " << vie.sys.spin << '.' << std::endl;
    vie.out   << "\nGeneral Infomation\n"
        << "Method: " << vie.sys.method << ", Basis Set: " << vie.sys.basis << ", Total Num: "
        << vie.sys.number << ", Charge: " << vie.sys.charge << ", Spin Multi: " << vie.sys.spin << '.' << std::endl;

    std::cout << "\nConfiguration R (a.u.)" << std::endl;
    for (auto it = vie.sys.atoms.begin(); it != vie.sys.atoms.end(); it++) {
        std::cout << (*it).sym;
        for (int i = 0; i < 3; i++) {
            if ((*it).R[i] >= 0) { std::cout << "    "; }
            else { std::cout << "   "; }
            std::cout << (*it).R[i];
        }
        std::cout << std::endl;
    }

    // count the number of electrons and
    // calulate the number of doubly occupied orbitals
    int nelectron = 0;
    for (auto i = 0; i < vie.sys.atoms.size(); i++)//
        nelectron += vie.sys.atoms[i].atom_num;//
    const int ndocc = nelectron / 2;

    /*** ============================ ***/
    /*** nuclear repulsion energy Vnn ***/
    /*** ============================ ***/

    double Vnn = 0.0;
    for (auto i = 0; i < vie.sys.atoms.size(); i++) {
        for (auto j = i + 1; j < vie.sys.atoms.size(); j++) {
            double Delta_Rx = vie.sys.atoms[i].R[0] - vie.sys.atoms[j].R[0];
            double Delta_Ry = vie.sys.atoms[i].R[1] - vie.sys.atoms[j].R[1];
            double Delta_Rz = vie.sys.atoms[i].R[2] - vie.sys.atoms[j].R[2];
            double Square_Delta_R = Delta_Rx * Delta_Rx + Delta_Ry * Delta_Ry + Delta_Rz * Delta_Rz;
            double Delta_R = sqrt(Square_Delta_R);
            Vnn += (double)vie.sys.atoms[i].atom_num * (double)vie.sys.atoms[j].atom_num / Delta_R;
        }
    }
    std::cout << "\nNuclear repulsion energy:\n"
        << "V_{nn} = " << Vnn << std::endl;
    vie.out   << "\nNuclear repulsion energy:\n"
        << "V_{nn} = " << Vnn << std::endl;

    /*** ===================== ***/
    /*** introduce a basis set ***/
    /*** ===================== ***/
    
    std::vector<libint2::Shell> shells = intro_basis(vie.sys.atoms);//
    size_t nao = 0;    // number of atomic orbitals (basis functions)
    for (auto s = 0; s < shells.size(); ++s)
        nao += shells[s].size();
    std::cout << "\nThe number of atomic orbitals is " << nao << '.' << std::endl;
    vie.out   << "\nThe number of atomic orbitals is " << nao << '.' << std::endl;

    /*** ===================== ***/
    /*** compute 1-e integrals ***/
    /*** ===================== ***/

    libint2::initialize();

    // overlap integrals
    auto S = Elevia::one_body_ints(shells, Operator::overlap);
    std::cout << "\nOverlap Integrals:\n"
        << "S =\n" << S << std::endl;
    vie.out   << "\nOverlap Integrals:\n"
        << "S =\n" << S << std::endl;

    // kinetic-energy integrals Te
    auto Te = Elevia::one_body_ints(shells, Operator::kinetic);
    std::cout << "\nKinetic-Energy Integrals:\n"
        << "T_{e} =\n" << Te << std::endl;
    vie.out   << "\nKinetic-Energy Integrals:\n"
        << "T_{e} =\n" << Te << std::endl;

    // nuclear-attraction integrals Vne
    Matrix Vne = Elevia::one_body_ints(shells, Operator::nuclear, vie.sys.atoms);//
    std::cout << "\nNuclear Attraction Integrals:\n"
        << "V_{ne} =\n" << Vne << std::endl;
    vie.out   << "\nNuclear Attraction Integrals:\n"
        << "V_{ne} =\n" << Vne << std::endl;

    // Core Hamiltonian Hcore = Te + Vne
    Matrix Hcore = Te + Vne;
    std::cout << "\nCore Hamiltonian:\n"
        << "H_{core} =\n" << Hcore << std::endl;
    vie.out   << "\nCore Hamiltonian:\n"
        << "H_{core} =\n" << Hcore << std::endl;

    // Te and Vne no longer needed, free up the memory
    Te.resize(0, 0);
    Vne.resize(0, 0);

    /*** ==================================== ***/
    /*** assum an initial bond-order matrix P ***/
    /*** ==================================== ***/

    Matrix P;

    // use Superposition-Of-Atomic-Densities (SOAD) guess for minimal basis set
    // use core Hamiltonian guass for other basis sets
    bool use_hcore_guess = true;
    std::string basis_sto3g = "STO-3G";
    if (vie.sys.basis == basis_sto3g) {
        use_hcore_guess = false;
    }
                                         
    if (use_hcore_guess) {
        // solve H C = e S C
        /*
         * Account for the following eigenvalues_module in Eigen
         *
         * \class GeneralizedSelfAdjointEigenSolver
         *
         * \brief Computes eigenvalues and eigenvectors of the generalized selfadjoint eigen problem
         *
         * This class solves the generalized eigenvalue problem
         * \f$ Av = \lambda Bv \f$. In this case, the matrix \f$ A \f$ should be
         * selfadjoint and the matrix \f$ B \f$ should be positive definite.
         */
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(Hcore, S);
        auto eps = gen_eig_solver.eigenvalues();
        auto C = gen_eig_solver.eigenvectors();
        std::cout << "\nInitial Coefficient Matrix:\n"
            << "C =\n" << C << std::endl;
        vie.out   << "\nInitial Coefficient Matrix:\n"
            << "C =\n" << C << std::endl;

        // P = C(occ) . C(occ)T
        auto C_occ = C.leftCols(ndocc);
        P = C_occ * C_occ.transpose();
    }
    else {  // SOAD as the guess density, assumes STO-nG basis
        P = Elevia::SAD(vie.sys.atoms);
    }
    
    std::cout << "\nInitial Density Matrix:\n"
        << "P =\n" << P << std::endl;
    vie.out   << "\nInitial Density Matrix:\n"
        << "P =\n" << P << std::endl;

    /*** ========================= ***/
    /*** main iterative loop (SCF) ***/
    /*** ========================= ***/

    const auto maxiter = 100;
    const real_t conv = 1e-12;
    auto iter = 0;
    real_t RMSD = 0.0;
    real_t Delta_E = 0.0;
    real_t Eelec = 0.0;
    do {
        const auto tstart = std::chrono::high_resolution_clock::now();
        ++iter;

        // Save a copy of the energy and the density
        auto Ehf_last = Eelec;
        auto P_last = P;

        // build a new Fock matrix
        auto F = Hcore;
        //F += compute_2body_fock_simple(shells, P);
        F += Elevia::two_body_fock(shells, P);

        if (iter == 1) {
            std::cout << "\nFock Matrix:\n"
                << "F =\n" << F << std::endl;
            vie.out   << "\nFock Matrix:\n"
                << "F =\n" << F << std::endl;
        }

        // solve F C = e S C
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
        auto eps = gen_eig_solver.eigenvalues();
        auto C = gen_eig_solver.eigenvectors();

        // compute density, D = C(occ) . C(occ)T
        auto C_occ = C.leftCols(ndocc);
        P = C_occ * C_occ.transpose();

        // compute HF energy
        Eelec = 0.0;
        for (auto i = 0; i < nao; i++)
            for (auto j = 0; j < nao; j++)
                Eelec += P(i, j) * (Hcore(i, j) + F(i, j));

        // compute difference with last iteration
        Delta_E = Eelec - Ehf_last;
        RMSD = (P - P_last).norm();

        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        if (iter == 1) {
            std::cout << "\n\nSCF Procedure:"
                "\nIter       E_{elec}            E_{tot}             \\Delta E              RMSD        Time(s)\n";
            vie.out   << "\n\nSCF Procedure:"
                "\nIter       E_{elec}            E_{tot}             \\Delta E              RMSD        Time(s)\n";
        }
        std::cout << std::setfill('0') << std::setw(2) << iter << std::setprecision(12)
            << std::setfill(' ') << std::setw(20) << Eelec << std::setw(20) << Eelec + Vnn
            << std::setw(20) << Delta_E << std::setw(20) << RMSD
            << std::setprecision(5) << std::setw(10) << time_elapsed.count() << std::endl;
        vie.out   << std::setfill('0') << std::setw(2) << iter << std::setprecision(12)
            << std::setfill(' ') << std::setw(20) << Eelec << std::setw(20) << Eelec + Vnn
            << std::setw(20) << Delta_E << std::setw(20) << RMSD
            << std::setprecision(5) << std::setw(10) << time_elapsed.count() << std::endl;

    } while (((fabs(Delta_E) > conv) || (fabs(RMSD) > conv)) && (iter < maxiter));

    libint2::finalize();

    std::cout << std::setprecision(8);
    std::cout << "\n**Hartree-Fock energy: HF = " << Eelec + Vnn << std::endl;
    std::cout << "\nNormal termination. Congratulations!" << std::endl;
    vie.out   << std::setprecision(8);
    vie.out   << "\n**Hartree-Fock energy: HF = " << Eelec + Vnn << std::endl;
    vie.out   << "\nNormal termination. Congratulations!" << std::endl;
}

