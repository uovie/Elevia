// standard C++ headers
#include <iostream>
#include <iomanip>

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
Matrix Elevia::one_body_ints(const std::vector<libint2::Shell>& shells, libint2::Operator obtype,
    const std::vector<Elevia::atom>& atoms = std::vector<Elevia::atom>());
Matrix Elevia::two_body_fock(const std::vector<libint2::Shell>& shells, const Matrix& P);

extern Elevia::fio vie;

void core(void)
{
    using libint2::Shell;
    using libint2::Engine;
    using libint2::Operator;

    std::string indent(4, ' ');

    std::vector<Elevia::atom> temp_atoms;
    for (auto i = 0; i < vie.sys.atoms.size(); i++) {
        temp_atoms.push_back({ vie.sys.atoms[i].sym, vie.sys.atoms[i].atom_num,
            {vie.sys.atoms[i].R[0], vie.sys.atoms[i].R[1], vie.sys.atoms[i].R[2]} });
    }

    // count the number of electrons and
    // calulate the number of doubly occupied orbitals
    int nelectron = 0;
    for (auto i = 0; i < temp_atoms.size(); i++)
        nelectron += temp_atoms[i].atom_num;
    const int ndocc = nelectron / 2;

    /*** ============================ ***/
    /*** nuclear repulsion energy Vnn ***/
    /*** ============================ ***/

    double Vnn = 0.0;
    for (auto i = 0; i < vie.sys.atoms.size(); i++) {
        for (auto j = i + 1; j < vie.sys.atoms.size(); j++) {
            Eigen::Vector3d R1, R2;
            R1 << vie.sys.atoms[i].R[0], vie.sys.atoms[i].R[1], vie.sys.atoms[i].R[2];
            R2 << vie.sys.atoms[j].R[0], vie.sys.atoms[j].R[1], vie.sys.atoms[j].R[2];
            double Delta_R = (R1 - R2).norm();
            Vnn += (double)vie.sys.atoms[i].atom_num * (double)vie.sys.atoms[j].atom_num / Delta_R;
        }
    }
    vie.out << "\nNuclear repulsion energy:\n" << indent
        << "V_{nn} = " << Vnn << std::endl;
    std::cout << "\nNuclear repulsion energy:\n" << indent
        << "V_{nn} = " << Vnn << std::endl;

    /*** ===================== ***/
    /*** introduce a basis set ***/
    /*** ===================== ***/
    
    std::vector<libint2::Shell> shells = intro_basis(temp_atoms);
    size_t nao = 0;    // number of atomic orbitals (basis functions)
    for (auto s = 0; s < shells.size(); ++s)
        nao += shells[s].size();

    /*** ===================== ***/
    /*** compute 1-e integrals ***/
    /*** ===================== ***/

    libint2::initialize();

    // overlap integrals
    auto S = Elevia::one_body_ints(shells, Operator::overlap);
    vie.out << "\nOverlap Integrals:\n" << indent
        << "S = " << S << std::endl;
    std::cout << "\nOverlap Integrals:\n" << indent
        << "S = " << S << std::endl;

    // kinetic-energy integrals Te
    auto Te = Elevia::one_body_ints(shells, Operator::kinetic);
    vie.out << "\nKinetic-Energy Integrals:\n" << indent
        << "T_{e} = " << Te << std::endl;
    std::cout << "\nKinetic-Energy Integrals:\n" << indent
        << "T_{e} = " << Te << std::endl;

    // nuclear-attraction integrals Vne
    Matrix Vne = Elevia::one_body_ints(shells, Operator::nuclear, temp_atoms);
    vie.out << "\nNuclear Attraction Integrals:\n" << indent
        << "V_{ne} = " << Vne << std::endl;
    std::cout << "\nNuclear Attraction Integrals:\n" << indent
        << "V_{ne} = " << Vne << std::endl;

    // Core Hamiltonian Hcore = Te + Vne
    Matrix Hcore = Te + Vne;
    vie.out << "\nCore Hamiltonian:\n" << indent
        << "H_{core} = " << Hcore << std::endl;
    std::cout << "\nCore Hamiltonian:\n" << indent
        << "H_{core} = " << Hcore << std::endl;

    // Te and Vne no longer needed, free up the memory
    Te.resize(0, 0);
    Vne.resize(0, 0);

    /*** ==================================== ***/
    /*** assum an initial bond-order matrix P ***/
    /*** ==================================== ***/

    Matrix P;

    // In the zeroth iteration, we put P = 0,
    // namely solve H C = e S C.

    /* Account for the following eigenvalues_module in Eigen
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
    vie.out << "\nInitial C Matrix:\n" << indent
        << "C = " << Hcore << std::endl;
    std::cout << "\nInitial C Matrix:\n" << indent
        << "C = " << Hcore << std::endl;

    // P = C(occ) . C(occ)T
    auto C_occ = C.leftCols(ndocc);
    P = C_occ * C_occ.transpose();

    /*** ========================= ***/
    /*** main iterative loop (SCF) ***/
    /*** ========================= ***/

    const auto maxiter = 100;
    const real_t conv = 1e-8;
    auto iter = 0;
    real_t RMSD = 0.0;
    real_t Ediff = 0.0;
    real_t Ehf = 0.0;
    do {
        const auto tstart = std::chrono::high_resolution_clock::now();
        ++iter;

        // Save a copy of the energy and the density
        auto Ehf_last = Ehf;
        auto P_last = P;

        // build a new Fock matrix
        auto F = Hcore;
        //F += compute_2body_fock_simple(shells, P);
        F += Elevia::two_body_fock(shells, P);

        if (iter == 1) {
            vie.out << "\nFock Matrix:\n" << indent
                << "F = " << F << std::endl;
            std::cout << "\nFock Matrix:\n" << indent
                << "F = " << F << std::endl;
        }

        // solve F C = e S C
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
        auto eps = gen_eig_solver.eigenvalues();
        auto C = gen_eig_solver.eigenvectors();

        // compute density, D = C(occ) . C(occ)T
        auto C_occ = C.leftCols(ndocc);
        P = C_occ * C_occ.transpose();

        // compute HF energy
        Ehf = 0.0;
        for (auto i = 0; i < nao; i++)
            for (auto j = 0; j < nao; j++)
                Ehf += P(i, j) * (Hcore(i, j) + F(i, j));

        // compute difference with last iteration
        Ediff = Ehf - Ehf_last;
        RMSD = (P - P_last).norm();

        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        if (iter == 1)
            std::cout <<
            "\n\n Iter        E(elec)              E(tot)               Delta(E)             RMS(D)         Time(s)\n";
        printf(" %02d %20.12f %20.12f %20.12f %20.12f %10.5lf\n", iter, Ehf, Ehf + Vnn,
            Ediff, RMSD, time_elapsed.count());

    } while (((fabs(Ediff) > conv) || (fabs(RMSD) > conv)) && (iter < maxiter));

    libint2::finalize();

    std::cout.precision(8);
    std::cout << std::fixed << "Hartree-Fock energy = " << Ehf + Vnn << std::endl;
}

