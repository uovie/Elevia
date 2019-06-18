// standard C++ headers
#include <iostream>

// Elevia headers
#include "Global.h"
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
Matrix Elevia::two_body_fock(const std::vector<libint2::Shell>& shells, const Matrix& D);

extern Elevia::fio vie;

void core(void)
{
    using libint2::Shell;
    using libint2::Engine;
    using libint2::Operator;

    /*** ===================== ***/
    /*** introduce a basis set ***/
    /*** ===================== ***/

    std::vector<libint2::Shell> shells = intro_basis(vie.sys.atoms);
    size_t n_bf = 0;    // number of basis functions
    for (auto s = 0; s < shells.size(); ++s)
        n_bf += shells[s].size();

    /*** ===================== ***/
    /*** compute 1-e integrals ***/
    /*** ===================== ***/

    libint2::initialize();

    // compute overlap integrals
    auto S = Elevia::one_body_ints(shells, Operator::overlap);
    std::cout << "\n\tOverlap Integrals:\n";
    std::cout << S << std::endl;

    // compute kinetic-energy integrals
    auto T = 1body_ints(shells, Operator::kinetic);
    cout << "\n\tKinetic-Energy Integrals:\n";
    cout << T << endl;

    // nuclear-attraction integrals Vne
    Matrix Vne = 1body_ints(shells, Operator::nuclear, atoms);
    cout << "\n\tNuclear Attraction Integrals:\n";
    cout << Vne << endl;

    // Core Hamiltonian = T + V
    Matrix H = T + V;
    cout << "\n\tCore Hamiltonian:\n";
    cout << H << endl;

    // T and V no longer needed, free up the memory
    T.resize(0, 0);
    V.resize(0, 0);

    // nuclear repulsion energy Vnn
    double Vnn = 0.0;
    for (auto i = 0; i < vie.sys.atoms.size(); i++)
        for (auto j = i + 1; j < vie.sys.atoms.size(); j++) {
            Eigen::Vector3d R1, R2;
            R1 << vie.sys.atoms[i].R[0], vie.sys.atoms[i].R[1], vie.sys.atoms[i].R[2];
            R2 << vie.sys.atoms[j].R[0], vie.sys.atoms[j].R[1], vie.sys.atoms[j].R[2];
            double Delta_R = (R1 - R2).norm();
            Vnn += (double)vie.sys.atoms[i].atom_num * (double)vie.sys.atoms[j].atom_num / Delta_R;
        }
    std::cout << "Nuclear repulsion energy = " << Vnn << std::endl;
}

