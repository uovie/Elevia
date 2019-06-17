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

using namespace Elevia;

extern fio vie;

void core(void)
{
    Eigen::Vector3d r_test(1,2,3);

    double zeta_test = 0.5;
    Eigen::Vector3d n_test(1,0,0);
    Eigen::Vector3d R_test(4,5,6);

    Elevia::PGF pgf_test = Elevia::PGF(r_test, zeta_test, n_test, R_test);
    std::cout << pgf_test.uvalue << '\n' << pgf_test.N << std::endl;

    /*** ===================== ***/
    /*** introduce a basis set ***/
    /*** ===================== ***/



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