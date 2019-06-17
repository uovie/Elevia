// standard C++ headers
#include <vector>

// Elevia headers
#include "Global.h"

// Libint Gaussian integrals library
#include <libint2.hpp>

using namespace Elevia;
using libint2::Shell;

extern Elevia::system sys;
std::vector<Shell> shells;

std::vector<libint2::Shell> intro_basis(const std::vector<atom>& atoms) {

    for (auto a = 0; a < atoms.size(); a++) {
        std::ifstream basis_data;
        basis_data.open("../dat/basis/" + sys.basis + ".g94");
        if (basis_data.fail()) {
            std::cout << "Can not open the basis file " << sys.basis << ".g94." << std::endl;
            exit(EXIT_FAILURE);
        }
        
        std::string temp_data;
        do {
            do {
                getline(basis_data, temp_data);
            } while (temp_data != "****");
            basis_data >> temp_data;
        } while (temp_data != atoms[a].sym);
        getline(basis_data, temp_data);

        std::string shell_sym;
        basis_data >> shell_sym;
        if (shell_sym.size() == 1) {
            if (shell_sym[0] == 'S') {
                int clen = 0;                       // contraction length
                basis_data >> clen;
                getline(basis_data, temp_data);     // read redundant data
                std::vector<double> zeta(clen);     // exponents of primitive Gaussians
                std::vector<double> coef(clen);     // constraction coefficient
                for (int i = 0; i < clen; i++) {
                    basis_data >> zeta[i] >> coef[i];
                }
                shells.push_back(
                    {
                      zeta, {{0, false, coef}},
                      {atoms[a].R[0], atoms[a].R[1], atoms[a].R[2]}
                    }
                );
            }
            else if (shell_sym[0] == 'P') {
                int clen = 0;
                basis_data >> clen;
                getline(basis_data, temp_data);
                std::vector<double> zeta(clen);
                std::vector<double> coef(clen);
                for (int i = 0; i < clen; i++) {
                    basis_data >> zeta[i] >> coef[i];
                }
                shells.push_back(
                    {
                      zeta, {{1, false, coef}},
                      {atoms[a].R[0], atoms[a].R[1], atoms[a].R[2]}
                    }
                );
            }
            else if (shell_sym[0] == 'D') {
                int clen = 0;
                basis_data >> clen;
                getline(basis_data, temp_data);
                std::vector<double> zeta(clen);
                std::vector<double> coef(clen);
                for (int i = 0; i < clen; i++) {
                    basis_data >> zeta[i] >> coef[i];
                }
                shells.push_back(
                    {
                      zeta, {{2, false, coef}},
                      {atoms[a].R[0], atoms[a].R[1], atoms[a].R[2]}
                    }
                );
            }
            else if (shell_sym[0] == 'F') {
                int clen = 0;
                basis_data >> clen;
                getline(basis_data, temp_data);
                std::vector<double> zeta(clen);
                std::vector<double> coef(clen);
                for (int i = 0; i < clen; i++) {
                    basis_data >> zeta[i] >> coef[i];
                }
                shells.push_back(
                    {
                      zeta, {{3, false, coef}},
                      {atoms[a].R[0], atoms[a].R[1], atoms[a].R[2]}
                    }
                );
            }
        }
        else if (shell_sym.size() == 2) {
            if (shell_sym[0] == 'S' && shell_sym[1] == 'P') {
                int clen = 0;
                basis_data >> clen;
                getline(basis_data, temp_data);
                std::vector<double> zeta(clen);
                std::vector<double> coef_s(clen), coef_p(clen);
                for (int i = 0; i < clen; i++) {
                    basis_data >> zeta[i] >> coef_s[i] >> coef_p[i];
                }
                shells.push_back(
                    {
                      zeta, {{0, false, coef_s}},
                      {atoms[a].R[0], atoms[a].R[1], atoms[a].R[2]}
                    }
                );
                shells.push_back(
                    {
                      zeta, {{1, false, coef_p}},
                      {atoms[a].R[0], atoms[a].R[1], atoms[a].R[2]}
                    }
                );
            }
            else
                throw "unknow shell symbol";
        }
    }
    return shells;
}