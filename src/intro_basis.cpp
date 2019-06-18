// standard C++ headers
#include <vector>

// Elevia headers
#include "Global.h"

// Libint Gaussian integrals library
#include <libint2.hpp>

using namespace Elevia;
using namespace libint2;

// function declarations
void shell_read(std::ifstream& data_file, std::array<Shell::real_t, 3> R, std::string shell_sym);

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

        std::array<Shell::real_t, 3> R;     // define new type for nuclear position
        R[0] = atoms[a].R[0]; R[1] = atoms[a].R[1]; R[2] = atoms[a].R[2];

        std::string shell_sym;
        basis_data >> shell_sym;
        shell_read(basis_data, R, shell_sym);

        basis_data.close();
    }
    return shells;
}

int shell_symbol_to_l(char sym) {
    switch (sym) {
    case 'S': return 0;
    case 'P': return 1;
    case 'D': return 2;
    case 'F': return 3;
    case 'G': return 4;
    case 'H': return 5;
    case 'I': return 6;
    case 'K': return 7;
    case 'M': return 8;
    case 'N': return 9;
    case 'O': return 10;
    case 'Q': return 11;
    case 'R': return 12;
    case 'T': return 13;
    case 'U': return 14;
    case 'V': return 15;
    case 'W': return 16;
    case 'X': return 17;
    case 'Y': return 18;
    case 'Z': return 19;
    default: throw "invalid angular momentum label";
    }
    return 0;
}

void shell_read(std::ifstream& data_file, std::array<Shell::real_t, 3> R, std::string shell_sym)
{
    int clen;                               // contraction length
    std::string temp_data;
    if (shell_sym.size() == 1) {
        int l = shell_symbol_to_l(shell_sym[0]);
        
        data_file >> clen;
        getline(data_file, temp_data);      // read redundant data

        svector<Shell::real_t> zeta(clen);  // exponents of primitive Gaussians
        svector<Shell::real_t> coef(clen);  // constraction coefficient
        svector<Shell::Contraction> contr;

        for (int i = 0; i < clen; i++) {
            data_file >> zeta[i] >> coef[i];
        }
        contr.push_back({ l, false, coef });
        shells.push_back({ zeta, contr, R });
    }
    else if (shell_sym.size() == 2) {
        int l1 = shell_symbol_to_l(shell_sym[0]);
        int l2 = shell_symbol_to_l(shell_sym[1]);
        
        data_file >> clen;
        getline(data_file, temp_data);

        svector<Shell::real_t> zeta(clen);
        svector<Shell::real_t> coef1(clen), coef2(clen);
        svector<Shell::Contraction> contr1, contr2;

        for (int i = 0; i < clen; i++) {
            data_file >> zeta[i] >> coef1[i] >> coef2[i];
        }
        contr1.push_back({ l1, false, coef1 });
        shells.push_back({ zeta, contr1, R });
        contr2.push_back({ l2, false, coef2 });
        shells.push_back({ zeta, contr2, R });
    }
    else
        throw "unknow shell symbol";
}