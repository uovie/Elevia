// The Input Function: Input the system information.
#include "Global.h"

using namespace Elevia;

Elevia::system sys;

void fio::read(int argc, char *argv[])
{
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << "filename" << std::endl;
		exit(EXIT_FAILURE);
	}

	in.open(argv[1]);

    if (in.fail()) {
        std::cout << "Can not open the file " << argv[1] << '.' << std::endl;
        exit(EXIT_FAILURE);
    }

    in >> sys.method >> sys.basis >> sys.number >> sys.charge >> sys.spin;

    atom atom_data;
    for (int i = 0; i < sys.number; i++) {
        in >> atom_data.sym >> atom_data.R[0] >> atom_data.R[1] >> atom_data.R[2];
        sys.component.push_back(atom_data);
    }

    const std::vector<std::string> element{ "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne",
    "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co",
    "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt",
    "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U" , "Np", "Pu", "Am",
    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };

    for (auto atom_it = sys.component.begin(); atom_it != sys.component.end(); atom_it++) {
        for (auto ele_it = element.begin(); ele_it != element.end(); ele_it++) {
            if ((*atom_it).sym == *ele_it) {
                (*atom_it).atom_num = (int)(ele_it - element.begin()) + 1;
            }
        }
        std::cout << (*atom_it).atom_num << std::endl;
    }
}

