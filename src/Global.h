#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace Elevia {

    struct atom {
        std::string sym;                  // atomic symbol
        int atom_num;                     // atomic number
        double R[3];                      // atomic position
    };

    struct system {
        std::vector<atom> atoms;
        std::string method;               // method
        std::string basis;                // basis set
        int number;                       // the total number of atoms
        int charge;                       // molecular charge
        int spin;                         // spin multiplicity
    };

    class fio {
    public:
        system sys;
        void read(int argc, char* argv[]);
        void print(int argc, char* argv[]);
    private:
        std::ifstream in;
        std::ofstream chk;
        std::ofstream out;
    };
}


