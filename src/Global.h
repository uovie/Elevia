// standard C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Eigen matrix algebra library
#include <Eigen/Dense>

// Libint Gaussian integrals library
#include <libint2.hpp>

/* a.u. of length: Bohr radius (bohr) */
// Unit Conversion: angstroms -> bohr
constexpr double angstrom_to_bohr = 1 / 0.52917721092; // 2010 CODATA value

using real_t = libint2::scalar_type;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

namespace Elevia {

    class atom {
    public:
        std::string sym;                  // atomic symbol
        int atom_num;                     // atomic number
        double R[3];                      // atomic position
    };

    class system {
    public:
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
        void open(int argc, char* argv[]);
        //void core();
        void close(int argc, char* argv[]);

        std::ifstream in;
        std::ofstream chk;
        std::ofstream out;
    };
}