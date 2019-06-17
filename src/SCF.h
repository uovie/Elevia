/*
 *  This header file is part of Elevia.
 *
 *  Self-Consistent Field
 *
 */

// standard C++ headers
#include <iostream>

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.hpp>

namespace Elevia {
    class SCF {
        friend std::istream& read();
        friend std::ostream& print();

    public:
        SCF() = default;
        SCF(const std::string& met, const std::string& bas):
            method(met), basis(bas), system(sys) { }


    private:
        std::string method;

    };

    std::istream& read(std::string& method, std::string& basis) {

    }
}
