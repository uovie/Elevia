// Self-Consistent Field
#include <iostream>
#include "Mol_Int.h"

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
