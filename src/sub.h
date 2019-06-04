#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

void read(int argc, char* argv[FILENAME_MAX]);
void core(void);
void print(int argc, char* argv[FILENAME_MAX]);

extern std::ifstream in;
extern std::ofstream chk;
extern std::ofstream out;

extern std::string method;      // method
extern std::string basis;       // basis set
extern int number;              // the total number of atoms
extern int charge;              // molecular charge
extern int spin;                // spin multiplicity

struct atom {
    std::string label;          // atomic label
    int atom_num;               // atomic number
    double r[3];                // atomic position
};

extern std::vector<atom> molecule;