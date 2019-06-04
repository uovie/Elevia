#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

void read(int argc, char *argv[FILENAME_MAX]);
void core(void);
void print(int argc, char *argv[FILENAME_MAX]);

std::ifstream in;
std::ofstream chk;
std::ofstream out;

std::string method;      // method
std::string basis;       // basis set
int number;              // the total number of atoms
int charge;              // molecular charge
int spin;                // spin multiplicity

struct atom {
    std::string label;   // atomic label
    int atom_num;        // atomic number
    double r[3];         // atomic position
};

std::vector<atom> molecule;