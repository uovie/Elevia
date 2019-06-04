// The Input Function: Input the system information.
#include "sub.h"

void read(int argc, char *argv[FILENAME_MAX])
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

    in >> method >> basis >> number >> charge >> spin;

    atom atom_data;
    for (int i = 0; i < number; i++) {
        in >> atom_data.label >> atom_data.r[0] >> atom_data.r[1] >> atom_data.r[2];
        molecule.push_back(atom_data);
    }

}

