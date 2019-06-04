// The Output Function: Desired information
#include "sub.h"

void print(int argc, char *argv[FILENAME_MAX])
{
	std::string argv1(argv[1]);
	std::string filename;
	for (auto c : argv1) {
		if (c != '.') {
			filename += c;
		} else { break;	}
	}
	chk.open(filename + ".chk");
    chk << "1. General Infomation\n";
    chk << "Method: " << method << ", Basis Set: " << basis << ", Total Num: "
        << number << ", Charge: " << charge << ", Spin Multi: " << spin << ".\n";
    chk << "\n2. Configuration Infomation\n";
    for (auto it = molecule.begin(); it != molecule.end(); it++) {
        chk << (*it).label;
        for (int i = 0; i < 3; i++) {
            if ((*it).r[i] >= 0) { chk << "    "; }
            else { chk << "   "; }
            chk << std::fixed << std::setprecision(8) << (*it).r[i];
        }
        chk << std::endl;
    }

	out.open(filename + ".out");

	in.close();
	chk.close();
	out.close();
}