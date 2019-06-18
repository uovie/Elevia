// fio class member function "print()"
#include "Global.h"
#include <iomanip>

using namespace Elevia;

extern Elevia::system sys;

void fio::print(int argc, char *argv[])
{
	std::string argv1(argv[1]);
	std::string filename;
	for (auto c : argv1) {
		if (c != '.') {
			filename += c;
		} else { break;	}
	}

    // check input
	chk.open(filename + ".chk");
    chk << "1. General Infomation\n";
    chk << "Method: " << sys.method << ", Basis Set: " << sys.basis << ", Total Num: "
        << sys.number << ", Charge: " << sys.charge << ", Spin Multi: " << sys.spin << ".\n";
    chk << "\n2. Configuration Infomation\n";
    for (auto it = sys.atoms.begin(); it != sys.atoms.end(); it++) {
        chk << (*it).sym;
        for (int i = 0; i < 3; i++) {
            if ((*it).R[i] >= 0) { chk << "    "; }
            else { chk << "   "; }
            chk << std::fixed << std::setprecision(8) << (*it).R[i];
        }
        chk << std::endl;
    }

    // open output file
	out.open(filename + ".out");

    // close files
	in.close();
	chk.close();
	out.close();
}