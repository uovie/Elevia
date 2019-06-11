#include "Global.h"

using namespace Elevia;

void core(void);

int main(int argc, char* argv[])
{
    fio file;
	file.read(argc, argv);
	core();
	file.print(argc, argv);
	return 0;
}