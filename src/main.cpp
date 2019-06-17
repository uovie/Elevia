#include "Global.h"

using namespace Elevia;

void core(void);

fio vie;

int main(int argc, char* argv[])
{
    
	vie.read(argc, argv);
	core();
	vie.print(argc, argv);
	return 0;
}