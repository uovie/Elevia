// The beauty of simplicity :)

#include "Global.h"
using namespace Elevia;
void core(void);
fio vie;

int main(int argc, char* argv[])
{
	vie.open(argc, argv);
	core();
	vie.close(argc, argv);
	return 0;
}