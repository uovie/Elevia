// fio class member function "close()"
#include "Global.h"
#include <iomanip>

using namespace Elevia;

void fio::close(int argc, char *argv[])
{
    // close files
	in.close();
	chk.close();
	out.close();
}