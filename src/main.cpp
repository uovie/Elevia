#include "main.h"

int main(int argc, char* argv[FILENAME_MAX])
{
	read(argc, argv);
	core();
	print(argc, argv);
	return 0;
}