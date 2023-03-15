#ifndef ATOMHIC_MYSTRUCTS_H
#define ATOMHIC_MYSTRUCTS_H

#include <iostream>


struct Position{
	double x;
	double y;
	double z;
};

struct Atom{
	unsigned int type_uint;
	Position pos;
};

#endif
