#ifndef ATOMHIC_MYSTRUCTS_H
#define ATOMHIC_MYSTRUCTS_H

#include <iostream>


struct Position{
	double x;
	double y;
	double z;
};

struct Atom{
	std::string type;
	unsigned int type_uint;
	Position pos;
	double charge; // charge
	double mass; // mass
};

#endif
