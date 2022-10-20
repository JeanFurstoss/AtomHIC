#include "MyStructs.h"
#include "Crystal.h"
#include "OrthorhombicCrystal.h"
#include "AtomHicConfig.h"
#include <iostream>
#include <dirent.h>
#include <vector>
#include <cmath>

using namespace std;

OrthorhombicCrystal::OrthorhombicCrystal(const string& crystalName) : Crystal(crystalName) {
	read_database();
	if( fabs(this->CellParameters[0]-this->CellParameters[1]) < 1e-3 || fabs(this->CellParameters[0]-this->CellParameters[2]) < 1e-3 || fabs(this->CellParameters[2]-this->CellParameters[1]) < 1e-3 ) cout << "The constructed crystal is not orthorhombic, consider changing of crystal class" << endl;
}

