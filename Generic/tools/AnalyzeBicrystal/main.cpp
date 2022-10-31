// AtomHic library files
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

int main()
{
	string filename;
	cout << "Atomic file name : " ;
	cin >> filename;
	Bicrystal MySystem(filename, "z", "Forsterite");
	cout << "output atomic file name : " ;
	cin >> filename;
	MySystem.printSystem_aux(filename, "Disorder");
	cout << "output GB profile file name : " ;
	cin >> filename;
	MySystem.Print1dDensity(filename, "GBProfile");
	cout << "output estimated GB profile file name : " ;
	cin >> filename;
	MySystem.Print1dDensity(filename, "GBProfile_Gauss");
	return 0;
}
