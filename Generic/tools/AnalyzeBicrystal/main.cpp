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
	cout << "output density file name : " ;
	cin >> filename;
	MySystem.Print1dDensity("Disorder", filename);
	return 0;
}
