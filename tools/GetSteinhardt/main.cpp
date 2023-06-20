// AtomHic library files
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"

using namespace std;

int main(int argc, char *argv[])
{
	vector<unsigned int> At_index;
	double x_lo, x_hi, y_lo, y_hi, z_lo, z_hi;
	string InputFilename, CrystalType, filename, filename_ids;
	// check the number of argument
	if( argc == 4 ){
		InputFilename = argv[1];
		filename_ids = argv[2];
		CrystalType = argv[3];
		unsigned int buffer_int;
		ifstream file(filename_ids, ios::in);
		string line;
		if(file){
        	        do{
        	                getline(file,line);
				istringstream text(line);
				text >> buffer_int;
				At_index.push_back(buffer_int);
			}while(file);
		}else{
			cerr << "Could not open " << filename_ids << " file" << endl;
		}
	}else if( argc == 9 ){
		InputFilename = argv[1];
		istringstream iss_x_lo(argv[2]);
		istringstream iss_x_hi(argv[3]);
		istringstream iss_y_lo(argv[4]);
		istringstream iss_y_hi(argv[5]);
		istringstream iss_z_lo(argv[6]);
		istringstream iss_z_hi(argv[7]);
		iss_x_lo >> x_lo;
		iss_x_hi >> x_hi;
		iss_y_lo >> y_lo;
		iss_y_hi >> y_hi;
		iss_z_lo >> z_lo;
		iss_z_hi >> z_hi;
		CrystalType = argv[8];
	}else{
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename AtomIndexListFileName CrystalType " << endl;
		cerr << "or: AnalyzeBicrystal_ARGS AtomicInputFilename x_lo x_hi y_lo y_hi z_lo z_hi CrystalType " << endl; 
		return EXIT_FAILURE;
	}
	Bicrystal MySystem(InputFilename, CrystalType);
	if( argc == 9 ) At_index = MySystem.selectAtomInBox(x_lo,x_hi,y_lo,y_hi,z_lo,z_hi);
	MySystem.get_CA()->PrintSteinhardtParam(At_index);
	return 0;
}
