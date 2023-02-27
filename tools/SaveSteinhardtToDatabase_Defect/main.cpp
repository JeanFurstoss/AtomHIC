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
	if( argc == 5 ){
		InputFilename = argv[1];
		CrystalType = argv[2];
		filename = argv[3];
		filename_ids = argv[4];
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
	}else if( argc == 10 ){
		InputFilename = argv[1];
		CrystalType = argv[2];
		filename = argv[3];
		istringstream iss_x_lo(argv[4]);
		istringstream iss_x_hi(argv[5]);
		istringstream iss_y_lo(argv[6]);
		istringstream iss_y_hi(argv[7]);
		istringstream iss_z_lo(argv[8]);
		istringstream iss_z_hi(argv[9]);
		iss_x_lo >> x_lo;
		iss_x_hi >> x_hi;
		iss_y_lo >> y_lo;
		iss_y_hi >> y_hi;
		iss_z_lo >> z_lo;
		iss_z_hi >> z_hi;
	}else{
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename CrystalType DefectName AtomIndexListFileName" << endl;
		cerr << "or: AnalyzeBicrystal_ARGS AtomicInputFilename CrystalType DefectName x_lo x_hi y_lo y_hi z_lo z_hi" << endl; 
		return EXIT_FAILURE;
	}
	Bicrystal MySystem(InputFilename, CrystalType);
	if( argc == 10 ) At_index = MySystem.selectAtomInBox(x_lo,x_hi,y_lo,y_hi,z_lo,z_hi);
	MySystem.get_CA()->SaveSteinhardtParamToDatabase_Defect(CrystalType,filename,At_index);
	cout << "Steinhardt parameters for the defect " << filename << " have been saved successfully !" << endl;
	cout << "You are invited to create the /data/Steinhardt/" << CrystalType << "/" << filename << "/ directory to store the atomic system and the atom index list or the box bounds you used for creating this reference defect" << endl;
	return 0;
}
