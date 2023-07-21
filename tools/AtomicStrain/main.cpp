// AtomHic library files
#include <AtomicSystem.h>
#include <ComputeAuxiliary.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"

using namespace std;

int main(int argc, char *argv[])
{
	string ReferenceFilename, ListOfFilename, buffer_string, line;
	double rc;
	if( argc == 4 ){
		ReferenceFilename = argv[1];
		ListOfFilename = argv[2];
		istringstream iss_rc(argv[3]);
		iss_rc >> rc;
	}else{
		cerr << "Usage: ./BondOriParam AtomicInputFilename CrystalType(The crystal type has to be defined in /data/Crystal/) outputFilename" << endl;
		cerr << "or: ./BondOriParam AtomicInputFilename outputFilename" << endl;
		cerr << "In the latter case, the crystal will be considered as a monosite crystal and the value of spherical harmonic degree and cutoff radius will be read in (/data/FixedParameters/FixedParameters.dat)" << endl; 
		cerr << "The first case should be used for multisite crystals for which the Steinhardt parameters have been stored in the AtomHic database using the SaveSteinhardtToDatabase_PerfectCrystal executable" << endl;
		return EXIT_FAILURE;
	}
	
	AtomicSystem ReferenceSystem(ReferenceFilename);
	ComputeAuxiliary CA(&ReferenceSystem);

	vector<string> AnalyzedFilenames;	
	
	ifstream file(ListOfFilename, ios::in);
	if(file){
                do{
        		getline(file,line);
			istringstream text(line);
			text >> buffer_string;
			AnalyzedFilenames.push_back(buffer_string);
		}while(file);
	}else{
		cerr << "Could not open " << ListOfFilename << " file" << endl;
	}
	file.close();

	string outfilename;
	string suf="AtomicStrain_";
	for(unsigned int i=0;i<AnalyzedFilenames.size();i++){
		AtomicSystem AnalyzedSystem(AnalyzedFilenames[i]);
		AnalyzedSystem.setAux_vec(CA.Compute_AtomicStrain(AnalyzedSystem,rc),9,"AtomicStrain");
		outfilename = suf+AnalyzedFilenames[i];
		AnalyzedSystem.printSystem_aux(outfilename,"AtomicStrain");
		cout << AnalyzedFilenames[i] << " done !" << endl;
		AnalyzedSystem.~AtomicSystem();
	}
	return 0;
}
