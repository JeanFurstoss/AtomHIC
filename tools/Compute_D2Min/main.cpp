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
	string ReferenceFilename, ListOfFilename, buffer_string, line, StrainFilename;
	double rc;
	if( argc == 4 ){
		ReferenceFilename = argv[1];
		ListOfFilename = argv[2];
		istringstream iss_rc(argv[3]);
		iss_rc >> rc;
	}else{
		cerr << "Usage: ./Compute_D2Min ReferenceFilename Filename_of_list_of_analyzed_filename cutoff_radius(not read in FixedParameters)" << endl;
		cerr << "The program will return the D2Min value according to the formulae in Delbecq et al. 2023" << endl;
		return EXIT_FAILURE;
	}
	
	AtomicSystem ReferenceSystem(ReferenceFilename);
	ComputeAuxiliary CA(&ReferenceSystem);

	unsigned int nbAt = ReferenceSystem.getNbAtom();

	vector<string> AnalyzedFilenames;	
	
	ifstream file(ListOfFilename, ios::in);
	if(file){
                while(file){
        		getline(file,line);
			if( !file ) break;
			istringstream text(line);
			text >> buffer_string;
			AnalyzedFilenames.push_back(buffer_string);
		}
	}else{
		cerr << "Could not open " << ListOfFilename << " file" << endl;
	}
	file.close();

	string outfilename;
	string suf="D2Min_";
	string ext = ".cfg";
	if( AnalyzedFilenames.size() != 0 ){
		AtomicSystem AnalyzedSystem(AnalyzedFilenames[0]);
		AnalyzedSystem.setAux(CA.Compute_D2Min(AnalyzedSystem,rc),"D2Min");
		outfilename = suf+AnalyzedFilenames[0].erase(AnalyzedFilenames[0].size()-4)+ext;
		AnalyzedSystem.printSystem_aux(outfilename,"D2Min");
		for(unsigned int i=1;i<AnalyzedFilenames.size();i++){
			AnalyzedSystem.FilenameConstructor(AnalyzedFilenames[i]);
			AnalyzedSystem.modifyAux_vec(CA.Compute_D2Min(AnalyzedSystem,rc),"D2Min");
			outfilename = suf+AnalyzedFilenames[i].erase(AnalyzedFilenames[i].size()-4)+ext;
			AnalyzedSystem.printSystem_aux(outfilename,"D2Min");
		}
	}
	return 0;
}
