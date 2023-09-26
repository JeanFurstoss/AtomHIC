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
	}else if( argc == 5 ){
		ReferenceFilename = argv[1];
		ListOfFilename = argv[2];
		istringstream iss_rc(argv[3]);
		iss_rc >> rc;
		StrainFilename = argv[4];
	}else{
		cerr << "Usage: ./AtomicStrain ReferenceFilename Filename_of_list_of_analyzed_filename cutoff_radius(not read in FixedParameters) Filename_of_dump_containing_strain_to_add(optional (TODO) )" << endl;
		cerr << "The program will return an 8 dimension vector per atom containing eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz, shear_invariant, hydrostatic_invariant following the work of Shimizu, Ogata, Li: Mater. Trans. 48 (2007), 2923" << endl;
		cerr << "if an other atomic strain from a dump file has to be add, the dump file must contains these 8 fields in the same ordering and having the name AtomicStrain" << endl;
		return EXIT_FAILURE;
	}
	
	AtomicSystem ReferenceSystem(ReferenceFilename);
	ComputeAuxiliary CA(&ReferenceSystem);

	unsigned int nbAt = ReferenceSystem.getNbAtom();
	double *strains;
	double *buffer_strains;
	if( argc == 5 ){
		cout << "Reading, atomic strain to add from file : " << StrainFilename << endl;
		strains = new double[nbAt*8];
		AtomicSystem AtStrainToAddSystem(StrainFilename);
		unsigned int size_s, ind_s;
		ind_s = AtStrainToAddSystem.getAuxIdAndSize("AtomicStrain",size_s);
		if( size_s != 8 ){
			cerr << "The provided file for atomic strain to add (" << StrainFilename << ") don't have the right number of field for AtomicStrain (8 needed)" << endl;
			return EXIT_FAILURE;
		}else{
			for(unsigned int i=0;i<nbAt;i++){
				for(unsigned int j=0;j<8;j++) strains[i*8+j] = AtStrainToAddSystem.getAux(ind_s)[i*8+j];
			}
		}
	}

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
	string suf="AtomicStrain_";
	string ext = ".cfg";
	if( AnalyzedFilenames.size() != 0 ){
		AtomicSystem AnalyzedSystem(AnalyzedFilenames[0]);
		if( argc == 5 ){
			buffer_strains = CA.Compute_AtomicStrain(AnalyzedSystem,rc);
			for(unsigned int i=0;i<nbAt;i++){
				for(unsigned int j=0;j<8;j++) buffer_strains[i*8+j] += strains[i*8+j];
			}
			AnalyzedSystem.setAux_vec(buffer_strains,8,"AtomicStrain");
		}else AnalyzedSystem.setAux_vec(CA.Compute_AtomicStrain(AnalyzedSystem,rc),8,"AtomicStrain");
		outfilename = suf+AnalyzedFilenames[0].erase(AnalyzedFilenames[0].size()-4)+ext;
		AnalyzedSystem.printSystem_aux(outfilename,"AtomicStrain");
		for(unsigned int i=1;i<AnalyzedFilenames.size();i++){
			AnalyzedSystem.FilenameConstructor(AnalyzedFilenames[i]);
			if( argc == 5 ){
				buffer_strains = CA.Compute_AtomicStrain(AnalyzedSystem,rc);
				for(unsigned int i=0;i<nbAt;i++){
					for(unsigned int j=0;j<8;j++) buffer_strains[i*8+j] += strains[i*8+j];
				}
				AnalyzedSystem.modifyAux_vec(buffer_strains,"AtomicStrain");
			}else{
				AnalyzedSystem.modifyAux_vec(CA.Compute_AtomicStrain(AnalyzedSystem,rc),"AtomicStrain");
			}
			outfilename = suf+AnalyzedFilenames[i].erase(AnalyzedFilenames[i].size()-4)+ext;
			AnalyzedSystem.printSystem_aux(outfilename,"AtomicStrain");
		}
	}
	if( argc == 5 ){
		delete[] strains;
	}
	return 0;
}
