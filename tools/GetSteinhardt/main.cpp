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
	unsigned int nbStruct;
	double x_lo, x_hi, y_lo, y_hi, z_lo, z_hi;
	string InputFilename, OutputFilename, CrystalType, struct_name, filename_ids, SteinhardtStyle;
	if( argc < 3 ){
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename number_of_struct number_of_struct*(AtomIndexListFileName name_of_struct) CrystalType OutputFilename SteinhardtStyle" << endl;
		cerr << "or: AnalyzeBicrystal_ARGS AtomicInputFilename number_of_struct number_of_struct*(x_lo x_hi y_lo y_hi z_lo z_hi name_of_struct) CrystalType OutputFilename SteinhardtStyle" << endl; 
		return EXIT_FAILURE;
	}
	
	InputFilename = argv[1];
	istringstream iss_nbS(argv[2]);
	iss_nbS >> nbStruct;
	if( argc != ( nbStruct*2 + 6 ) && argc != ( nbStruct*7 + 6 ) ){
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename number_of_struct number_of_struct*(AtomIndexListFileName name_of_struct) CrystalType OutputFilename SteinhardtStyle" << endl;
		cerr << "or: AnalyzeBicrystal_ARGS AtomicInputFilename number_of_struct number_of_struct*(x_lo x_hi y_lo y_hi z_lo z_hi name_of_struct) CrystalType OutputFilename SteinhardtStyle" << endl; 
		return EXIT_FAILURE;
	}else if( argc == ( nbStruct*2 + 6 ) ){
		CrystalType = argv[(nbStruct*2)+3];
		Bicrystal MySystem(InputFilename, CrystalType);
		OutputFilename = argv[(nbStruct*2)+4];
		SteinhardtStyle = argv[(nbStruct*2)+5];
		unsigned int buffer_int;
		string line;
		for(unsigned int i=0;i<nbStruct;i++){
			filename_ids = argv[3+(i*2)];
			cout << "Using " << filename_ids << " as index list of analyzed ions" << endl;
			struct_name = argv[4+(i*2)];
			At_index.clear();
			ifstream file(filename_ids, ios::in);
			if(file){
        		        do{
        		                getline(file,line);
					istringstream text(line);
					text >> buffer_int;
					At_index.push_back(buffer_int);
				}while(file);
				MySystem.get_CA()->PrintSteinhardtParam(At_index, struct_name, SteinhardtStyle);
			}else{
				cerr << "Could not open " << filename_ids << " file" << endl;
			}
			file.close();
		}
		MySystem.printSystem_aux(OutputFilename,"Q");
	}else if( argc == ( nbStruct*7 + 6 ) ){
		CrystalType = argv[(nbStruct*7)+3];
		OutputFilename = argv[(nbStruct*7)+4];
		SteinhardtStyle = argv[(nbStruct*7)+5];
		Bicrystal MySystem(InputFilename, CrystalType);
		for(unsigned int i=0;i<nbStruct;i++){
			istringstream iss_x_lo(argv[(i*7)+3]);
			istringstream iss_x_hi(argv[(i*7)+4]);
			istringstream iss_y_lo(argv[(i*7)+5]);
			istringstream iss_y_hi(argv[(i*7)+6]);
			istringstream iss_z_lo(argv[(i*7)+7]);
			istringstream iss_z_hi(argv[(i*7)+8]);
			struct_name = argv[9+(i*7)];
			iss_x_lo >> x_lo;
			iss_x_hi >> x_hi;
			iss_y_lo >> y_lo;
			iss_y_hi >> y_hi;
			iss_z_lo >> z_lo;
			iss_z_hi >> z_hi;
			At_index = MySystem.selectAtomInBox(x_lo,x_hi,y_lo,y_hi,z_lo,z_hi);
			MySystem.get_CA()->PrintSteinhardtParam(At_index, struct_name, SteinhardtStyle);
		}
		MySystem.printSystem_aux(OutputFilename,"Q");
	}else{
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename number_of_struct number_of_struct*(AtomIndexListFileName name_of_struct) CrystalType OutputFilename SteinhardtStyle" << endl;
		cerr << "or: AnalyzeBicrystal_ARGS AtomicInputFilename number_of_struct number_of_struct*(x_lo x_hi y_lo y_hi z_lo z_hi name_of_struct) CrystalType OutputFilename SteinhardtStyle" << endl; 
		return EXIT_FAILURE;
	}
	return 0;
}
