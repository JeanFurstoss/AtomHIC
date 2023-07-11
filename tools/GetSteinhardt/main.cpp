// AtomHic library files
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"

using namespace std;

void error_msg(){
	cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename number_of_struct number_of_struct*(AtomIndexListFileName name_of_struct) CrystalType OutputFilename SteinhardtStyle AveStyle" << endl;
	cerr << "or: AnalyzeBicrystal_ARGS AtomicInputFilename number_of_struct number_of_struct*(x_lo x_hi y_lo y_hi z_lo z_hi name_of_struct) CrystalType OutputFilename SteinhardtStyle AveStyle" << endl; 
	cerr << "number_of_struct can be 0" << endl;
	cerr << endl;
	cerr << "Possible Steinhardt styles : " << endl;
        cerr << "Mono : Classic Steinhardt parameters using ions ofthe same species" << endl;
        cerr << "Multi : Classic Steinhardt parameters using all ion types" << endl;
        cerr << "Filtered : Steinhardt parameters filtered during their computation using ion type (dim = l*(nbAtomType+1))" << endl;
        cerr << "Possible average styles : " << endl;
	cerr << endl;
        cerr << "none : no averaging" << endl;
        cerr << "Mono : averaging over ions ofthe same species" << endl;
        cerr << "Multi : averaging over all ion types" << endl;
        cerr << "Filtered : filter during the averaging using ion type (dim = l*(nbAtomType+1))" << endl;
}

int main(int argc, char *argv[])
{
	vector<unsigned int> At_index;
	unsigned int nbStruct;
	double x_lo, x_hi, y_lo, y_hi, z_lo, z_hi;
	string InputFilename, OutputFilename, CrystalType, struct_name, filename_ids, SteinhardtStyle, AveStyle;
	if( argc < 3 ){
		error_msg();
		return EXIT_FAILURE;
	}
	
	InputFilename = argv[1];
	istringstream iss_nbS(argv[2]);
	iss_nbS >> nbStruct;
	if( argc != ( nbStruct*2 + 7 ) && argc != ( nbStruct*7 + 7 ) ){
		error_msg();
		return EXIT_FAILURE;
	}else if( argc == ( nbStruct*2 + 7 ) ){
		CrystalType = argv[(nbStruct*2)+3];
		Bicrystal MySystem(InputFilename, CrystalType);
		OutputFilename = argv[(nbStruct*2)+4];
		SteinhardtStyle = argv[(nbStruct*2)+5];
		AveStyle = argv[(nbStruct*2)+6];
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
				MySystem.get_CA()->PrintSteinhardtParam(At_index, struct_name, SteinhardtStyle, AveStyle);
			}else{
				cerr << "Could not open " << filename_ids << " file" << endl;
			}
			file.close();
		}
		double rc = MySystem.get_rcut();
		unsigned int l_sph = MySystem.get_lsph();
		unsigned int size;
		if( SteinhardtStyle == "Filtered" || AveStyle == "Filtered" ) size = (l_sph+1)*(MySystem.getCrystal()->getNbAtomType()+1);
		else size = l_sph+1;
		MySystem.setAux_vec(MySystem.get_CA()->ComputeSteinhardtParameters(rc,l_sph,SteinhardtStyle,AveStyle),size,"Q");
		MySystem.printSystem_aux(OutputFilename,"Q");
	}else if( argc == ( nbStruct*7 + 7 ) ){
		CrystalType = argv[(nbStruct*7)+3];
		OutputFilename = argv[(nbStruct*7)+4];
		SteinhardtStyle = argv[(nbStruct*7)+5];
		AveStyle = argv[(nbStruct*7)+6];
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
			MySystem.get_CA()->PrintSteinhardtParam(At_index, struct_name, SteinhardtStyle, AveStyle);
		}
		double rc = MySystem.get_rcut();
		unsigned int l_sph = MySystem.get_lsph();
		unsigned int size;
		if( SteinhardtStyle == "Filtered" || AveStyle == "Filtered" ) size = (l_sph+1)*MySystem.getCrystal()->getNbAtomType();
		else size = l_sph+1;
		MySystem.setAux_vec(MySystem.get_CA()->ComputeSteinhardtParameters(rc,l_sph,SteinhardtStyle,AveStyle),size,"Q");
		MySystem.printSystem_aux(OutputFilename,"Q");
	}else{
		error_msg();
		return EXIT_FAILURE;
	}
	return 0;
}
