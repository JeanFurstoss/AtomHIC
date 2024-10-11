#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "Descriptors.h"
#include "SteinhardtDescriptors.h"
#include "GaussianMixtureModel.h"
#include "AtomicSystem.h"
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
	if( argc < 4 ){
		cerr << "Usage: ./GMMClassification InputFilename (NameOfDescriptor) NameOfDatabase OutputFilename" << endl;
		cerr << "InputFilename could be an atomic system or a file directly containing the descriptor values, in the latter case NameOfDescriptor argument should not be used" << endl;
		cerr << "The argument NameOfDatabse should be a directory in /data/MachineLearningModels/GaussianMixtureModel/ which could be generated using the ./FitAndSaveGMM executable" << endl;
		cerr << "In the case of an atomic system the obtained output file contains, the computed descriptors if the descriptors are not provided, the index of the labels (Struct[1]) and the maximum likelihood classifier (Struct[2])" << endl;
		cerr << "In the case of a simple file, the descriptors are simply printed again as well as Struct[1] and Struct[2]" << endl;
		GaussianMixtureModel GMM;
		vector<string> basis = GMM.getAvailableDatabases();
		cerr << "Available GMM databases:" << endl;
		for(unsigned int i=0;i<basis.size();i++){
			cerr << basis[i] << endl;
		}
		return EXIT_FAILURE;
	}
	
	cout << "Calculation running using " << omp_get_max_threads() << " threads" << endl;
	
	auto start = high_resolution_clock::now();

	if( argc == 4 ){	
		string InputFilename = argv[1];
		string DatabaseFilename = argv[2];
		string OutputFilename = argv[3];
		GaussianMixtureModel GMM;
		GMM.ReadModelParamFromDatabase(DatabaseFilename);
		AtomicSystem MySystem;
		if( MySystem.FilenameConstructor(InputFilename) ){ // in this case, the provided file is an atomic file without descriptors, thus compute the descriptors needed from the database
			string DescriptorName, buffer_s;
			size_t pos_name;
			for(unsigned int s=0;s<GMM.getDescriptorProperties().size();s++){
				pos_name = GMM.getDescriptorProperties()[s].find("DESCRIPTOR_NAME");
				if( pos_name!=string::npos ){
					istringstream text(GMM.getDescriptorProperties()[s]);
					text >> buffer_s >> DescriptorName;
				}
			}
			if( DescriptorName == "Steinhardt" ){
				SteinhardtDescriptors MyDescriptors(&MySystem,GMM.getDescriptorProperties());
				GMM.setDescriptors(&MyDescriptors);
				GMM.LabelClassification();
				MySystem.setAux_vec(MyDescriptors.getDescriptors(),MyDescriptors.getDim(),"Steinhardt");
				MySystem.setAux_vec(GMM.getClassificator(),2,"Struct");
				MySystem.printSystem_aux(OutputFilename,"Steinhardt Struct");
				return 0;
			}else{ // other developped descriptors could be put here
				cerr << "The descriptor name does not correspond to a descriptor that AtomHIC can compute, aborting" << endl;
				return EXIT_FAILURE;
			}
		}else{
			Descriptors MyDescriptors(InputFilename);
			GMM.setDescriptors(&MyDescriptors);
			GMM.LabelClassification();
			GMM.PrintClassifiedData(OutputFilename);
		}
	}else{
		string InputFilename = argv[1];
		string DescriptorName = argv[2];
		string DatabaseFilename = argv[3];
		string OutputFilename = argv[4];
		GaussianMixtureModel GMM;
		GMM.ReadModelParamFromDatabase(DatabaseFilename);
		string ftype = GMM.getFilteringType();
		cout << endl;
		cout << "WARNING !!! From the database we expect that the provided descriptors (" << DescriptorName << ") are computed using the following properties: " << endl;
		for(unsigned int s=0;s<GMM.getDescriptorProperties().size();s++) cout << GMM.getDescriptorProperties()[s] << endl;
		cout << endl;
		AtomicSystem MySystem(InputFilename);
		Descriptors MyDescriptors(&MySystem,DescriptorName,ftype);
		GMM.setDescriptors(&MyDescriptors);
		GMM.LabelClassification();
		MySystem.setAux_vec(GMM.getClassificator(),2,"Struct");
		MySystem.printSystem_aux(OutputFilename,"Struct");
	}

	auto end = high_resolution_clock::now();	
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Resolution time = " << duration.count() << endl;
}
