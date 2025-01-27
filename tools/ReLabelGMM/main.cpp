#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "Descriptors.h"
#include "GaussianMixtureModel.h"
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
	if( argc < 5 ){
		cerr << "Usage: ./FitAndSaveGMM InputDirectory (NameOfDescriptor) n_clust_min n_clust_max NameOfDatabase" << endl;
		cerr << "InputDirectory should contains subdirectories having the names of the labels and containing either: " << endl;
	        cerr << "\t- atomic files (in .cfg format such as printed by AtomHIC) where the name of the auxiliary property to use as descriptor for GMM has to be provided in the NameOfDescriptor argument, in this case it is recommanded to provide the properties of the descriptor in a \"InputDirectory/DescriptorProperties.ath\" file (an example file could be found in /data/ExampleFiles/), the descriptors can be filtered (by type or element for instance) by setting the value of DESCRIPTORS_FILTERING_TYPE in /data/FixedParameters/FixedParameters.dat file" << endl;
		cerr << "\t- or files listing only the descriptors" << endl;
		cerr << "n_clust_min and n_clust_max are the range for searching the optimal GMM number of cluster" << endl;
		cerr << "Once fitted and labelled, the obtained GMM parameters will be stored in a database which could then be used to classify data using the ./GMMClassification executable" << endl;
		return EXIT_FAILURE;
	}
	
	cout << "Calculation running using " << omp_get_max_threads() << " threads" << endl;
	
	auto start = high_resolution_clock::now();

	if( argc == 5 ){	
		string InputDir = argv[1];
		istringstream iss_nmin(argv[2]);
		unsigned int nmin;
		iss_nmin >> nmin;
		istringstream iss_nmax(argv[3]);
		unsigned int nmax;
		iss_nmax >> nmax;
		string DatabaseFilename = argv[4];
		Descriptors MyDescriptors(InputDir);
		GaussianMixtureModel GMM;
		GMM.ReadModelParamFromDatabase(DatabaseFilename);
		GMM.setDescriptors(&MyDescriptors);
		GMM.Labelling();
		GMM.PrintToDatabase(DatabaseFilename);
	}else{
		string InputDir = argv[1];
		string DescriptorName = argv[2];
		istringstream iss_nmin(argv[3]);
		unsigned int nmin;
		iss_nmin >> nmin;
		istringstream iss_nmax(argv[4]);
		unsigned int nmax;
		iss_nmax >> nmax;
		string DatabaseFilename = argv[5];
		Descriptors MyDescriptors(InputDir,DescriptorName);
		GaussianMixtureModel GMM;
		GMM.ReadModelParamFromDatabase(DatabaseFilename);
		GMM.setDescriptors(&MyDescriptors);
		GMM.Labelling();
		GMM.PrintToDatabase(DatabaseFilename);
	}

	auto end = high_resolution_clock::now();	
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Resolution time = " << duration.count() << endl;
}
