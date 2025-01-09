#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "Descriptors.h"
#include "DBScan.h"
#include "AtomicSystem.h"
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
	if( argc != 3 ){
		cerr << "Usage: ./DBScanClustering AtomicInputFilename OutputFilename" << endl;
		cerr << "AtomicInputFilename should be an atomic system" << endl;
		cerr << "The program will return an ovitio output file with the id of cluster (parameters of DBScan are read from FixedParameters.dat)" << endl;
		cerr << "For the moment ions cannot be filtered" << endl;
		return EXIT_FAILURE;
	}
	
	cout << "Calculation running using " << omp_get_max_threads() << " threads" << endl;
	
	auto start = high_resolution_clock::now();

	string InputFilename = argv[1];
	string OutputFilename = argv[2];
	AtomicSystem MySystem(InputFilename);
	unsigned int nbAt = MySystem.getNbAtom();
	double *aux = new double[nbAt];
	unsigned int Struct_ind, size_Struct;
	Struct_ind = MySystem.getAuxIdAndSize("Struct",size_Struct);
	for(unsigned int i=0;i<nbAt;i++){
		if( MySystem.getAux(Struct_ind)[i*size_Struct] == 0 || MySystem.getAux(Struct_ind)[i*size_Struct] == 5 ) aux[i] = 1;
		else aux[i] = 0;
	}
	MySystem.setAux(aux,"Crystal");
	Descriptors MyDes(&MySystem,"Position","Crystal");
	DBScan MyDB;
	MyDB.setDescriptors(&MyDes);
	MyDB.TrainModel(to_string(0.));
	MySystem.setAux_vec(MyDB.getClassificator(),2,"ClusterId");
	MySystem.printSystem_aux(OutputFilename,"ClusterId");

	auto end = high_resolution_clock::now();	
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Resolution time = " << duration.count() << endl;

}
