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
	if( argc < 5 ){
		cerr << "Usage: ./ClusterAnalysis AtomicInputFilename n_clust_min n_clust_max OutputFilename" << endl;
		cerr << "AtomicInputFilename should be an atomic system" << endl;
		cerr << "The executable returns two auxialiary properties Clust[1] containing the id of cluster which the ions belongs to and Clust[2] containing the MLC" << endl;
		cerr << "For the moment ions cannot be filtered" << endl;
		return EXIT_FAILURE;
	}
	
	cout << "Calculation running using " << omp_get_max_threads() << " threads" << endl;
	
	auto start = high_resolution_clock::now();

	string InputFilename = argv[1];
	unsigned int nclust_min, nclust_max;
	istringstream iss_min(argv[2]);
	iss_min >> nclust_min;
	istringstream iss_max(argv[3]);
	iss_max >> nclust_max;
	string OutputFilename = argv[4];
	AtomicSystem MySystem(InputFilename);
	double *pos = new double[MySystem.getNbAtom()*3];
	for(unsigned int i=0;i<MySystem.getNbAtom();i++){
		pos[i*3] = MySystem.getAtom(i).pos.x;
		pos[i*3+1] = MySystem.getAtom(i).pos.y;
		pos[i*3+2] = MySystem.getAtom(i).pos.z;
	}
	Descriptors MyDes(pos,MySystem.getNbAtom(),3);
	GaussianMixtureModel GMM;
	GMM.setDescriptors(&MyDes);
	GMM.fitOptimalGMM(nclust_min,nclust_max);
	GMM.Classify();
	MySystem.setAux_vec(GMM.getClassificator(),2,"Clust");
	MySystem.printSystem_aux(OutputFilename,"Clust");
	delete[] pos;
}
