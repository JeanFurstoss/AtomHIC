// AtomHic library files
#include <AtomicSystem.h>
#include <SteinhardtDescriptors.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

void error_msg(){
}

int main(int argc, char *argv[])
{
	if( argc < 3 ){
		cerr << "Usage: ./ComputeDescriptors AtomicInputFilename NameOfDescriptor OutputFilename" << endl;
		cerr << "The descriptor properties will be read from /data/FixedParameters/FixedParameters.dat" << endl;
		cerr << "In addition to the output atomic file, this executable will generate a DescriptorProperties.ath file containing the descriptors properties and which can be used for fitting a ML for instance" << endl;
		return EXIT_FAILURE;
	}
	
	auto start = high_resolution_clock::now();	
	
	string InputFilename, OutputFilename, DescriptorName;
	InputFilename = argv[1];
	DescriptorName = argv[2];
	OutputFilename = argv[3];

	AtomicSystem MySystem(InputFilename);

	if( DescriptorName == "Steinhardt" ){
		// Compute the descriptor
		SteinhardtDescriptors MyDescriptors(&MySystem);
		// Set the auxiliary property
		MySystem.setAux_vec(MyDescriptors.getDescriptors(),MyDescriptors.getDim(),DescriptorName);
		// Print descriptor parameters
		ofstream writefile("DescriptorProperties.ath");
		MyDescriptors.printDescriptorsPropToDatabase(writefile);
		writefile.close();
	}else{ // other developped descriptors could be put here
		cerr << "The descriptor name does not correspond to a descriptor that AtomHIC can compute, aborting" << endl;
		return EXIT_FAILURE;
	}

	MySystem.printSystem_aux(OutputFilename,DescriptorName);
	
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Execution time : " << duration.count() << endl;
	
	return 0;
}
