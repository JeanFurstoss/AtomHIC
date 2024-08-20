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

int main(int argc, char *argv[])
{
		GaussianMixtureModel GMM;
		GMM.ReadModelParamFromDatabase("AluminaPhase");
		AtomicSystem MySystem("Input.cfg");
		string DescriptorName = "Steinhardt";
		string ftype="element";
		Descriptors MyDescriptors(&MySystem,DescriptorName,ftype);
		GMM.setDescriptors(&MyDescriptors);
		GMM.Classify();
		MySystem.setAux_vec(GMM.getClassificator(),2,"Struct");
		MySystem.printSystem_aux("Output.cfg","Struct");
}
