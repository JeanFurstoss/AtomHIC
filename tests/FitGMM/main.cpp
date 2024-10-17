#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "Descriptors.h"
#include "GaussianMixtureModel.h"
#include <chrono>

using namespace std;

int main(int argc, char *argv[])
{
	Descriptors MyDescriptors("ForTest");
	GaussianMixtureModel GMM;
	vector<string> Properties;
	Properties.push_back("GMM_TOL_LKH_EM 1e-3");
	Properties.push_back("GMM_MAX_ITER_EM 100");
	Properties.push_back("GMM_NB_MAX_CLUSTER 50");
	Properties.push_back("GMM_ELBOW_FACTOR 0.1");
	Properties.push_back("GMM_NB_INIT 1");
	GMM.ReadProperties(Properties);
	vector<string> KMProp;
	KMProp.push_back("KMEANS_NB_MAX_CLUSTER 100");
	KMProp.push_back("KMEANS_TOL 1e-5");
	KMProp.push_back("KMEANS_MAX_ITER 1000");
	KMProp.push_back("KMEANS_NB_INIT 100");
	GMM.SetKMeansProperties(KMProp);
	GMM.setDescriptors(&MyDescriptors);
	unsigned int nmin= 3;
	unsigned int nmax= 15;
	GMM.fitOptimalGMM(nmin,nmax);
	GMM.Labelling();
	GMM.PrintModelParams("output.dat");
}
