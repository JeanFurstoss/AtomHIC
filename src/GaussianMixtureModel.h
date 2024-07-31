#ifndef GAUSSIANMIXTUREMODEL_H
#define GAUSSIANMIXTUREMODEL_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "MachineLearningModel.h"

class ATOMHIC_EXPORT GaussianMixtureModel : public MachineLearningModel {
private:
	unsigned int nbClust;
	unsigned int dim;
	unsigned int dim2;
	unsigned int nbDat;
	double *weights; // weights[k] weight of cluster k
	double *mu; // mu[k*dim+d] d component of the mean of cluster k 
	long double *V; // V[k*dim2+d1*dim+d2] d1,d2 component of variance of cluster k
	long double *V_inv;
	long double *det_V;
	double *weights_old;
	double *mu_old;
	long double *V_old;
	long double *V_inv_old;
	long double *det_V_old;
	long double LogLikelihood;
	double BIC;

	// buffers for prob
	double *buffer_vec_1_dim;
	double *buffer_vec_2_dim;

	// buffers for EM	
	double *D_i;
	double *C_di;
	double *buffer_di;
	double *E_d;

	// Parameters to put in FixedParameters
	unsigned int nbMaxClusters = 50;
	double tol_Lkh_EM = 1e-5;
	unsigned int MaxIter_EM = 1000;

public:
	// constructors
	GaussianMixtureModel();
	// methods
	void setDescriptors(Descriptors *D);
	void TrainModel();
	void UpdateParams();
	void ComputeLogLikelihood();
	void ComputeBIC();
	void PrintModelParams();
	void EM(unsigned int _nbClust);
	double Prob_Cluster(unsigned int index_cluster, unsigned int DescriptorIndex);
	//void ReadModelParamFromDatabase();
	// destructor
	~GaussianMixtureModel();
	
};

#endif
