#ifndef GAUSSIANMIXTUREMODEL_H
#define GAUSSIANMIXTUREMODEL_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "MachineLearningModel.h"
#include "KMeans.h"

class ATOMHIC_EXPORT GaussianMixtureModel : public MachineLearningModel {
private:
	unsigned int nbClust;
	double *weights; // weights[k] weight of cluster k
	double *mu; // mu[k*dim+d] d component of the mean of cluster k 
	long double *V; // V[k*dim2+d1*dim+d2] d1,d2 component of variance of cluster k
	long double *V_inv;
	long double *det_V;
	long double LogLikelihood;
	double BIC;
	
	// old variables for EM iterations
	double *weights_old;
	double *mu_old;
	long double *V_old;
	long double *V_inv_old;
	long double *det_V_old;
	
	// saved variables for training with multiple number of clusters
	unsigned int *saved_nbClust;
	double *saved_bic;
	double *saved_weights; // saved_weight[c*nbMaxClusters+k] = weight of cluster k for run number c
	double *saved_mu; // saved_mu[c*nbMaxClusters*dim+k*dim+d]
	long double *saved_V; // saved_V[c*nbMaxClusters*dim2+k*dim2+d1*dim+d2]
	long double *saved_V_inv;
	long double *saved_det_V;
	bool SavedVariablesInitialized = false;

	// buffers for EM	
	double *D_i;
	double *C_di;
	double *buffer_di;
	double *E_d;

	// Parameters to put in FixedParameters
	unsigned int nbMaxClusters = 50;
	double tol_Lkh_EM = 1e-3;
	unsigned int MaxIter_EM = 500;
	double fac_elbow = .1; // reduction factor for considering that its is a real elbow

	bool IsKMeans = false;
	KMeans *MyKM;

public:
	// constructors
	GaussianMixtureModel();
	// methods
	void setDescriptors(Descriptors *D);
	void TrainModel(unsigned int &_nbClust);
	void TrainModel(unsigned int &_nbClust_min, unsigned int &_nbClust_max);
	void InitFromKMeans(unsigned int &_nbClust);
	void SaveVariables(unsigned int &current_nbClust);
	void SetOptimalModel(unsigned int &opt_index);
	void UpdateParams();
	void ComputeLogLikelihood();
	void ComputeBIC();
	void PrintModelParams(std::string filename);
	void EM();
	long double Prob_Cluster(unsigned int index_cluster, unsigned int DescriptorIndex);
	//void ReadModelParamFromDatabase();
	// getters
	double getBIC(){ return BIC; }
	double getLogLikelihood(){ return LogLikelihood; }
	KMeans *getKMeans(){ return MyKM; }
	long double *getCov(){ return V_inv; }
	double *getMu(){ return mu; }
	// destructor
	~GaussianMixtureModel();
};

#endif
