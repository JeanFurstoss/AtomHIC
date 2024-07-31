#include "GaussianMixtureModel.h"
#include "AtomHicConfig.h"
#include <cmath>

using namespace std; 

GaussianMixtureModel::GaussianMixtureModel(){
	this->name = "GaussianMixtureModel"; // TODO the user could affect name to the model or it could be read from database
	buffer_vec_1_dim = new double[dim]; //TODO change place
	buffer_vec_2_dim = new double[dim]; //TODO change place
}

void GaussianMixtureModel::setDescriptors(Descriptors *D){
	MachineLearningModel::setDescriptors(D);
	dim = _MyDescriptors->getDim();
	dim2 = dim*dim;
	nbDat = _MyDescriptors->getNbDat();
}

void GaussianMixtureModel::TrainModel(){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, training aborted" << endl;
		exit(EXIT_FAILURE);
	}
	weights = new double[nbMaxClusters];
	mu = new double[dim*nbMaxClusters];
	V = new long double[dim*dim*nbMaxClusters];
	V_inv = new long double[dim*dim*nbMaxClusters];
	det_V = new long double[nbMaxClusters];
	
	D_i = new double[nbDat];
	C_di = new double[nbMaxClusters*nbDat];
	buffer_di = new double[nbDat*nbDat];
	E_d = new double[nbMaxClusters];
	
	weights_old = new double[nbMaxClusters];
	mu_old = new double[dim*nbMaxClusters];
	V_old = new long double[dim*dim*nbMaxClusters];
	V_inv_old = new long double[dim*dim*nbMaxClusters];
	det_V_old = new long double[nbMaxClusters];

	// testing	
	weights[0] = 0.6;
	mu[0] = 0.1;
        mu[1] = 2.9;
        mu[2] = 4.1;
        V[0] = 0.5;
        V[4] = 0.12;
        V[8] = 0.22;
        V[1] = 0.083;
        V[2] = 0.014;
        V[5] = 0.092;
        V[3] = 0.083;
        V[6] = 0.014;
        V[7] = 0.092;
	unsigned int index = 0;
	MT->invMat_LU(V,V_inv,dim,index,det_V[index]);

	weights[1] = 0.4;
	mu[3+0] = 1.1;
        mu[3+1] = 1.9;
        mu[3+2] = 3.1;
	V[9+0] = 0.1;
	V[9+4] = 0.22;
	V[9+8] = 0.32;
	V[9+1] = 0.023;
	V[9+2] = 0.034;
	V[9+5] = 0.042;
	V[9+3] = 0.023;
	V[9+6] = 0.034;
	V[9+7] = 0.042;
	index = 1;
	MT->invMat_LU(V,V_inv,dim,index,det_V[index]);

	nbClust = 2;

	UpdateParams();
	ComputeLogLikelihood();
	long double Lkh_old = LogLikelihood;
	double eps;
	unsigned int iter = 0;
	do{
		EM(2);
		UpdateParams();
		ComputeLogLikelihood();
		eps = LogLikelihood - Lkh_old;
		Lkh_old = LogLikelihood;
		iter++;
	}while( eps > tol_Lkh_EM && iter < MaxIter_EM );

	if( eps < tol_Lkh_EM ){
		cout << "The EM algorithm converged in " << iter << " iteration" << endl;
		cout << "Final log likelihood " << LogLikelihood << endl;
		PrintModelParams();
	}
}

void GaussianMixtureModel::UpdateParams(){
	for(unsigned int d=0;d<nbClust;d++){
		weights_old[d] = weights[d];
		det_V_old[d] = det_V[d];
		for(unsigned int x1=0;x1<dim;x1++){
			mu_old[d*dim+x1] = mu[d*dim+x1];
			for(unsigned int x2=0;x2<dim;x2++){
				V_old[d*dim2+x1*dim+x2] = V[d*dim2+x1*dim+x2];
				V_inv_old[d*dim2+x1*dim+x2] = V_inv[d*dim2+x1*dim+x2];
			}
		}
	}
}

void GaussianMixtureModel::EM(unsigned int _nbClust){
	if( _nbClust > nbMaxClusters ){
		cerr << "The number of cluster for the expectation maximization algorithm is higher than the maximum number of cluster, aborting.." << endl;
		exit(EXIT_FAILURE);
	}
	nbClust = _nbClust;
	// Expectation
	for(unsigned int d=0;d<nbClust;d++) E_d[d] = 0.;
	for(unsigned int i=0;i<nbDat;i++){
		D_i[i] = 0.;
		for(unsigned int d=0;d<nbClust;d++){
			buffer_di[(d*nbDat)+i] = Prob_Cluster(d,i);
			D_i[i] += buffer_di[(d*nbDat)+i];
		}
		for(unsigned int d=0;d<nbClust;d++) C_di[(d*nbDat)+i] = buffer_di[(d*nbDat)+i] / D_i[i];
		for(unsigned int d=0;d<nbClust;d++) E_d[d] += C_di[(d*nbDat)+i];
	}
	// Maximization
	for(unsigned int d=0;d<nbClust;d++){
		// new weights
		weights[d] = E_d[d]/nbDat;
		// new means
		for(unsigned int x=0;x<dim;x++){
			mu[d*dim+x] = 0.;
			for(unsigned int i=0;i<nbDat;i++) mu[d*dim+x] += C_di[(d*nbDat)+i]*_MyDescriptors->getDescriptors()[i*dim+x];
			mu[d*dim+x] /= E_d[d];
		}
		// new variances
		for(unsigned int x1=0;x1<dim;x1++){
			for(unsigned int x2=x1;x2<dim;x2++){
				V[d*dim2+x1*dim+x2] = 0.;
				for(unsigned int i=0;i<nbDat;i++) V[d*dim2+x1*dim+x2] += C_di[(d*nbDat)+i]*( (_MyDescriptors->getDescriptors()[i*dim+x1] - mu[d*dim+x1]) * (_MyDescriptors->getDescriptors()[i*dim+x2] - mu[d*dim+x2]) );
				V[d*dim2+x1*dim+x2] /= E_d[d];
			}
			// symmetric part
			for(unsigned int x2=x1+1;x2<dim;x2++){
				V[d*dim2+x2*dim+x1] = V[d*dim2+x1*dim+x2];
			}
		}
		// invert variances
		MT->invMat_LU(V,V_inv,dim,d,det_V[d]);
	}
}

double GaussianMixtureModel::Prob_Cluster(unsigned int index_cluster, unsigned int DescriptorIndex){
	double sp = 0.;
	for(unsigned int j=0;j<dim;j++) buffer_vec_2_dim[j] = (_MyDescriptors->getDescriptors()[DescriptorIndex*dim+j]-mu[index_cluster*dim+j]);
	for(unsigned int i=0;i<dim;i++){
		buffer_vec_1_dim[i] = 0.;
		for(unsigned int j=0;j<dim;j++) buffer_vec_1_dim[i] += V_inv[index_cluster*dim2+i*dim+j]*buffer_vec_2_dim[j];
		sp += buffer_vec_1_dim[i]*buffer_vec_2_dim[i];
	}
	return ( weights[index_cluster] / ( pow(2.*M_PI, dim/2.) * sqrt(det_V[index_cluster]) ) ) * exp( -.5*sp );
}

void GaussianMixtureModel::ComputeLogLikelihood(){
	LogLikelihood = 0.;
	for(unsigned int d=0;d<nbClust;d++) LogLikelihood += Prob_Cluster(d,0);
	LogLikelihood = log(LogLikelihood);
	for(unsigned int i=1;i<nbDat;i++){
		long double sum = 0.;
		for(unsigned int d=0;d<nbClust;d++) sum += Prob_Cluster(d,i);
		LogLikelihood += log(sum);
	}
}

void GaussianMixtureModel::ComputeBIC(){
	unsigned int NbIndepParams = ( ( nbClust * ( ( dim*( dim + 1)/2) + dim + 1. ) ) - 1. );
	BIC = -( 2.*LogLikelihood ) + ( log((double) nbDat) * (double) NbIndepParams );
}

void GaussianMixtureModel::PrintModelParams(){
	cout << "The Gaussian Mixture Model is composed by " << nbClust << " clusters" << endl;
        for(unsigned int k=0;k<nbClust;k++){
                cout << "Cluster " << k << endl;
                cout << "weight = " << weights[k] << endl;
                cout << "mean = ";
                for(unsigned int i=0;i<dim;i++) cout << mu[k*dim+i] << " ";
                cout << endl;
                cout << "variance = " << endl;
                for(unsigned int i=0;i<dim;i++){
                        for(unsigned int j=0;j<dim;j++) cout << V[k*dim2+i*dim+j] << " ";
                        cout << endl;
                }
        }
	cout << endl;
}

GaussianMixtureModel::~GaussianMixtureModel(){
	delete[] buffer_vec_1_dim; 
	delete[] buffer_vec_2_dim; 
}
