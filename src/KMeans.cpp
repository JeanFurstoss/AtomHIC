#include "KMeans.h"
#include "AtomHicConfig.h"
#include <cmath>

using namespace std; 

KMeans::KMeans(){
	this->name = "KMeans"; // TODO the user could affect name to the model or it could be read from database
}

void KMeans::setDescriptors(Descriptors *D){
	MachineLearningModel::setDescriptors(D);

	buffer_k = new double[nbClustMax];
	nbDat2Cluster = new unsigned int[nbClustMax];
	centroids = new double[nbClustMax*dim];
	V = new long double[nbClustMax*dim*dim];
	V_inv = new long double[nbClustMax*dim*dim];
	det_V = new long double[nbClustMax];
	centroids_old = new double[nbClustMax*dim];

	buffer_dat = new double[nbDat];
	Data2Cluster = new unsigned int[nbDat];
}

double KMeans::SquareEuclidianDistance(const unsigned int &DescriptorIndex, const unsigned int &ClusterIndex){
	double dist = 0.;
	for(unsigned int d=0;d<dim;d++) dist += ( _MyDescriptors->getDescriptors()[DescriptorIndex*dim+d] - centroids[ClusterIndex*dim+d] ) * ( _MyDescriptors->getDescriptors()[DescriptorIndex*dim+d] - centroids[ClusterIndex*dim+d] );
	return dist;
}

void KMeans::AffectData2Cluster(){
	for(unsigned int k=0;k<nbClust;k++) nbDat2Cluster[k] = 0;
	for(unsigned int i=0;i<nbDat;i++){
		for(unsigned int k=0;k<nbClust;k++) buffer_k[k] = SquareEuclidianDistance(i,k);
		Data2Cluster[i] = MT->min_p_ind(buffer_k,nbClust);
		nbDat2Cluster[Data2Cluster[i]] += 1;
	}
}

void KMeans::KMeansPPInitialization(unsigned int &_nbClust){
	if( _nbClust > nbClustMax ){
		cerr << "The number of cluster for KMeans is higher than the maximum number of cluster, aborting" << endl;
		exit(EXIT_FAILURE);
	}

	nbClust = _nbClust;

	srand(time(0));
	unsigned int rand_index = rand() % ( nbDat + 1 );
	// Initialize the first centroid from data points randomly
	for(unsigned int d=0;d<dim;d++) centroids[d] = _MyDescriptors->getDescriptors()[rand_index*dim+d];
	for(unsigned int k1=1;k1<nbClust;k1++){
		for(unsigned int p=0;p<nbDat;p++){
			// find out the closest cluster to the point
			for(unsigned int k2=0;k2<k1;k2++) buffer_k[k2] = SquareEuclidianDistance(p,k2);
			buffer_dat[p] = buffer_k[MT->min_p_ind(buffer_k,k1)];
		}
		// Initialize the next centroid with the data point with largest distance
		unsigned int ind = MT->max_p_ind(buffer_dat,nbDat);
		for(unsigned int d=0;d<dim;d++) centroids[k1*dim+d] = _MyDescriptors->getDescriptors()[MT->max_p_ind(buffer_dat,nbDat)*dim+d];
	}
}

void KMeans::TrainModel(unsigned int &_nbClust){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, training aborted" << endl;
		exit(EXIT_FAILURE);
	}

        if( !SavedVariablesInitialized ){
                saved_LogLikelihood = new long double[nbInit];
                saved_nbDat2Cluster = new unsigned int[nbClustMax*nbInit];
                saved_centroids = new double[nbClustMax*nbInit*dim];
                saved_V = new long double[dim2*nbClustMax*nbInit];
                saved_V_inv = new long double[dim2*nbClustMax*nbInit];
                saved_det_V = new long double[nbClustMax*nbInit];
                SavedVariablesInitialized = true;
        }

	for(unsigned int nb=0;nb<nbInit;nb++){
		KMeansPPInitialization(_nbClust);
		double res;
		unsigned int iter = 0;
		// Classical KMeans algorithm
		do{
			AffectData2Cluster();
			for(unsigned int k=0;k<nbClust;k++){
				for(unsigned int d=0;d<dim;d++){
					centroids_old[k*dim+d] = centroids[k*dim+d];
					centroids[k*dim+d] = 0.;
				}
			}
			for(unsigned int i=0;i<nbDat;i++){
				for(unsigned int d=0;d<dim;d++) centroids[Data2Cluster[i]*dim+d] += _MyDescriptors->getDescriptors()[i*dim+d];
			}
			res = 0.;
			for(unsigned int k=0;k<nbClust;k++){
				for(unsigned int d=0;d<dim;d++){
					centroids[k*dim+d] /= nbDat2Cluster[k];
					res += (centroids[k*dim+d]-centroids_old[k*dim+d]) * (centroids[k*dim+d]-centroids_old[k*dim+d]);
				}
			}
			iter++;
		}while( iter < MaxIter_KMeans && res > tol_KMeans );
			
		AffectData2Cluster();
		ComputeFullVariances();
		ComputeLogLikelihood();
		ComputeBIC();
		SaveVariables(nb);

		if( res > tol_KMeans ){
			cout << "KMeans clustering did not converged after " << MaxIter_KMeans << " iterations" << endl;
		        cout << "Maybe increase number of iteration or the tolerance for KMeans in /data/FixedParameters/FixedParameters.dat" << endl;
		}
	}

	unsigned int index_opt = MT->max_p_ind(saved_LogLikelihood,nbInit);
	LogLikelihood = saved_LogLikelihood[index_opt];
	for(unsigned int k=0;k<nbClust;k++){
		nbDat2Cluster[k]  = saved_nbDat2Cluster[index_opt*nbClustMax+k];
		det_V[k] = saved_det_V[index_opt*nbClustMax+k];
		for(unsigned int d1=0;d1<dim;d1++){
			centroids[k*dim+d1] = saved_centroids[index_opt*nbClustMax*dim+k*dim+d1];
			for(unsigned int d2=0;d2<dim;d2++){
				V[k*dim2+d1*dim+d2] = saved_V[index_opt*nbClustMax*dim2+k*dim2+d1*dim+d2];
				V_inv[k*dim2+d1*dim+d2] = saved_V_inv[index_opt*nbClustMax*dim2+k*dim2+d1*dim+d2];
			}
		}
	}
	ComputeBIC();
}

void KMeans::SaveVariables(unsigned int &current){
	saved_LogLikelihood[current] = LogLikelihood;
	for(unsigned int k=0;k<nbClust;k++){
		saved_nbDat2Cluster[current*nbClustMax+k] = nbDat2Cluster[k];
		saved_det_V[current*nbClustMax+k] = det_V[k];
		for(unsigned int d1=0;d1<dim;d1++){
			saved_centroids[current*nbClustMax*dim+k*dim+d1] = centroids[k*dim+d1];
			for(unsigned int d2=0;d2<dim;d2++){
				saved_V[current*nbClustMax*dim2+k*dim2+d1*dim+d2] = V[k*dim2+d1*dim+d2];
				saved_V_inv[current*nbClustMax*dim2+k*dim2+d1*dim+d2] = V_inv[k*dim2+d1*dim+d2];
			}
		}
	}
}

void KMeans::ComputeFullVariances(){
	for(unsigned int k=0;k<nbClust;k++){
		for(unsigned int d1=0;d1<dim;d1++){
			for(unsigned int d2=d1;d2<dim;d2++){
				V[k*dim2+d1*dim+d2] = 0.;
				for(unsigned int i=0;i<nbDat;i++) V[k*dim2+d1*dim+d2] += ( _MyDescriptors->getDescriptors()[i*dim+d1] - centroids[k*dim+d1] ) * ( _MyDescriptors->getDescriptors()[i*dim+d2] - centroids[k*dim+d2] );
				V[k*dim2+d1*dim+d2] /= nbDat;
			}
			for(unsigned int d2=d1+1;d2<dim;d2++) V[k*dim2+d2*dim+d1] = V[k*dim2+d1*dim+d2];
		}
		MT->invMat_LU(V,V_inv,dim,k,det_V[k]);
	}
}

double KMeans::ComputeGaussianProb(unsigned int &DescriptorIndex){
        double sp = 0.;
        for(unsigned int j=0;j<dim;j++) buffer_vec_2_dim[j] = (_MyDescriptors->getDescriptors()[DescriptorIndex*dim+j]-centroids[Data2Cluster[DescriptorIndex]*dim+j]);
        for(unsigned int i=0;i<dim;i++){
                buffer_vec_1_dim[i] = 0.;
                for(unsigned int j=0;j<dim;j++) buffer_vec_1_dim[i] += V_inv[Data2Cluster[DescriptorIndex]*dim2+i*dim+j]*buffer_vec_2_dim[j];
                sp += buffer_vec_1_dim[i]*buffer_vec_2_dim[i];
        }
        return ( ( (double) nbDat2Cluster[Data2Cluster[DescriptorIndex]] / (double) nbDat ) / ( pow(2.*M_PI, (double) dim/2.) * sqrt(det_V[Data2Cluster[DescriptorIndex]]) ) ) * exp( -.5*sp );
}

// Here we compute the log likelihood by considering a mixture of gaussian distribution
void KMeans::ComputeLogLikelihood(){
	LogLikelihood = 0.;
	for(unsigned int i=1;i<nbDat;i++) LogLikelihood += log(ComputeGaussianProb(i));
}

void KMeans::ComputeBIC(){
	BIC = -( 2.*LogLikelihood ) + ( log((double) nbDat) * (double) dim * (double) nbClust );
}

void KMeans::PrintModelParams(){
	cout << "The KMean clustering leads to " << nbClust << " clusters" << endl;
        for(unsigned int k=0;k<nbClust;k++){
                cout << "Cluster " << k << endl;
		cout << "Weight = " << (double) nbDat2Cluster[k] / nbDat << endl;
                cout << "centroid = ";
                for(unsigned int i=0;i<dim;i++) cout << centroids[k*dim+i] << " ";
                cout << endl;
        }
}

KMeans::~KMeans(){
	if( IsDescriptor ){
		delete[] buffer_k;
		delete[] nbDat2Cluster;
		delete[] centroids;
		delete[] V;
		delete[] V_inv;
		delete[] det_V;
		delete[] centroids_old;
		delete[] buffer_dat;
		delete[] Data2Cluster;
	}
        if( !SavedVariablesInitialized ){
                delete[] saved_LogLikelihood;
                delete[] saved_nbDat2Cluster;
                delete[] saved_centroids;
                delete[] saved_V;
                delete[] saved_V_inv;
                delete[] saved_det_V;
        }
}
