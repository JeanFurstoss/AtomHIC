#include "KMeans.h"
#include "AtomHicConfig.h"
#include <cmath>

using namespace std; 

KMeans::KMeans(){
	this->name = "KMeans"; // TODO the user could affect name to the model or it could be read from database
}

void KMeans::setDescriptors(Descriptors *D){
	MachineLearningModel::setDescriptors(D);

	BIC = new double[nbFilter];
	LogLikelihood = new long double[nbFilter];
	nbClust = new unsigned int[nbFilter];

	centroids = new double[nbClustMax*dim*nbFilter];
	V = new long double[nbClustMax*dim*dim*nbFilter];
	V_inv = new long double[nbClustMax*dim*dim*nbFilter];
	det_V = new long double[nbClustMax*nbFilter];
	centroids_old = new double[nbClustMax*dim*nbFilter];

	nbDat2Cluster = new unsigned int[nbClustMax*nbFilter];
	Data2Cluster = new unsigned int[nbDatMax*nbFilter];
	
	buffer_k = new double[nbClustMax];
	buffer_dat = new double[nbDatMax];
	nbDat = new unsigned int[nbFilter];
	for(unsigned int f=0;f<nbFilter;f++) nbDat[f] = _MyDescriptors->getNbDat(f);
}

double KMeans::SquareEuclidianDistance(const unsigned int &DescriptorIndex, const unsigned int &ClusterIndex, unsigned int &filter_value){
	double dist = 0.;
	for(unsigned int d=0;d<dim;d++) dist += ( _MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+DescriptorIndex*dim+d] - centroids[ClusterIndex*dim*nbFilter+d*nbFilter+filter_value] ) * ( _MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+DescriptorIndex*dim+d] - centroids[ClusterIndex*dim*nbFilter+d*nbFilter+filter_value] );
	return dist;
}

void KMeans::AffectData2Cluster(unsigned int &filter_value){
	for(unsigned int k=0;k<nbClust[filter_value];k++) nbDat2Cluster[k*nbFilter+filter_value] = 0;
	for(unsigned int i=0;i<nbDat[filter_value];i++){
		for(unsigned int k=0;k<nbClust[filter_value];k++) buffer_k[k] = SquareEuclidianDistance(i,k,filter_value);
		Data2Cluster[i*nbFilter+filter_value] = MT->min_p_ind(buffer_k,nbClust[filter_value]);
		nbDat2Cluster[Data2Cluster[i*nbFilter+filter_value]*nbFilter+filter_value] += 1;
	}
}

void KMeans::KMeansPPInitialization(unsigned int &_nbClust, unsigned int &filter_value){
	if( _nbClust > nbClustMax ){
		cerr << "The number of cluster for KMeans is higher than the maximum number of cluster, aborting" << endl;
		exit(EXIT_FAILURE);
	}

	nbClust[filter_value] = _nbClust;

	srand(time(0));
	unsigned int rand_index = rand() % ( nbDat[filter_value] + 1 );
	// Initialize the first centroid from data points randomly
	for(unsigned int d=0;d<dim;d++) centroids[d*nbFilter+filter_value] = _MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+rand_index*dim+d];
	for(unsigned int k1=1;k1<nbClust[filter_value];k1++){
		for(unsigned int p=0;p<nbDat[filter_value];p++){
			// find out the closest cluster to the point
			for(unsigned int k2=0;k2<k1;k2++) buffer_k[k2] = SquareEuclidianDistance(p,k2,filter_value);
			buffer_dat[p] = buffer_k[MT->min_p_ind(buffer_k,k1)];
		}
		// Initialize the next centroid with the data point with largest distance
		unsigned int ind = MT->max_p_ind(buffer_dat,nbDat[filter_value]);
		for(unsigned int d=0;d<dim;d++) centroids[k1*dim*nbFilter+d*nbFilter+filter_value] = _MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+MT->max_p_ind(buffer_dat,nbDat[filter_value])*dim+d];
	}
}

void KMeans::TrainModel(unsigned int &_nbClust, unsigned int &filter_value){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, training aborted" << endl;
		exit(EXIT_FAILURE);
	}

        if( !SavedVariablesInitialized ){
                saved_LogLikelihood = new long double[nbInit*nbFilter];
                saved_nbDat2Cluster = new unsigned int[nbClustMax*nbInit*nbFilter];
                saved_centroids = new double[nbClustMax*nbInit*dim*nbFilter];
                saved_V = new long double[dim2*nbClustMax*nbInit*nbFilter];
                saved_V_inv = new long double[dim2*nbClustMax*nbInit*nbFilter];
                saved_det_V = new long double[nbClustMax*nbInit*nbFilter];
                SavedVariablesInitialized = true;
        }

	for(unsigned int nb=0;nb<nbInit;nb++){
		KMeansPPInitialization(_nbClust,filter_value);
		double res;
		unsigned int iter = 0;
		// Classical KMeans algorithm
		do{
			AffectData2Cluster(filter_value);
			for(unsigned int k=0;k<nbClust[filter_value];k++){
				for(unsigned int d=0;d<dim;d++){
					centroids_old[k*dim*nbFilter+d*nbFilter+filter_value] = centroids[k*dim*nbFilter+d*nbFilter+filter_value];
					centroids[k*dim*nbFilter+d*nbFilter+filter_value] = 0.;
				}
			}
			for(unsigned int i=0;i<nbDat[filter_value];i++){
				for(unsigned int d=0;d<dim;d++) centroids[Data2Cluster[i*nbFilter+filter_value]*dim*nbFilter+d*nbFilter+filter_value] += _MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+i*dim+d];
			}
			res = 0.;
			for(unsigned int k=0;k<nbClust[filter_value];k++){
				for(unsigned int d=0;d<dim;d++){
					centroids[k*dim*nbFilter+d*nbFilter+filter_value] /= nbDat2Cluster[k*nbFilter+filter_value];
					res += (centroids[k*dim*nbFilter+d*nbFilter+filter_value]-centroids_old[k*dim*nbFilter+d*nbFilter+filter_value]) * (centroids[k*dim*nbFilter+d*nbFilter+filter_value]-centroids_old[k*dim*nbFilter+d*nbFilter+filter_value]);
				}
			}
			iter++;
		}while( iter < MaxIter_KMeans && res > tol_KMeans );
			
		AffectData2Cluster(filter_value);
		ComputeFullVariances(filter_value);
		ComputeLogLikelihood(filter_value);
		ComputeBIC(filter_value);
		SaveVariables(nb,filter_value);

		if( res > tol_KMeans ){
			cout << "KMeans clustering did not converged after " << MaxIter_KMeans << " iterations" << endl;
		        cout << "Maybe increase number of iteration or the tolerance for KMeans in /data/FixedParameters/FixedParameters.dat" << endl;
		}
	}
	
	unsigned int index_opt = 0;
	for(unsigned int i=1;i<nbInit;i++){
		if( saved_LogLikelihood[i*nbFilter+filter_value] > saved_LogLikelihood[index_opt*nbFilter+filter_value] ) index_opt = i;
	}
	LogLikelihood[filter_value] = saved_LogLikelihood[index_opt*nbFilter+filter_value];
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		nbDat2Cluster[k*nbFilter+filter_value]  = saved_nbDat2Cluster[index_opt*nbClustMax*nbFilter+k*nbFilter+filter_value];
		det_V[k*nbFilter+filter_value] = saved_det_V[index_opt*nbClustMax*nbFilter+k*nbFilter+filter_value];
		for(unsigned int d1=0;d1<dim;d1++){
			centroids[k*dim*nbFilter+d1*nbFilter+filter_value] = saved_centroids[index_opt*nbClustMax*dim*nbFilter+k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = saved_V[index_opt*nbClustMax*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
				V_inv[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = saved_V_inv[index_opt*nbClustMax*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
			}
		}
	}
	ComputeBIC(filter_value);
}

void KMeans::SaveVariables(unsigned int &current, unsigned int &filter_value){
	saved_LogLikelihood[current*nbFilter+filter_value] = LogLikelihood[filter_value];
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		saved_nbDat2Cluster[current*nbClustMax*nbFilter+k*nbFilter+filter_value] = nbDat2Cluster[k*nbFilter+filter_value];
		saved_det_V[current*nbClustMax*nbFilter+k*nbFilter+filter_value] = det_V[k*nbFilter+filter_value];
		for(unsigned int d1=0;d1<dim;d1++){
			saved_centroids[current*nbClustMax*dim*nbFilter+k*dim*nbFilter+d1*nbFilter+filter_value] = centroids[k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				saved_V[current*nbClustMax*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
				saved_V_inv[current*nbClustMax*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = V_inv[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
			}
		}
	}
}

void KMeans::ComputeFullVariances(unsigned int &filter_value){
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		for(unsigned int d1=0;d1<dim;d1++){
			for(unsigned int d2=d1;d2<dim;d2++){
				V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = 0.;
				for(unsigned int i=0;i<nbDat[filter_value];i++) V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] += ( _MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+i*dim+d1] - centroids[k*dim*nbFilter+d1*nbFilter+filter_value] ) * ( _MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+i*dim+d2] - centroids[k*dim*nbFilter+d2*nbFilter+filter_value] );
				V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] /= nbDat[filter_value];
			}
			for(unsigned int d2=d1+1;d2<dim;d2++) V[k*dim2*nbFilter+d2*dim*nbFilter+d1*nbFilter+filter_value] = V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
		}
		MT->invMat_LU(V,V_inv,dim,k,nbFilter,filter_value,det_V[k*nbFilter+filter_value]);
	}
}

double KMeans::ComputeGaussianProb(unsigned int &DescriptorIndex, unsigned int &filter_value){
        double sp = 0.;
        for(unsigned int j=0;j<dim;j++) buffer_vec_2_dim[j] = (_MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+DescriptorIndex*dim+j]-centroids[Data2Cluster[DescriptorIndex*nbFilter+filter_value]*dim*nbFilter+j*nbFilter+filter_value]);
        for(unsigned int i=0;i<dim;i++){
                buffer_vec_1_dim[i] = 0.;
                for(unsigned int j=0;j<dim;j++) buffer_vec_1_dim[i] += V_inv[Data2Cluster[DescriptorIndex*nbFilter+filter_value]*dim2*nbFilter+i*dim*nbFilter+j*nbFilter+filter_value]*buffer_vec_2_dim[j];
                sp += buffer_vec_1_dim[i]*buffer_vec_2_dim[i];
        }
        return ( ( (double) nbDat2Cluster[Data2Cluster[DescriptorIndex*nbFilter+filter_value]*nbFilter+filter_value] / (double) nbDat[filter_value] ) / ( pow(2.*M_PI, (double) dim/2.) * sqrt(det_V[Data2Cluster[DescriptorIndex*nbFilter+filter_value]*nbFilter+filter_value]) ) ) * exp( -.5*sp );
}

// Here we compute the log likelihood by considering a mixture of gaussian distribution
void KMeans::ComputeLogLikelihood(unsigned int &filter_value){
	LogLikelihood[filter_value] = 0.;
	for(unsigned int i=1;i<nbDat[filter_value];i++) LogLikelihood[filter_value] += log(ComputeGaussianProb(i,filter_value));
}

void KMeans::ComputeBIC(unsigned int &filter_value){
	BIC[filter_value] = -( 2.*LogLikelihood[filter_value] ) + ( log((double) nbDat[filter_value]) * (double) dim * (double) nbClust[filter_value] );
}

void KMeans::PrintModelParams(unsigned int &filter_value){
	cout << "The KMean clustering leads to " << nbClust[filter_value] << " clusters" << endl;
        for(unsigned int k=0;k<nbClust[filter_value];k++){
                cout << "Cluster " << k << endl;
		cout << "Weight = " << (double) nbDat2Cluster[k*nbFilter+filter_value] / nbDat[filter_value] << endl;
                cout << "centroid = ";
                for(unsigned int i=0;i<dim;i++) cout << centroids[k*dim*nbFilter+i*nbFilter+filter_value] << " ";
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
		delete[] BIC;
		delete[] LogLikelihood;
		delete[] nbClust;
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
