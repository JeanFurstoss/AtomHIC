//**********************************************************************************
//*   KMeans.cpp                                                                   *
//**********************************************************************************
//* This file contains the implementation of the KMeans class                      *
//**********************************************************************************
//* (C) Jan 2025 - Jean Furstoss                                                   *
//*     Université de Poitiers, Institut PPRIME                                    *
//*     UPR CNRS 3346, 86360 Chasseuneuil-du-Poitou, France                        *
//*     jean.furstoss@univ-poitiers.fr                                             *
//* Last modification: J. Furstoss - 28 Janv 2025                                  *
//**********************************************************************************
//* This program is free software: you can redistribute it and/or modify           *
//* it under the terms of the GNU General Public License as published by           *
//* the Free Software Foundation, either version 3 of the License, or              *
//* (at your option) any later version.                                            *
//*                                                                                *
//* This program is distributed in the hope that it will be useful,                *
//* but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *
//* GNU General Public License for more details.                                   *
//*                                                                                *
//* You should have received a copy of the GNU General Public License              *
//* along with this program.  If not, see <http://www.gnu.org/licenses/>.          *
//**********************************************************************************


#include "KMeans.h"
#include "AtomHicConfig.h"
#include <cmath>

using namespace std; 

KMeans::KMeans(){
	this->name = "KMeans"; // TODO the user could affect name to the model or it could be read from database
	readFixedParams();
}

void KMeans::setDescriptors(Descriptors *D){
	MachineLearningModel::setDescriptors(D);
	if( this->IsDescriptor ){
		delete[] BIC;
		delete[] LogLikelihood;
		delete[] nbClust;
		delete[] centroids;
		delete[] V;
		delete[] V_inv;
		delete[] det_V;
		delete[] centroids_old;
		delete[] nbDat2Cluster;
		delete[] Data2Cluster;
		delete[] buffer_k;
		delete[] buffer_dat;
	}
	
	this->IsDescriptor = true;
	
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
}

double KMeans::SquareEuclidianDistance(const unsigned int &DescriptorIndex, const unsigned int &ClusterIndex, unsigned int &filter_value){
	double dist = 0.;
	double temp;
	for(unsigned int d=0;d<dim;d++){
		temp = ( _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+DescriptorIndex)*dim+d] - centroids[ClusterIndex*dim*nbFilter+d*nbFilter+filter_value] );
		dist += temp*temp;
	}
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
	for(unsigned int d=0;d<dim;d++) centroids[d*nbFilter+filter_value] = _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+rand_index)*dim+d];
	for(unsigned int k1=1;k1<nbClust[filter_value];k1++){
		for(unsigned int p=0;p<nbDat[filter_value];p++){
			// find out the closest cluster to the point
			for(unsigned int k2=0;k2<k1;k2++) buffer_k[k2] = SquareEuclidianDistance(p,k2,filter_value);
			buffer_dat[p] = buffer_k[MT->min_p_ind(buffer_k,k1)];
		}
		// Initialize the next centroid with the data point with largest distance
		unsigned int ind = MT->max_p_ind(buffer_dat,nbDat[filter_value]);
		for(unsigned int d=0;d<dim;d++) centroids[k1*dim*nbFilter+d*nbFilter+filter_value] = _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+MT->max_p_ind(buffer_dat,nbDat[filter_value]))*dim+d];
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
					unsigned int ind = k*dim*nbFilter+d*nbFilter+filter_value;
					centroids_old[ind] = centroids[ind];
					centroids[ind] = 0.;
				}
			}
			for(unsigned int i=0;i<nbDat[filter_value];i++){
				for(unsigned int d=0;d<dim;d++) centroids[Data2Cluster[i*nbFilter+filter_value]*dim*nbFilter+d*nbFilter+filter_value] += _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+i)*dim+d];
			}
			res = 0.;
			for(unsigned int k=0;k<nbClust[filter_value];k++){
				double norm(0.),res_temp(0.);
				for(unsigned int d=0;d<dim;d++){
					unsigned int ind = k*dim*nbFilter+d*nbFilter+filter_value;
					centroids[ind] /= nbDat2Cluster[k*nbFilter+filter_value];
					res_temp += (centroids[ind]-centroids_old[ind]) * (centroids[ind]-centroids_old[ind]);
					norm += centroids[ind]*centroids[ind];
				}
				if( norm != 0. ) res += res_temp/norm;
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
			cout << "Final residual : " << res << endl;
		        cout << "Maybe increase number of iteration or the tolerance for KMeans in /data/FixedParameters/FixedParameters.dat" << endl;
		}
	}
	
	unsigned int index_opt = 0;
	for(unsigned int i=1;i<nbInit;i++){
		if( saved_LogLikelihood[i*nbFilter+filter_value] > saved_LogLikelihood[index_opt*nbFilter+filter_value] ) index_opt = i;
	}
	LogLikelihood[filter_value] = saved_LogLikelihood[index_opt*nbFilter+filter_value];
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		unsigned int ind = index_opt*nbClustMax*nbFilter+k*nbFilter+filter_value;
		nbDat2Cluster[k*nbFilter+filter_value]  = saved_nbDat2Cluster[ind];
		det_V[k*nbFilter+filter_value] = saved_det_V[ind];
		for(unsigned int d1=0;d1<dim;d1++){
			centroids[k*dim*nbFilter+d1*nbFilter+filter_value] = saved_centroids[index_opt*nbClustMax*dim*nbFilter+k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				unsigned int ind1 = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				unsigned int ind2 = index_opt*nbClustMax*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				V[ind1] = saved_V[ind2];
				V_inv[ind1] = saved_V_inv[ind2];
			}
		}
	}
	ComputeBIC(filter_value);
	//cout << "AfterKMEAN" << endl;
	//for(unsigned int k1=0;k1<nbClust[filter_value];k1++){
	//	for(unsigned int d=0;d<dim;d++) cout << centroids[k1*dim*nbFilter+d*nbFilter+filter_value] << " ";
	//	cout << endl;
	//}
}

void KMeans::SaveVariables(unsigned int &current, unsigned int &filter_value){
	saved_LogLikelihood[current*nbFilter+filter_value] = LogLikelihood[filter_value];
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		unsigned int ind = current*nbClustMax*nbFilter+k*nbFilter+filter_value;
		saved_nbDat2Cluster[ind] = nbDat2Cluster[k*nbFilter+filter_value];
		saved_det_V[ind] = det_V[k*nbFilter+filter_value];
		for(unsigned int d1=0;d1<dim;d1++){
			saved_centroids[current*nbClustMax*dim*nbFilter+k*dim*nbFilter+d1*nbFilter+filter_value] = centroids[k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				unsigned int ind1 = current*nbClustMax*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				unsigned int ind2 = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				saved_V[ind1] = V[ind2];
				saved_V_inv[ind1] = V_inv[ind2];
			}
		}
	}
}

void KMeans::ComputeFullVariances(unsigned int &filter_value){
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		for(unsigned int d1=0;d1<dim;d1++){
			for(unsigned int d2=d1;d2<dim;d2++) V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = 0.;
		}
	}
	for(unsigned int d1=0;d1<dim;d1++){
		for(unsigned int d2=d1;d2<dim;d2++){
			for(unsigned int i=0;i<nbDat[filter_value];i++){
				V[Data2Cluster[i*nbFilter+filter_value]*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] += ( _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+i)*dim+d1] - centroids[Data2Cluster[i*nbFilter+filter_value]*dim*nbFilter+d1*nbFilter+filter_value] ) * ( _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+i)*dim+d2] - centroids[Data2Cluster[i*nbFilter+filter_value]*dim*nbFilter+d2*nbFilter+filter_value] );
			}
		}
	}
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		for(unsigned int d1=0;d1<dim;d1++){
			for(unsigned int d2=d1;d2<dim;d2++){
				unsigned int ind = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				if( nbDat2Cluster[k*nbFilter+filter_value] != 0 ) V[ind] /= nbDat2Cluster[k*nbFilter+filter_value];
				else V[ind] = 0.;
				V[k*dim2*nbFilter+d2*dim*nbFilter+d1*nbFilter+filter_value] = V[ind];
			}
		}
		MT->invMat_LU(V,V_inv,dim,k,nbFilter,filter_value,det_V[k*nbFilter+filter_value]);
	}
}

double KMeans::ComputeGaussianProb(unsigned int &DescriptorIndex, unsigned int &filter_value){
        double sp = 0.;
        for(unsigned int j=0;j<dim;j++) buffer_vec_2_dim[j] = (_MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+DescriptorIndex)*dim+j]-centroids[Data2Cluster[DescriptorIndex*nbFilter+filter_value]*dim*nbFilter+j*nbFilter+filter_value]);
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

void KMeans::readFixedParams(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_nbCMax, pos_tol, pos_maxIter, pos_nbInit;
	string buffer_s, line;
	unsigned int ReadOk(0);
	if(file){
		while(file){
			getline(file,line);
			pos_nbCMax=line.find("KMEANS_NB_MAX_CLUSTER");
			if(pos_nbCMax!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nbClustMax;
				ReadOk++;
			}
			pos_tol=line.find("KMEANS_TOL");
			if(pos_tol!=string::npos){
				istringstream text(line);
				text >> buffer_s >> tol_KMeans;
				ReadOk++;
			}
			pos_maxIter=line.find("KMEANS_MAX_ITER");
			if(pos_maxIter!=string::npos){
				istringstream text(line);
				text >> buffer_s >> MaxIter_KMeans;
				ReadOk++;
			}
			pos_nbInit=line.find("KMEANS_NB_INIT");
			if(pos_nbInit!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nbInit;
				ReadOk++;
			}
		}
	}else{
		cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
		exit(EXIT_FAILURE);
	}
	file.close();
	if( ReadOk != 4 ){
		cerr << "Error during reading of FixedParameters.dat for KMeans, aborting" << endl;
		exit(EXIT_FAILURE);
	}
}

void KMeans::ReadProperties(vector<string> Properties){
	size_t pos_nbCMax, pos_tol, pos_maxIter, pos_nbInit;
	string buffer_s;
	for(unsigned int i=0;i<Properties.size();i++){
		pos_nbCMax=Properties[i].find("KMEANS_NB_MAX_CLUSTER");
		if(pos_nbCMax!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> nbClustMax;
		}
		pos_tol=Properties[i].find("KMEANS_TOL");
		if(pos_tol!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> tol_KMeans;
		}
		pos_maxIter=Properties[i].find("KMEANS_MAX_ITER");
		if(pos_maxIter!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> MaxIter_KMeans;
		}
		pos_nbInit=Properties[i].find("KMEANS_NB_INIT");
		if(pos_nbInit!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> nbInit;
		}
	}
}

KMeans::~KMeans(){
	if( this->IsDescriptor ){
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
