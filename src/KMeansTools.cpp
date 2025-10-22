//**********************************************************************************
//*   GaussianMixtureModel.cpp                                                     *
//**********************************************************************************
//* This file contains the implementation of the GaussianMixtureModel class        *
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


#include "KMeansTools.h"
#include "AtomHicConfig.h"
#include <cmath>
#include <filesystem>
#include <dirent.h>
#include <random>

using namespace std; 

// constructor without descriptors
KMeansTools::KMeansTools(unsigned int &nbClust, MatrixXd *dataMat, unsigned int &nbDat, unsigned int &dim):_nbClust(nbClust), _dataMat(dataMat), _nbDat(nbDat), _dim(dim){
	
	if( _nbClust > _nbDat ){
		cerr << "The number of cluster of KMeans tools is higher than the number of data, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	
	readFixedParams();
	
	InitializeKMeansVariables();
	InitializeDataVariables();

}

KMeansTools::KMeansTools(unsigned int &nbClust, unsigned int &dim, vector<vector<double>> &centroids):_nbClust(nbClust), _dim(dim){
	readFixedParams();
	InitializeKMeansVariables();
	setKMeansVariables(centroids);
}

KMeansTools::KMeansTools(unsigned int &nbClust, MatrixXd *dataMat, unsigned int &nbDat, unsigned int &dim, vector<vector<double>> &centroids):_nbClust(nbClust), _dataMat(dataMat), _nbDat(nbDat), _dim(dim){
	readFixedParams();
	InitializeKMeansVariables();
	InitializeDataVariables();
	setKMeansVariables(centroids);
}

void KMeansTools::InitializeDataVariables(){
	_Data2Cluster = new unsigned int[_nbDat];
}

void KMeansTools::InitializeKMeansVariables(){
	_centroids = MatrixXd(_nbClust,_dim);
	_centroids_old = MatrixXd(_nbClust,_dim);
	_optimal_centroids = MatrixXd(_nbClust,_dim);
	_V = MatrixXd(_nbClust * _dim, _dim);
	_nbDat2Cluster = new unsigned int[_nbClust];
}

void KMeansTools::setDataMat(MatrixXd *dataMat, unsigned int &nbDat){
	_dataMat = dataMat;
	_nbDat = nbDat;
	if( _Data2Cluster ) delete[] _Data2Cluster;
	_Data2Cluster = new unsigned int[_nbDat];
}

void KMeansTools::setKMeansVariables(vector<vector<double>> &centroids){
	unsigned int nbClust = centroids.size();
	if( nbClust != _nbClust ){
		cerr << "The number of cluster of the KMeansTools is not consistent with the provided variables, aborting.." << endl;
		exit(EXIT_FAILURE);
	}
	unsigned int dim = centroids[0].size();
	if( dim != _dim ){
		cerr << "The number of dimension of the KMeansTools is not consistent with the provided variables, aborting.." << endl;
		exit(EXIT_FAILURE);
	}
	for(unsigned int k=0;k<_nbClust;k++)
		for(unsigned int d1=0;d1<_dim;d1++)
			_centroids(k,d1) = centroids[k][d1];
}

double KMeansTools::ComputeSilhouette(){
	vector<vector<unsigned int>> cluster_members(_nbClust);
	
	for(unsigned int i=0;i<_nbDat;i++) cluster_members[_Data2Cluster[i]].push_back(i);
	
	double silhouette_sum = 0.;

	#pragma omp parallel for	
	for(unsigned int i=0;i<_nbDat;i++){
		unsigned int cluster_i = _Data2Cluster[i];
		
		// Compute a(i)
		double a = 0.;
		if(cluster_members[cluster_i].size() > 1){
			for(unsigned int j : cluster_members[cluster_i])
		    		if (i != j) a += (_dataMat->row(i) - _dataMat->row(j)).norm();
		    	a /= (cluster_members[cluster_i].size() - 1);
		}else a = 0.; // singleton cluster
		
		// Compute b(i)
		double b = std::numeric_limits<double>::max();
		for(unsigned int k=0;k<_nbClust;k++){
			if(k == cluster_i || cluster_members[k].empty()) continue;
			double dist = 0.;
			for(unsigned int j : cluster_members[k]) dist += (_dataMat->row(i) - _dataMat->row(j)).norm();
			dist /= cluster_members[k].size();
			if (dist < b) b = dist;
		}
		
		double s = 0.;
		if(max(a, b) > 0.0){
		    s = (b - a) / max(a, b);
		}
		silhouette_sum += s;
	}
	
	return silhouette_sum / _nbDat;
}

void KMeansTools::AffectData2Cluster(){
	_inertia = 0.;
	MatrixXd dists(_nbClust,_nbDat);
	for(unsigned int k=0;k<_nbClust;k++) dists.row(k) = (_dataMat->rowwise() - _centroids.row(k)).rowwise().squaredNorm();
	VectorXd minVals(_nbDat);
	
	for(unsigned int i=0;i<_nbDat;i++){
		Index minIdx;
		minVals[i] = dists.col(i).minCoeff(&minIdx);
		_Data2Cluster[i] = static_cast<unsigned int>(minIdx);
	}
	
	_inertia = minVals.sum();
	
	for(unsigned int k=0;k<_nbClust;k++) _nbDat2Cluster[k] = 0;
	for(unsigned int i=0;i<_nbDat;i++) _nbDat2Cluster[_Data2Cluster[i]]++;
}

void KMeansTools::setSeed(unsigned int &seed){
        FixedSeed = true;
        this->seed = seed;
}

void KMeansTools::KMeansPPInitialization(){
	random_device rd;
	mt19937 gen;
    	if( FixedSeed ) gen.seed(seed);
	else gen.seed(rd());
	uniform_int_distribution<> distrib(0, _nbDat-1);
	unsigned int rand_index = distrib(gen);
	_centroids.row(0) = _dataMat->row(rand_index);
	_optimal_centroids.row(0) = _centroids.row(0);
	VectorXd min_distances(_nbDat);
	for(unsigned int k=1;k<_nbClust;k++){
		VectorXd new_dists = (_dataMat->rowwise() - _centroids.row(k-1)).rowwise().squaredNorm();
		if( k == 1 ) min_distances = new_dists;
		else min_distances = min_distances.cwiseMin(new_dists);
		
		Index max_idx;
		min_distances.maxCoeff(&max_idx);
		
		_centroids.row(k) = _dataMat->row(max_idx);
		_optimal_centroids.row(k) = _dataMat->row(max_idx);
	}
}

void KMeansTools::fit(){
	_optimal_inertia = numeric_limits<double>::max();
	for(unsigned int nb=0;nb<nbInit;nb++){
		if( nbInit > 1 ) cout << "\r" << nb+1 << "/" << nbInit << " KMeans initializations " << flush;
		KMeansPPInitialization();
		double res;
		unsigned int iter = 0;
		// Classical KMeans algorithm
		do{
			AffectData2Cluster();

			_centroids_old = _centroids;
			_centroids.setZero();
			
			for(unsigned int i=0;i<_nbDat;i++) _centroids.row(_Data2Cluster[i]) += _dataMat->row(i);
			
			for(unsigned int k=0;k<_nbClust;k++){
				if(_nbDat2Cluster[k] > 0) _centroids.row(k) /= _nbDat2Cluster[k];
			}
			
			res = 0.;
			for(unsigned int k=0;k<_nbClust;k++){
				double norm = _centroids.row(k).squaredNorm();
				if(norm > 0.){
					VectorXd diff = _centroids.row(k) - _centroids_old.row(k);
					res += diff.squaredNorm() / norm;
				}
			}
			iter++;
		}while( iter < MaxIter_KMeans && res > tol_KMeans );
			
		AffectData2Cluster();

		if( res > tol_KMeans ){
			cout << "KMeans clustering did not converged after " << MaxIter_KMeans << " iterations" << endl;
			cout << "Final residual : " << res << endl;
		        cout << "Maybe increase number of iteration or the tolerance for KMeans in /data/FixedParameters/FixedParameters.dat" << endl;
		}
		if( _inertia < _optimal_inertia ){
			_optimal_inertia = _inertia;
			_optimal_centroids = _centroids;
		}
	}
	if( nbInit > 1 ) cout << endl;
	ComputeFullVariances();
}

void KMeansTools::ComputeFullVariances(){

	for(unsigned int k=0;k<_nbClust;k++){
		if(_nbDat2Cluster[k] == 0) continue;
		
		vector<int> indices;
		for(unsigned int i=0;i<_nbDat;i++){ 
			if(_Data2Cluster[i] == k) indices.push_back(i);
		}
		
		MatrixXd cluster_data(indices.size(), _dim);
		for(size_t row=0;row<indices.size();row++) cluster_data.row(row) = _dataMat->row(indices[row]);
		
		MatrixXd centered = cluster_data.rowwise() - _optimal_centroids.row(k);
		
		MatrixXd cov = centered.transpose() * centered / _nbDat2Cluster[k];
		
		_V.block(k*_dim,0,_dim,_dim) = cov.selfadjointView<Upper>();
	}
}

void KMeansTools::readFixedParams(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_tol, pos_maxIter, pos_nbInit;
	string buffer_s, line;
	unsigned int ReadOk(0);
	if(file){
		while(file){
			getline(file,line);
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
	if( ReadOk != 3 ){
		cerr << "Error during reading of FixedParameters.dat for KMeansTools, aborting" << endl;
		exit(EXIT_FAILURE);
	}
}

void KMeansTools::ReadProperties(vector<string> Properties){
	size_t pos_tol, pos_maxIter, pos_nbInit;
	string buffer_s;
	for(unsigned int i=0;i<Properties.size();i++){
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


KMeansTools::~KMeansTools(){
	delete[] _Data2Cluster;
	delete[] _nbDat2Cluster;
}
