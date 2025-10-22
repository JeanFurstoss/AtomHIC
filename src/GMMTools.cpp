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


#include "GMMTools.h"
#include "AtomHicConfig.h"
#include <cmath>
#include <filesystem>
#include <dirent.h>

using namespace std; 
using namespace std::chrono;

GMMTools::GMMTools(unsigned int &nbClust, MatrixXd *dataMat, unsigned int &nbDat, unsigned int &dim):_nbClust(nbClust), _dataMat(dataMat), _nbDat(nbDat), _dim(dim){
	
	if( _nbClust > _nbDat ){
		cerr << "The number of cluster of GMM tools is higher than the number of data, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	
	readFixedParams();
	InitializeGMMVariables();
	InitializeDataVariables();

	if( InitMethod == "KMEANS" || InitMethod == "KMEANSPP" ){
		_MyKMeansTools = new KMeansTools(_nbClust, _dataMat, _nbDat, _dim);
		IsKMeans = true;
	}
}

GMMTools::GMMTools(unsigned int &nbClust, unsigned int &dim, vector<double> &weights, vector<vector<double>> &mu, vector<vector<double>> &V):_nbClust(nbClust), _dim(dim){
	readFixedParams();
	InitializeGMMVariables();
	setGMMVariables(weights,mu,V);
}

GMMTools::GMMTools(unsigned int &nbClust, MatrixXd *dataMat, unsigned int &nbDat, unsigned int &dim, vector<double> &weights, vector<vector<double>> &mu, vector<vector<double>> &V):_nbClust(nbClust), _dataMat(dataMat), _nbDat(nbDat), _dim(dim){
	readFixedParams();
	InitializeGMMVariables();
	InitializeDataVariables();
	setGMMVariables(weights,mu,V);
}

void GMMTools::InitializeDataVariables(){
	_WeightedLogProb = VectorXd(_nbDat);
	_WeightedClusterLogProb = MatrixXd(_nbClust,_nbDat);
	_Resp = MatrixXd(_nbClust,_nbDat);
}

void GMMTools::InitializeGMMVariables(){
	_weights = VectorXd(_nbClust);
	_mu = MatrixXd(_nbClust,_dim);
	_V = MatrixXd(_nbClust * _dim ,_dim);
	_V_ldlt = vector<LDLT<MatrixXd>>(_nbClust);
	_det_V = vector<long double>(_nbClust,0.);

	_optimal_weights = VectorXd(_nbClust);
	_optimal_mu = MatrixXd(_nbClust,_dim);
	_optimal_V = MatrixXd(_nbClust * _dim ,_dim);
	_optimal_V_ldlt = vector<LDLT<MatrixXd>>(_nbClust);
	_optimal_det_V = vector<long double>(_nbClust,0.);

	_nbDat2Cluster = VectorXd(_nbClust);
	
	logFacNorm = -.5*_dim*log(2.*M_PI);
	regcovar = 1e-6 * MatrixXd::Identity(_dim, _dim);
}

void GMMTools::setGMMVariables(vector<double> &weights, vector<vector<double>> &mu, vector<vector<double>> &V){
	unsigned int nbClust = weights.size();
	if( nbClust != _nbClust ){
		cerr << "The number of cluster of the GMMTools is not consistent with the provided variables, aborting.." << endl;
		exit(EXIT_FAILURE);
	}
	unsigned int dim = mu[0].size();
	if( dim != _dim ){
		cerr << "The number of dimension of the GMMTools is not consistent with the provided variables, aborting.." << endl;
		exit(EXIT_FAILURE);
	}
	for(unsigned int k=0;k<_nbClust;k++){
		_weights[k] = weights[k];
		MatrixXd tempmat(_dim,_dim);
		for(unsigned int d1=0;d1<_dim;d1++){
			_mu(k,d1) = mu[k][d1];
			for(unsigned int d2=0;d2<_dim;d2++) tempmat(d1,d2) = V[k][d1*_dim+d2];
		}
		tempmat += regcovar;
		_V.block(k * _dim, 0, _dim, _dim) = tempmat;
		_V_ldlt[k] = tempmat.ldlt();
		_det_V[k] = _V_ldlt[k].vectorD().array().log().sum();
	}

}

void GMMTools::setDataMat(MatrixXd *dataMat, unsigned int &nbDat){
	_dataMat = dataMat;
	_nbDat = nbDat;
	
	_WeightedLogProb = VectorXd(_nbDat);
	_WeightedClusterLogProb = MatrixXd(_nbClust,_nbDat);
	_Resp = MatrixXd(_nbClust,_nbDat);
}

void GMMTools::fit(){
	_optimal_LogLikelihood = -numeric_limits<double>::max();
	bool converged;
	double old_lkh;
	double residual;
	for(unsigned int init=0;init<nbInit;init++){
		cout << "\r" << init+1 << "/" << nbInit << " GMM initializations " << flush;
		Initialize();
		converged = false;
		old_lkh = expectation();
		unsigned int niter = 0;
		for(unsigned int n=0;n<MaxIter_EM;n++){
			maximization();
			_LogLikelihood = expectation();
			residual = fabs(old_lkh - _LogLikelihood);
			if( residual < tol_Lkh_EM ){
			       converged = true;
			       niter = n+1;
			       break;
		       	}else old_lkh = _LogLikelihood;
		}
		if( converged ){
			//cout << "EM converged in " << niter << " steps" << endl;
			if( _LogLikelihood > _optimal_LogLikelihood ) SaveOptimalRun();
		}else cout << "EM algo does not converge (residual = " << residual << ")" << endl;
	}
	cout << endl;
	SetOptimalRun();
	//for(unsigned int k=0;k<_nbClust;k++){
	//	cout << "Cluster : " << k+1 << endl;
	//	cout << "weight = " << _optimal_weights[k] << endl;
	//	cout << "Mu = ";
	//	for(unsigned int d=0;d<_dim;d++) cout << _optimal_mu(k,d) << " ";
	//	cout << endl;
	//	MatrixXd cov = _optimal_V.block(k*_dim,0,_dim,_dim);
	//	cout << "Cov = " << endl;
	//	cout << cov << endl;
	//}
	//cout << "Time means = " << time_means << endl;
	//cout << "Time weighted cluster prob = " << time_weighted_cluster_prob << endl;
	//cout << "Time weighted prob = " << time_weighted_prob << endl;
	//cout << "Time resp = " << time_resp << endl;
	//cout << "Time variances = " << time_variances << endl;
	//cout << "Time nbdat = " << time_nbdat << endl;
	//cout << "Time kmeans = " << time_kmeans << endl;
}

double GMMTools::expectation(){
	ComputeWeightedClusterLogProb();
	ComputeWeightedLogProb();
	ComputeResp();
	return AverageWeightedLogProb();
}

void GMMTools::maximization(){
	ComputeNbDat2Cluster();
	ComputeWeights();
	ComputeMeans();
	ComputeVariances();
}

void GMMTools::ComputeWeights(){
	_weights = _nbDat2Cluster / (double) _nbDat;
}

void GMMTools::ComputeMeans(){
	//time_beg = high_resolution_clock::now();

	MatrixXd weighted_sum = _Resp * (*(_dataMat));
	for(unsigned int k=0;k<_nbClust;k++) _mu.row(k) = (weighted_sum.row(k).array() / _nbDat2Cluster[k]).matrix();
	
	//time_end = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(time_end - time_beg);
	//time_means += duration.count();
}

void GMMTools::ComputeVariances(){
	//time_beg = high_resolution_clock::now();
	
	for(unsigned int k=0;k<_nbClust;k++){
		MatrixXd centered_data = _dataMat->rowwise() - _mu.row(k);
		VectorXd resp_k = _Resp.row(k);
		MatrixXd weighted_cov = ( (centered_data.transpose() * resp_k.asDiagonal() * centered_data) / _nbDat2Cluster[k] ) + regcovar;
		_V.block(k * _dim, 0, _dim, _dim) = weighted_cov;
		_V_ldlt[k] = weighted_cov.ldlt();
		_det_V[k] = _V_ldlt[k].vectorD().array().log().sum();
	}
	
	//time_end = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(time_end - time_beg);
	//time_variances += duration.count();
}

void GMMTools::ComputeNbDat2Cluster(){
	//time_beg = high_resolution_clock::now();

	_nbDat2Cluster = _Resp.rowwise().sum();

	//time_end = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(time_end - time_beg);
	//time_nbdat += duration.count();
}

void GMMTools::ComputeWeightedClusterLogProb(){ // compute the log of the weighted probability in each cluster for each datapoint
	//time_beg = high_resolution_clock::now();
	
	for(unsigned int k=0;k<_nbClust;k++){
		MatrixXd X_centered = _dataMat->rowwise() - _mu.row(k);
        	MatrixXd temp = _V_ldlt[k].solve(X_centered.transpose());
        	VectorXd dist = (X_centered.array() * temp.transpose().array()).rowwise().sum();
        	_WeightedClusterLogProb.row(k) = logFacNorm - 0.5 * ( _det_V[k] + dist.array());
	}

	//time_end = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(time_end - time_beg);
	//time_weighted_cluster_prob += duration.count();
}

void GMMTools::ComputeWeightedLogProb(){
	//time_beg = high_resolution_clock::now();

	_WeightedLogProb = _WeightedClusterLogProb.array().exp().colwise().sum().log();

	//time_end = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(time_end - time_beg);
	//time_weighted_prob += duration.count();
}

double GMMTools::AverageWeightedLogProb(){
	return _WeightedLogProb.mean();
}

void GMMTools::ComputeResp(){ // compute responsibilities
	//time_beg = high_resolution_clock::now();
	
	for(unsigned int k=0;k<_nbClust;k++){
		_Resp.row(k) = (_WeightedClusterLogProb.row(k).array() - _WeightedLogProb.transpose().array()).exp();
	}
	
	//time_end = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(time_end - time_beg);
	//time_resp += duration.count();
}

double GMMTools::ComputeBIC(){
	_LogLikelihood = expectation();
	unsigned int NbIndepParams = ( ( _nbClust * ( ( _dim * ( _dim + 1.)/2.) + _dim + 1. ) ) - 1. );
	return -( 2.*_nbDat*_LogLikelihood ) + ( log((double) _nbDat) * (double) NbIndepParams ); // BIC
	//return -( 2.*_nbDat*_LogLikelihood ) + ( (double) NbIndepParams ); // AIC
}

void GMMTools::ComputeMLC(MatrixXd &MLC){
	ComputeWeightedClusterLogProb();
	ComputeWeightedLogProb();
	ComputeResp();
	MLC = _Resp;
}

void GMMTools::SaveOptimalRun(){
	_optimal_LogLikelihood = _LogLikelihood;
	_optimal_weights = _weights;
	_optimal_mu = _mu;
	_optimal_V = _V;
	for(unsigned int k=0;k<_nbClust;k++){
		_optimal_det_V[k] = _det_V[k];
		_optimal_V_ldlt[k] = _V_ldlt[k];
	}
}

void GMMTools::SetOptimalRun(){
	_LogLikelihood = _optimal_LogLikelihood;
	_weights = _optimal_weights;
	_mu = _optimal_mu;
	_V = _optimal_V;
	for(unsigned int k=0;k<_nbClust;k++){
		_det_V[k] = _optimal_det_V[k];
		_V_ldlt[k] = _optimal_V_ldlt[k];
	}
}


void GMMTools::Initialize(){
	if( InitMethod == "KMEANS" || InitMethod == "KMEANSPP" ) InitFromKMeans();
	else if( InitMethod == "RANDOM" ) RandomInit();
	else{
		cerr << "The initialization method provided for GMM fitting is unknown, aborting" << endl;
		exit(EXIT_FAILURE);
	}
}

void GMMTools::setSeed(unsigned int &seed){
        FixedSeed = true;
        this->seed = seed;
}

void GMMTools::InitFromKMeans(){
	if( !IsKMeans ){
		_MyKMeansTools = new KMeansTools(_nbClust, _dataMat, _nbDat, _dim);
		IsKMeans = true;
	}
	//time_beg = high_resolution_clock::now();
	_MyKMeansTools->ReadProperties(current_Properties);
	if( FixedSeed ) _MyKMeansTools->setSeed(seed);
	
	if( InitMethod == "KMEANS" ) _MyKMeansTools->fit();
	else if( InitMethod == "KMEANSPP" ){
		_MyKMeansTools->KMeansPPInitialization();
		_MyKMeansTools->AffectData2Cluster();
		_MyKMeansTools->ComputeFullVariances();
	}else{
		cerr << "Error in kmeans initialization of GMM" << endl;
		exit(EXIT_FAILURE);
	}
	//time_end = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(time_end - time_beg);
	//time_kmeans += duration.count();

	MatrixXd tempmat1(_dim,_dim);
	_mu = _MyKMeansTools->getCentroids();
	_V = _MyKMeansTools->getV();
	for(unsigned int k=0;k<_nbClust;k++){
		_weights[k] = (double) _MyKMeansTools->getNbDat2Cluster(k) / (double) _nbDat;
		tempmat1 = _V.block(k * _dim, 0, _dim, _dim);
		tempmat1 += regcovar;
		_V.block(k * _dim, 0, _dim, _dim) = tempmat1;
		_V_ldlt[k] = tempmat1.ldlt();
		_det_V[k] = _V_ldlt[k].vectorD().array().log().sum();

	}

}

void GMMTools::RandomInit(){
	// Initialize randomly by picking random points in the dataset and affecting the same variance to all cluster given by the variance of the whole dataset
// TODO
cerr << "not implemented for the moment" << endl;
exit(EXIT_FAILURE);
	//	vector<unsigned int> already_stored;
//	bool store;
//	unsigned int rand_index;
//	unsigned int ind, ind1, ind2;
//	for(unsigned int k=0;k<_nbClust;k++){
//		weights[k] = 1./_nbClust;
//		store = false;
//		while( !store ){
//			srand(time(0));
//			rand_index = rand() % ( _nbDat + 1 );
//			store = true;
//			for(unsigned int i=0;i<already_stored.size();i++){
//				if( already_stored[i] == rand_index ){
//				       store = false;
//			      		break;	       
//				}
//			}
//		}
//		already_stored.push_back(rand_index);
//		ind = rand_index*_dim;
//		for(unsigned int d=0;d<_dim;d++) mu[k][d] = _datapoints[ind+d];
//	}
//	double *mean = new double[_dim];
//	for(unsigned int d=0;d<_dim;d++) mean[d] = 0.;
//	for(unsigned int i=0;i<_nbDat;i++){
//		ind = i*_dim;
//		for(unsigned int d=0;d<_dim;d++) mean[d] += _datapoints[ind+d];
//	}
//	for(unsigned int d=0;d<_dim;d++) mean[d] /= _nbDat;
//	for(unsigned int k=0;k<_nbClust;k++){
//		for(unsigned int d1=0;d1<_dim;d1++){
//			for(unsigned int d2=d1;d2<_dim;d2++){
//				ind = d1*_dim+d2;
//				V[k][ind] = 0.;
//				for(unsigned int i=0;i<_nbDat;i++) V[k][ind] += ( _datapoints[i*dim+d1] - mean[d1] ) * ( _datapoints[i*dim+d2] - mean[d2] );
//				V[k][ind] /= _nbDat;
//			}
//			for(unsigned int d2=d1+1;d2<dim;d2++) V[k][d2*dim+d1] = V[k][d1*dim+d2];
//		}
//		MT->invMat_LU(V,V_inv,dim,0,nbFilter,filter_value,det_V[filter_value]);
//	}

}



void GMMTools::readFixedParams(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_tol, pos_maxIter, pos_nbInit, pos_InitMethod;
	string buffer_s, line;
	unsigned int ReadOk(0);
	if(file){
		while(file){
			getline(file,line);
			pos_tol=line.find("GMM_TOL_LKH_EM");
			if(pos_tol!=string::npos){
				istringstream text(line);
				text >> buffer_s >> tol_Lkh_EM;
				ReadOk++;
			}
			pos_maxIter=line.find("GMM_MAX_ITER_EM");
			if(pos_maxIter!=string::npos){
				istringstream text(line);
				text >> buffer_s >> MaxIter_EM;
				ReadOk++;
			}
			pos_nbInit=line.find("GMM_NB_INIT");
			if(pos_nbInit!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nbInit;
				ReadOk++;
			}
			pos_InitMethod=line.find("GMM_INIT_METHOD");
			if(pos_InitMethod!=string::npos){
				istringstream text(line);
				text >> buffer_s >> InitMethod;
				ReadOk++;
			}
		}
	}else{
		cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
		exit(EXIT_FAILURE);
	}
	file.close();
	if( ReadOk != 4 ){
		cerr << "Error during reading of FixedParameters.dat for GaussianMixtureModel, aborting" << endl;
		exit(EXIT_FAILURE);
	}
}

void GMMTools::ReadProperties(vector<string> Properties){
	size_t pos_tol, pos_maxIter, pos_nbInit, pos_InitMethod;
	string buffer_s;
	size_t pos_kmtol, pos_kmmaxIter, pos_kmnbInit;
	for(unsigned int i=0;i<Properties.size();i++){
		pos_tol=Properties[i].find("GMM_TOL_LKH_EM");
		if(pos_tol!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> tol_Lkh_EM;
			// search if this properties is already stored in current_Properties and remove it to store the new value
			size_t thispos;
			for(unsigned int p=0;p<current_Properties.size();p++){
				thispos = current_Properties[p].find("GMM_TOL_LKH_EM");
				if( thispos!=string::npos ) current_Properties.erase(current_Properties.begin()+p);
			}
			current_Properties.push_back(Properties[i]);
		}
		pos_maxIter=Properties[i].find("GMM_MAX_ITER_EM");
		if(pos_maxIter!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> MaxIter_EM;
			size_t thispos;
			for(unsigned int p=0;p<current_Properties.size();p++){
				thispos = current_Properties[p].find("GMM_MAX_ITER_EM");
				if( thispos!=string::npos ) current_Properties.erase(current_Properties.begin()+p);
			}
			current_Properties.push_back(Properties[i]);
		}
		pos_nbInit=Properties[i].find("GMM_NB_INIT");
		if(pos_nbInit!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> nbInit;
			size_t thispos;
			for(unsigned int p=0;p<current_Properties.size();p++){
				thispos = current_Properties[p].find("GMM_NB_INIT");
				if( thispos!=string::npos ) current_Properties.erase(current_Properties.begin()+p);
			}
			current_Properties.push_back(Properties[i]);
		}
		pos_InitMethod=Properties[i].find("GMM_INIT_METHOD");
		if(pos_InitMethod!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> InitMethod;
			size_t thispos;
			for(unsigned int p=0;p<current_Properties.size();p++){
				thispos = current_Properties[p].find("GMM_INIT_METHOD");
				if( thispos!=string::npos ) current_Properties.erase(current_Properties.begin()+p);
			}
			current_Properties.push_back(Properties[i]);
		}
		// Read also KMeans properties
		pos_kmtol=Properties[i].find("KMEANS_TOL");
		if(pos_kmtol!=string::npos){
			size_t thispos;
			for(unsigned int p=0;p<current_Properties.size();p++){
				thispos = current_Properties[p].find("KMEANS_TOL");
				if( thispos!=string::npos ) current_Properties.erase(current_Properties.begin()+p);
			}
			current_Properties.push_back(Properties[i]);
		}
		pos_kmmaxIter=Properties[i].find("KMEANS_MAX_ITER");
		if(pos_kmmaxIter!=string::npos){
			size_t thispos;
			for(unsigned int p=0;p<current_Properties.size();p++){
				thispos = current_Properties[p].find("KMEANS_MAX_ITER");
				if( thispos!=string::npos ) current_Properties.erase(current_Properties.begin()+p);
			}
			current_Properties.push_back(Properties[i]);
		}
		pos_kmnbInit=Properties[i].find("KMEANS_NB_INIT");
		if(pos_kmnbInit!=string::npos){
			size_t thispos;
			for(unsigned int p=0;p<current_Properties.size();p++){
				thispos = current_Properties[p].find("KMEANS_NB_INIT");
				if( thispos!=string::npos ) current_Properties.erase(current_Properties.begin()+p);
			}
			current_Properties.push_back(Properties[i]);
		}
	}

	if( ( InitMethod == "KMEANS" || InitMethod == "KMEANSPP" ) && !IsKMeans ){
		_MyKMeansTools = new KMeansTools(_nbClust, _dataMat, _nbDat, _dim);
		IsKMeans = true;
	}

}

GMMTools::~GMMTools(){
	if( IsKMeans ) delete _MyKMeansTools;
}
