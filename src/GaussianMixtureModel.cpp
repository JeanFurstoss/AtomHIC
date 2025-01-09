#include "GaussianMixtureModel.h"
#include "AtomHicConfig.h"
#include <cmath>
#include <filesystem>
#include <dirent.h>

using namespace std; 

GaussianMixtureModel::GaussianMixtureModel(){
	this->name = "GaussianMixtureModel"; 
	readFixedParams();
}

void GaussianMixtureModel::setDescriptors(Descriptors *D){
	MachineLearningModel::setDescriptors(D);

	if( !IsRead ){
		weights = new double[nbMaxClusters*nbFilter];
		mu = new double[dim*nbMaxClusters*nbFilter];
		V_inv = new long double[dim2*nbMaxClusters*nbFilter];
		det_V = new long double[nbMaxClusters*nbFilter];
		nbClust = new unsigned int[nbFilter];
	}
	
	if( IsFilterIndexModified ) ChangeFilterIndex();
	
	if( this->IsDescriptor ){
		delete[] V;
		delete[] D_i;
		delete[] C_di;
		delete[] buffer_di;
		delete[] E_d; 
		delete[] weights_old;
		delete[] mu_old;
		delete[] V_old;
		delete[] V_inv_old;
		delete[] det_V_old;
		delete[] BIC;
		delete[] LogLikelihood;
	}
	this->IsDescriptor = true;
	V = new long double[dim2*nbMaxClusters*nbFilter];
	D_i = new double[nbDatMax];
	C_di = new double[nbMaxClusters*nbDatMax];
	buffer_di = new double[nbMaxClusters*nbDatMax];
	E_d = new double[nbMaxClusters];
	
	weights_old = new double[nbMaxClusters*nbFilter];
	mu_old = new double[dim*nbMaxClusters*nbFilter];
	V_old = new long double[dim2*nbMaxClusters*nbFilter];
	V_inv_old = new long double[dim2*nbMaxClusters*nbFilter];
	det_V_old = new long double[nbMaxClusters*nbFilter];

	BIC = new double[nbFilter];
	LogLikelihood = new long double[nbFilter];
}

void GaussianMixtureModel::fitOptimalGMM(unsigned int &_nbClust_min, unsigned int &_nbClust_max){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, training aborted" << endl;
		exit(EXIT_FAILURE);
	}
	if( _nbClust_max > nbMaxClusters ){
		cerr << "The number of cluster max for GMM training is higher than the maximum number of cluster (see /data/FixedParameters/FixedParameters.dat to modify), aborting" << endl;
		exit(EXIT_FAILURE);
	}
	
	unsigned int nbRun = _nbClust_max-_nbClust_min+1;	
	if( nbRun < 3 ){ //TODO allow for less than 3 clusters
		cerr << "The training of GMM for an unknown N number of cluster is not possible with N < 3, aborting" << endl;
		exit(EXIT_FAILURE);
	}

	cout << "Training the GMM for number of cluster between " << _nbClust_min << " and " << _nbClust_max << endl;

	if( !SavedVariablesInitialized ){
		saved_bic = new double[nbRun*nbFilter];
		saved_nbClust = new unsigned int[nbRun*nbFilter];
		saved_weights = new double[nbMaxClusters*nbRun*nbFilter];
		saved_mu = new double[dim*nbMaxClusters*nbRun*nbFilter];
		saved_V = new long double[dim2*nbMaxClusters*nbRun*nbFilter];
		saved_V_inv = new long double[dim2*nbMaxClusters*nbRun*nbFilter];
		saved_det_V = new long double[nbMaxClusters*nbRun*nbFilter];
		SavedVariablesInitialized = true;
	}

	for(unsigned int f=0;f<nbFilter;f++){
		cout << "Training for filter " << _MyDescriptors->getFilterValue(f) << endl;
		cout << "Training for = " << endl;
		cout << _nbClust_min << " clusters" << endl;

		TrainModel(_nbClust_min,f);

		cout << "Training done !" << endl;
		unsigned int current_save_index = 0;
		SaveVariables(current_save_index,f);
		
		bool OptModelFound = false;
		for(unsigned int k=_nbClust_min+1;k<_nbClust_max+1;k++){
			cout << k << " clusters" << endl;
			current_save_index++;
			TrainModel(k,f);
			cout << "Training done !" << endl;
			SaveVariables(current_save_index,f);
			if( current_save_index >= nb_bic_increase ){
				OptModelFound = true;
				for(unsigned int b=0;b<nb_bic_increase;b++){
					if( saved_bic[(current_save_index-b)*nbFilter+f] < saved_bic[(current_save_index-nb_bic_increase)*nbFilter+f] ){
						OptModelFound = false;
						break;
					}
				}
				if( OptModelFound ){
				       cout << "We find a minimum in BIC, we then stop the computation and consider that the optimal GMM is for a numer of cluster = " << k-nb_bic_increase << endl;
				       OptModelFound = true;
				       unsigned int opt = current_save_index-nb_bic_increase;
				       SetOptimalModel(opt,f);
				       break;
				}
			}
		}
		if( !OptModelFound ){
			cout << "We didn't find a local minimum of the BIC reflecting an optimal GMM" << endl;
			// Elbow method, we search the intersection between the two lines formed by the two first and the two last moints (maybe refine later)
			double diff_beg = saved_bic[nbFilter+f]-saved_bic[f];
			double diff_end = saved_bic[(nbRun-1)*nbFilter+f]-saved_bic[(nbRun-2)*nbFilter+f];
			unsigned int opt;
			if( diff_end < diff_beg*fac_elbow ){
				cout << "The BIC vs number of cluster curve does not really present an elbow" << endl;
				if( after_elbow_choice == "Max" ){
					cout << "We then selected N = number_of_cluster_max according to /data/FixedParameters/FixedParameters.dat" << endl;
					opt = nbRun-1;
				}else if( after_elbow_choice == "Min" ){
					cout << "We then selected N = number_of_cluster_min according to /data/FixedParameters/FixedParameters.dat" << endl;
					opt = 0.;
				}else{
					unsigned int nbclust_to_choose;
					istringstream iss_N(after_elbow_choice);
					iss_N >> nbclust_to_choose;
					if( nbclust_to_choose < _nbClust_min ){
						cout << "The provided number of cluster by GMM_NO_BIC_MIN_NO_ELBOW_CHOICE in /data/FixedParameters/FixedParameters.dat is lower than nbClust_min, we then selected N = nbClust_max" << endl;
						opt = nbRun-1;
					}else if( nbclust_to_choose > _nbClust_max ){
						cout << "The provided number of cluster by GMM_NO_BIC_MIN_NO_ELBOW_CHOICE in /data/FixedParameters/FixedParameters.dat is higher than nbClust_max, we then selected N = nbClust_max" << endl;
						opt = nbRun-1;
					}else{
						cout << "We then selected N = " << nbclust_to_choose << " according to /data/FixedParameters/FixedParameters.dat" << endl;
						for(unsigned int nsaved=0;nsaved<nbRun;nsaved++){
							if( saved_nbClust[nsaved*nbFilter+f] == nbclust_to_choose ){
								opt = nsaved;
								break;
							}
						}
					}
				}
				cout << "Consider increasing number_of_cluster_max for GMM this training" << endl;
			}else{
				unsigned int opt_N = round( ( saved_bic[(nbRun-1)*nbFilter+f] - saved_bic[f] + ( (double) saved_nbClust[f] * diff_beg ) - ( (double) saved_nbClust[(nbRun-1)*nbFilter+f] * diff_end ) ) / ( diff_beg - diff_end ) );
				if( opt_N < 1 ) opt_N = 1;
				for(unsigned int i=0;i<nbRun;i++){
					if( saved_nbClust[i*nbFilter+f] == opt_N ){
						opt = i;
						break;
					}
				}
				if( opt > _nbClust_max ){
					opt = nbRun-1;
					cout << "The Elbow method as implemented here does not found an optimal N number of cluster, we then set N = number_of_cluster_max" << endl;
				}else{
					cout << "From the Elbow method based on BIC, we infer an optimal number of cluster = " << opt_N << endl;
				}
			}
			SetOptimalModel(opt,f);
		}
		ofstream writefile("BIC_vs_N_"+_MyDescriptors->getFilterValue(f)+".dat");
		writefile << "Number_of_Cluster BIC" << endl;
		for(unsigned int i=0;i<current_save_index+1;i++) writefile << saved_nbClust[i*nbFilter+f] << " " << saved_bic[i*nbFilter+f] << endl;
		writefile.close();
		cout << "Yo may give a look on the BIC_vs_N_" << _MyDescriptors->getFilterValue(f) << ".dat file to see if the optimal number of cluster suites to you" << endl;
	}
}

void GaussianMixtureModel::TrainModel(unsigned int &_nbClust, unsigned int &filter_value){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, training aborted" << endl;
		exit(EXIT_FAILURE);
	}
	if( _nbClust > nbMaxClusters ){
		cerr << "The number of cluster for GMM training is higher than the maximum number of cluster (see /data/FixedParameters/FixedParameters.dat to modify), aborting" << endl;
		exit(EXIT_FAILURE);
	}

	if( !TrainSavedVariablesInitialized ){
		train_saved_bic = new double[nbInit*nbFilter];
		train_saved_weights = new double[nbMaxClusters*nbInit*nbFilter];
		train_saved_mu = new double[dim*nbMaxClusters*nbInit*nbFilter];
		train_saved_V = new long double[dim2*nbMaxClusters*nbInit*nbFilter];
		train_saved_V_inv = new long double[dim2*nbMaxClusters*nbInit*nbFilter];
		train_saved_det_V = new long double[nbMaxClusters*nbInit*nbFilter];
		TrainSavedVariablesInitialized = true;
	}
	nbClust[filter_value] = _nbClust;
	for(unsigned int i=0;i<nbInit;i++){	
		if( InitMethod == "KMEANS" ) InitFromKMeans(nbClust[filter_value],filter_value);
		else if( InitMethod == "KMEANSPP" ) InitFromKMeansPP(nbClust[filter_value],filter_value);
		else if( InitMethod == "RANDOM" ) RandomInit(nbClust[filter_value],filter_value);
		else{
			cerr << "The initialization method provided for GMM fitting is unknown, aborting" << endl;
			exit(EXIT_FAILURE);
		}
		UpdateParams(filter_value);
		ComputeLogLikelihood(filter_value);
		long double Lkh_old = LogLikelihood[filter_value];
		double eps;
		unsigned int iter = 0;
		do{
			EM(filter_value);
			UpdateParams(filter_value);
			ComputeLogLikelihood(filter_value);
			eps = LogLikelihood[filter_value] - Lkh_old;
			Lkh_old = LogLikelihood[filter_value];
			cout << "LogLkh = << " << Lkh_old << endl;
			iter++;
		}while( eps > tol_Lkh_EM && iter < MaxIter_EM );
		
		ComputeBIC(filter_value);
		SaveTrainVariables(i,filter_value);
		if( eps > tol_Lkh_EM ){
			cout << "The EM algorythm does not converge until tolerance for a number of cluster = " << nbClust[filter_value] << ", final residual = " << eps << endl;
			cout << "Maybe increase number of iteration or the tolerance for EM in /data/FixedParameters/FixedParameters.dat" << endl;
		}
	}
	// at the end we keep the one with the lowest BIC
	cout << "Saved BIC:" << endl;
	for(unsigned int i=0;i<nbInit;i++) cout << train_saved_bic[i*nbFilter+filter_value] << " ";
	cout << endl;
	unsigned int index_opt = MT->min_p_ind(train_saved_bic,nbInit,nbFilter,filter_value);
	cout << "MIN = " << train_saved_bic[index_opt*nbFilter+filter_value] << endl;
	SetOptimalTrainModel(index_opt,filter_value);
}

void GaussianMixtureModel::SetOptimalTrainModel(unsigned int &opt_index, unsigned int &filter_value){
	BIC[filter_value] = train_saved_bic[opt_index*nbFilter+filter_value];
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		unsigned int ind1 = k*nbFilter+filter_value;
		unsigned int ind2 = opt_index*nbMaxClusters*nbFilter+k*nbFilter+filter_value;
		weights[ind1] = train_saved_weights[ind2];
		det_V[ind1] = train_saved_det_V[ind2];
		for(unsigned int d1=0;d1<dim;d1++){
			mu[k*dim*nbFilter+d1*nbFilter+filter_value] = train_saved_mu[opt_index*nbMaxClusters*dim*nbFilter+k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				unsigned int indb1 = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				unsigned int indb2 = opt_index*nbMaxClusters*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				V[indb1] = train_saved_V[indb2];
				V_inv[indb1] = train_saved_V_inv[indb2];
			}
		}
	}
}

void GaussianMixtureModel::SetOptimalModel(unsigned int &opt_index, unsigned int &filter_value){
	BIC[filter_value] = saved_bic[opt_index*nbFilter+filter_value];
	nbClust[filter_value] = saved_nbClust[opt_index*nbFilter+filter_value];
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		unsigned int ind1 = k*nbFilter+filter_value;
		unsigned int ind2 = opt_index*nbMaxClusters*nbFilter+k*nbFilter+filter_value;
		weights[ind1] = saved_weights[ind2];
		det_V[ind1] = saved_det_V[ind2];
		for(unsigned int d1=0;d1<dim;d1++){
			mu[k*dim*nbFilter+d1*nbFilter+filter_value] = saved_mu[opt_index*nbMaxClusters*dim*nbFilter+k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				unsigned int indb1 = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				unsigned int indb2 = opt_index*nbMaxClusters*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				V[indb1] = saved_V[indb2];
				V_inv[indb1] = saved_V_inv[indb2];
			}
		}
	}
}

void GaussianMixtureModel::SaveTrainVariables(unsigned int &current, unsigned int &filter_value){
	if( !TrainSavedVariablesInitialized ){
		cerr << "Saved variables for GMM optimization have not been initialized, something went wrong, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	if( BIC[filter_value] == BIC[filter_value] ) train_saved_bic[current*nbFilter+filter_value] = BIC[filter_value];
	else train_saved_bic[current*nbFilter+filter_value] = 0.;
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		unsigned int ind1 = current*nbMaxClusters*nbFilter+k*nbFilter+filter_value;
		unsigned int ind2 = k*nbFilter+filter_value;
		train_saved_weights[ind1] = weights[ind2];
		train_saved_det_V[ind1] = det_V[ind2];
		for(unsigned int d1=0;d1<dim;d1++){
			train_saved_mu[current*nbMaxClusters*dim*nbFilter+k*dim*nbFilter+d1*nbFilter+filter_value] = mu[k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				unsigned int indb1 = current*nbMaxClusters*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				unsigned int indb2 = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				train_saved_V[indb1] = V[indb2];
				train_saved_V_inv[indb1] = V_inv[indb2];
			}
		}
	}
}

void GaussianMixtureModel::SaveVariables(unsigned int &current, unsigned int &filter_value){
	if( !SavedVariablesInitialized ){
		cerr << "Saved variables for GMM optimization have not been initialized, something went wrong, aborting" << endl;
		exit(EXIT_FAILURE);
	}

	saved_bic[current*nbFilter+filter_value] = BIC[filter_value];
	saved_nbClust[current*nbFilter+filter_value] = nbClust[filter_value];
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		unsigned int ind1 = current*nbMaxClusters*nbFilter+k*nbFilter+filter_value;
		unsigned int ind2 = k*nbFilter+filter_value;
		saved_weights[ind1] = weights[ind2];
		saved_det_V[ind1] = det_V[ind2];
		for(unsigned int d1=0;d1<dim;d1++){
			saved_mu[current*nbMaxClusters*dim*nbFilter+k*dim*nbFilter+d1*nbFilter+filter_value] = mu[k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				unsigned int indb1 = current*nbMaxClusters*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				unsigned int indb2 = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				saved_V[indb1] = V[indb2];
				saved_V_inv[indb1] = V_inv[indb2];
			}
		}
	}
}

void GaussianMixtureModel::RandomInit(unsigned int &_nbClust, unsigned int &filter_value){
	// Initialize randomly by picking random points in the dataset and affecting the same variance to all cluster given by the variance of the whole dataset
	vector<unsigned int> already_stored;
	bool store;
	unsigned int rand_index;
	if( _nbClust > nbDat[filter_value] ){
		cerr << "The number of cluster s higher than the number of data, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	for(unsigned int k=0;k<_nbClust;k++){
		weights[k*nbFilter+filter_value] = 1./_nbClust;
		store = false;
		while( !store ){
			srand(time(0));
			rand_index = rand() % ( nbDat[filter_value] + 1 );
			store = true;
			for(unsigned int i=0;i<already_stored.size();i++){
				if( already_stored[i] == rand_index ){
				       store = false;
			      		break;	       
				}
			}
		}
		already_stored.push_back(rand_index);
		for(unsigned int d=0;d<dim;d++) mu[k*dim*nbFilter+d*nbFilter+filter_value] = _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+rand_index)*dim+d];
	}
	double *mean = new double[dim];
	for(unsigned int d=0;d<dim;d++) mean[d] = 0.;
	for(unsigned int i=0;i<nbDat[filter_value];i++){	
		for(unsigned int d=0;d<dim;d++) mean[d] += _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+i)*dim+d];
	}
	for(unsigned int d=0;d<dim;d++) mean[d] /= nbDat[filter_value];
	for(unsigned int d1=0;d1<dim;d1++){
		for(unsigned int d2=d1;d2<dim;d2++){
			unsigned int ind = 0*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
			V[ind] = 0.;
			for(unsigned int i=0;i<nbDat[filter_value];i++) V[ind] += ( _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+i)*dim+d1] - mean[d1] ) * ( _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+i)*dim+d2] - mean[d2] );
			V[ind] /= nbDat[filter_value];
		}
		for(unsigned int d2=d1+1;d2<dim;d2++) V[0*dim2*nbFilter+d2*dim*nbFilter+d1*nbFilter+filter_value] = V[0*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
	}
	MT->invMat_LU(V,V_inv,dim,0,nbFilter,filter_value,det_V[filter_value]);
	for(unsigned int k=1;k<_nbClust;k++){
		det_V[k*nbFilter+filter_value] = det_V[filter_value];
		for(unsigned int d1=0;d1<dim;d1++){
			for(unsigned int d2=0;d2<dim;d2++){
				unsigned int ind = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				unsigned int ind0 = 0*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				V[ind] = V[ind0];
				V_inv[ind] = V_inv[ind0];
			}
		}
	}
	delete[] mean;
}


void GaussianMixtureModel::InitFromKMeansPP(unsigned int &_nbClust, unsigned int &filter_value){
	if( !IsKMeans ){
		MyKM = new KMeans();
		if( IsKMeansProperties ) MyKM->ReadProperties(KMeansProperties);
		MyKM->setDescriptors(_MyDescriptors);
	}

	MyKM->KMeansPPInitialization(nbClust[filter_value],filter_value);
	MyKM->AffectData2Cluster(filter_value);
	MyKM->ComputeFullVariances(filter_value);
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		unsigned int ind1 = k*nbFilter+filter_value;
		weights[ind1] = (double) MyKM->getNbDat2Cluster()[ind1] / (double) nbDat[filter_value];
		det_V[ind1] = MyKM->getDet_V()[ind1];
		for(unsigned int d1=0;d1<dim;d1++){
			unsigned int ind2 = k*dim*nbFilter+d1*nbFilter+filter_value;
			mu[ind2] = MyKM->getCentroids()[ind2];
			for(unsigned int d2=0;d2<dim;d2++){
				unsigned int ind3 = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				V[ind3] = MyKM->getV()[ind3];
				V_inv[ind3] = MyKM->getV_inv()[ind3];
			}
		}
	}

}

void GaussianMixtureModel::InitFromKMeans(unsigned int &_nbClust, unsigned int &filter_value){
	if( !IsKMeans ){
		MyKM = new KMeans();
		if( IsKMeansProperties ) MyKM->ReadProperties(KMeansProperties);
		MyKM->setDescriptors(_MyDescriptors);
	}

	MyKM->TrainModel(nbClust[filter_value],filter_value);
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		unsigned int ind1 = k*nbFilter+filter_value;
		weights[ind1] = (double) MyKM->getNbDat2Cluster()[ind1] / (double) nbDat[filter_value];
		det_V[ind1] = MyKM->getDet_V()[ind1];
		for(unsigned int d1=0;d1<dim;d1++){
			unsigned int ind2 = k*dim*nbFilter+d1*nbFilter+filter_value;
			mu[ind2] = MyKM->getCentroids()[ind2];
			for(unsigned int d2=0;d2<dim;d2++){
				unsigned int ind3 = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value;
				V[ind3] = MyKM->getV()[ind3];
				V_inv[ind3] = MyKM->getV_inv()[ind3];
			}
		}
	}
}

void GaussianMixtureModel::UpdateParams(unsigned int &filter_value){
	for(unsigned int d=0;d<nbClust[filter_value];d++){
		unsigned int ind1 = d*nbFilter+filter_value;
		weights_old[ind1] = weights[ind1];
		det_V_old[ind1] = det_V[ind1];
		for(unsigned int x1=0;x1<dim;x1++){
			unsigned int ind2 = d*dim*nbFilter+x1*nbFilter+filter_value;
			mu_old[ind2] = mu[ind2];
			for(unsigned int x2=0;x2<dim;x2++){
				unsigned int ind3 = d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value;
				V_old[ind3] = V[ind3];
				V_inv_old[ind3] = V_inv[ind3];
			}
		}
	}
}

void GaussianMixtureModel::EM(unsigned int &filter_value){
	// Expectation
	for(unsigned int d=0;d<nbClust[filter_value];d++) E_d[d] = 0.;
	for(unsigned int i=0;i<nbDat[filter_value];i++){
		D_i[i] = 0.;
		for(unsigned int d=0;d<nbClust[filter_value];d++){
			buffer_di[(d*nbDat[filter_value])+i] = Prob_Cluster(d,i,filter_value);
			D_i[i] += buffer_di[(d*nbDat[filter_value])+i];
		}
		for(unsigned int d=0;d<nbClust[filter_value];d++) C_di[(d*nbDat[filter_value])+i] = buffer_di[(d*nbDat[filter_value])+i] / D_i[i];
		for(unsigned int d=0;d<nbClust[filter_value];d++) E_d[d] += C_di[(d*nbDat[filter_value])+i];
	}
	// Maximization
	#pragma omp parallel for
	for(unsigned int d=0;d<nbClust[filter_value];d++){
		// new weights
		weights[d*nbFilter+filter_value] = E_d[d]/nbDat[filter_value];
		// new means
		for(unsigned int x=0;x<dim;x++){
			unsigned int ind = d*dim*nbFilter+x*nbFilter+filter_value;
			mu[ind] = 0.;
			for(unsigned int i=0;i<nbDat[filter_value];i++) mu[ind] += C_di[(d*nbDat[filter_value])+i]*_MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+i)*dim+x];
			mu[ind] /= E_d[d];
		}
		// new variances
		for(unsigned int x1=0;x1<dim;x1++){
			for(unsigned int x2=x1;x2<dim;x2++){
				unsigned int ind = d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value;
				unsigned int ind1 = d*dim*nbFilter+x1*nbFilter+filter_value;
				unsigned int ind2 = d*dim*nbFilter+x2*nbFilter+filter_value;
				V[ind] = 0.;
				for(unsigned int i=0;i<nbDat[filter_value];i++) V[ind] += C_di[(d*nbDat[filter_value])+i]*( (_MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+i)*dim+x1] - mu[ind1]) * (_MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+i)*dim+x2] - mu[ind2]) );
				V[ind] /= E_d[d];
			}
			// symmetric part
			for(unsigned int x2=x1+1;x2<dim;x2++){
				V[d*dim2*nbFilter+x2*dim*nbFilter+x1*nbFilter+filter_value] = V[d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value];
			}
		}
	}
	for(unsigned int d=0;d<nbClust[filter_value];d++){
		// invert variances
		MT->invMat_LU(V,V_inv,dim,d,nbFilter,filter_value,det_V[d*nbFilter+filter_value]);
	}
}

void GaussianMixtureModel::ComputeLogLikelihood(unsigned int &filter_value){
	LogLikelihood[filter_value] = 0.;
	unsigned int zero_val = 0;
	for(unsigned int d=0;d<nbClust[filter_value];d++) LogLikelihood[filter_value] += Prob_Cluster(d,zero_val,filter_value);
	LogLikelihood[filter_value] = log(LogLikelihood[filter_value]);
	for(unsigned int i=1;i<nbDat[filter_value];i++){
		long double sum = 0.;
		for(unsigned int d=0;d<nbClust[filter_value];d++) sum += Prob_Cluster(d,i,filter_value);
		LogLikelihood[filter_value] += log(sum);
	}
	// Test normalization TODO
	LogLikelihood[filter_value] /= nbDat[filter_value];

}

void GaussianMixtureModel::ComputeBIC(unsigned int &filter_value){
	unsigned int NbIndepParams = ( ( nbClust[filter_value] * ( ( dim*( dim + 1)/2) + dim + 1. ) ) - 1. );
	BIC[filter_value] = -( 2.*LogLikelihood[filter_value] ) + ( log((double) nbDat[filter_value]) * (double) NbIndepParams );
}

// Label the GMM by affecting to each cluster the label which has the highest probability in its and return the second highest value to see if there is overlapping between labels 
void GaussianMixtureModel::Labelling(){
	nbLabel = _MyDescriptors->getNbLabels();
	if( nbLabel < 2 ){
		cerr << "The descriptors are not labelled, or the label has only one value, we then cannot label the GMM, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	unsigned int zero_lab = 0;
	for(unsigned int f=0;f<nbFilter;f++){
		unsigned int max = _MyDescriptors->getLabelsSize(f,zero_lab);
		unsigned int min = _MyDescriptors->getLabelsSize(f,zero_lab);
		for(unsigned int l=1;l<nbLabel;l++){
			 if( _MyDescriptors->getLabelsSize(f,l) > max ) max = _MyDescriptors->getLabelsSize(f,l);
			 else if( _MyDescriptors->getLabelsSize(f,l) < min ) min = _MyDescriptors->getLabelsSize(f,l);
		}
		if( max*tolLabelSize > min ) cout << "Warning the \"" << _MyDescriptors->getFilterValue(f) << "\" descriptors are not homogeneously distributed within the labels (high difference between number of descriptors in each labels), which could biased the results obtained with the fitted GMM" << endl;
	}

	if( !IsLabelled ){
		ClusterLabel = new unsigned int[nbFilter*nbMaxClusters];
		AveLabelProb = new long double[nbFilter*nbMaxClusters*nbLabel];
		IsLabelled = true;
	}

	// Compute average probability of labelled descriptors in each clusters
	double *SecondLabelProb = new double[nbFilter*nbMaxClusters];
	
	for(unsigned int ind=0;ind<nbFilter*nbLabel*nbMaxClusters;ind++) AveLabelProb[ind] = 0.;
	
	for(unsigned int f=0;f<nbFilter;f++){
		for(unsigned int i=0;i<nbDat[f];i++){
			//for(unsigned int k=0;k<nbClust[f];k++) AveLabelProb[k*nbFilter*nbLabel+_MyDescriptors->getLabels_uint()[f*nbDatMax+i]*nbFilter+f] += Prob_Cluster(k,i,f);
			for(unsigned int k=0;k<nbClust[f];k++) AveLabelProb[k*nbFilter*nbLabel+_MyDescriptors->getLabels_uint()[f*nbDatMax+i]*nbFilter+f] += MaximumLikelihoodClassifier(k,i,f); // TODO warning with change in the MLC function
		}
	}
	
	for(unsigned int f=0;f<nbFilter;f++){
		for(unsigned int k=0;k<nbClust[f];k++){
			for(unsigned int l=0;l<nbLabel;l++) AveLabelProb[k*nbFilter*nbLabel+l*nbFilter+f] /= _MyDescriptors->getLabelsSize(f,l);
		}
	}

	// Search which label has the highest and second highest probability in each cluster
	for(unsigned int f=0;f<nbFilter;f++){
		for(unsigned int k=0;k<nbClust[f];k++){
			long double max1 = AveLabelProb[k*nbFilter*nbLabel+f];
			unsigned int ind1 = 0;
			for(unsigned int l=1;l<nbLabel;l++){
				if( AveLabelProb[k*nbFilter*nbLabel+l*nbFilter+f] > max1 ){
					max1 = AveLabelProb[k*nbFilter*nbLabel+l*nbFilter+f];
					ind1 = l;
				}
			}
			unsigned int ind2;
			for(unsigned int l=0;l<nbLabel;l++){
				if( l != ind1 ){
					ind2 = l;
					break;
				}
			}
			long double max2 = AveLabelProb[k*nbFilter*nbLabel+ind2*nbFilter+f];
			for(unsigned int l=0;l<nbLabel;l++){
				if( l != ind1 && l != ind2 && AveLabelProb[k*nbFilter*nbLabel+l*nbFilter+f] > max2 ){
					max2 = AveLabelProb[k*nbFilter*nbLabel+l*nbFilter+f];
					ind2 = l;
				}
			}
			ClusterLabel[k*nbFilter+f] = ind1;
			SecondLabelProb[k*nbFilter+f] = max2;
		}
	}

	for(unsigned int f=0;f<nbFilter;f++){
		double sum2ndProb = 0.;
		for(unsigned int k=0;k<nbClust[f];k++) sum2ndProb += SecondLabelProb[k*nbFilter+f];
		if( sum2ndProb > tol2ndProb ){
			cout << "Warning, the sum of the second highest probability associated to each cluster with filter \"" << _MyDescriptors->getFilterValue(f) << "\" indicate that some labels are not well separated within descriptor space" << endl;
			cout << "Sum of second average probability over all clusters = " << sum2ndProb << endl;
			cout << "You should have a look on the \"labelling.out\" file to see the precise result of the labeling" << endl;
			cout << "We have anyway labelled the GMM with the labels having the highest probability in each cluster" << endl;
		}else cout << "The labelling of the GMM with filter \"" << _MyDescriptors->getFilterValue(f) << "\" have been successfully achieved" << endl;
	}
	
	delete[] SecondLabelProb;

}

long double GaussianMixtureModel::Prob_Cluster(unsigned int &index_cluster, unsigned int &DescriptorIndex, unsigned int &filter_value){
	double sp = 0.;
	for(unsigned int j=0;j<dim;j++){
		unsigned int index = _MyDescriptors->getFilterIndex(filter_value*nbDatMax+DescriptorIndex);
		buffer_vec_2_dim[j] = (_MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+DescriptorIndex)*dim+j]-mu[index_cluster*dim*nbFilter+j*nbFilter+filter_value]);
	}
	for(unsigned int i=0;i<dim;i++){
		buffer_vec_1_dim[i] = 0.;
		for(unsigned int j=0;j<dim;j++) buffer_vec_1_dim[i] += V_inv[index_cluster*dim2*nbFilter+i*dim*nbFilter+j*nbFilter+filter_value]*buffer_vec_2_dim[j];
		sp += buffer_vec_1_dim[i]*buffer_vec_2_dim[i];
	}
	return ( weights[index_cluster*nbFilter+filter_value] / ( pow(2.*M_PI, (double) dim/2.) * sqrt(det_V[index_cluster*nbFilter+filter_value]) ) ) * exp( -.5*sp );
}

double GaussianMixtureModel::MaximumLikelihoodClassifier(unsigned int &index_cluster, unsigned int &DescriptorIndex, unsigned int &filter_value){
	long double sumProbs = 0.;
	for(unsigned int k=0;k<nbClust[filter_value];k++) sumProbs += Prob_Cluster(k,DescriptorIndex,filter_value);
	return Prob_Cluster(index_cluster,DescriptorIndex,filter_value) / sumProbs;
}

// Classify the data using the maximum likelihood classifier
void GaussianMixtureModel::LabelClassification(){
	if( !IsLabelled && !IsRead ){
		cerr << "The GMM is not labelled, we then cannot classify the data, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	if( !IsDescriptor ){
		cerr << "We dont have descriptor to classify, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Classifying the descriptors.." << endl;
	unsigned int nbDatTot = 0;
	for(unsigned int f=0;f<nbFilter;f++) nbDatTot += nbDat[f];
	if( !IsClassified ){
		IsClassified = true;
		Classificator = new double[2*nbDatTot];
	}
	
	long double *LabelProb = new long double[nbLabel];
	for(unsigned int f=0;f<nbFilter;f++){
		for(unsigned int j=0;j<nbDat[f];j++){
			double sum = 0.;
			for(unsigned int l=0;l<nbLabel;l++) LabelProb[l] = 0.;
			for(unsigned int k=0;k<nbClust[f];k++){
				LabelProb[ClusterLabel[k*nbFilter+f]] += Prob_Cluster(k,j,f);
			}
			for(unsigned int l=0;l<nbLabel;l++) sum += LabelProb[l];
			if( sum != 0. ){
				unsigned int index_maxp = MT->max_p_ind(LabelProb,nbLabel);
				Classificator[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*2] = index_maxp;
				Classificator[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*2+1] = LabelProb[index_maxp] / sum;
			}else{
				Classificator[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*2] = nbLabel; 
				Classificator[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*2+1] = 0.;
			}	
		}
	}

        // write the StructureIndex.txt file
        ofstream writefile("StructureIndex.txt");
        writefile << "Structures index used for the Gaussian Mixture Model structural analysis" << endl;
        for(unsigned int l=0;l<nbLabel;l++) writefile << l << " " << Labels[l] << endl;
        writefile << nbLabel << " Not identified" << endl;
        writefile.close();

	delete[] LabelProb;
	cout << "Done" << endl;
	cout << "The StructureIndex.txt file containing the names of the labels has been printed" << endl;
}

// Classify the data using the maximum likelihood classifier
void GaussianMixtureModel::Classify(){
	if( !IsDescriptor ){
		cerr << "We dont have descriptor to classify, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Classifying the descriptors.." << endl;
	unsigned int nbDatTot = 0;
	for(unsigned int f=0;f<nbFilter;f++) nbDatTot += nbDat[f];
	if( !IsClassified ){
		IsClassified = true;
		Classificator = new double[2*nbDatTot];
	}
	
	for(unsigned int f=0;f<nbFilter;f++){
		long double *ClusterProb = new long double[nbClust[f]];
		for(unsigned int j=0;j<nbDat[f];j++){
			double sum = 0.;
			for(unsigned int k=0;k<nbClust[f];k++){
				ClusterProb[k] = Prob_Cluster(k,j,f);
				sum += ClusterProb[k];
			}
			if( sum != 0. ){
				unsigned int index_maxp = MT->max_p_ind(ClusterProb,nbClust[f]);
				Classificator[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*2] = index_maxp;
				Classificator[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*2+1] = ClusterProb[index_maxp] / sum;
			}else{
				Classificator[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*2] = nbClust[f]; 
				Classificator[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*2+1] = 0.;
			}	
		}
		delete[] ClusterProb;
	}
}

void GaussianMixtureModel::ChangeFilterIndex(){
	unsigned int *nbClust_tmp = new unsigned int[nbFilter];
	string *FilterValue_tmp = new string[nbFilter];
	double *weights_tmp = new double[nbMaxClusters*nbFilter];
	long double *det_V_tmp = new long double[nbMaxClusters*nbFilter];
	unsigned int *ClusterLabel_tmp = new unsigned int[nbMaxClusters*nbFilter];
	double *mu_tmp = new double[dim*nbMaxClusters*nbFilter];
	long double *V_inv_tmp = new long double[dim2*nbMaxClusters*nbFilter];

	// Copy data and reinitialize them
	for(unsigned int f=0;f<nbFilter;f++){
		nbClust_tmp[f] = nbClust[f];
		nbClust[f] = 0;
		FilterValue_tmp[f] = FilterValue[f];
		FilterValue[f] = 0.;
		for(unsigned int k=0;k<nbClust_tmp[f];k++){
			weights_tmp[k*nbFilter+f] = weights[k*nbFilter+f];
			weights[k*nbFilter+f] = 0.;
			det_V_tmp[k*nbFilter+f] = det_V[k*nbFilter+f];
			det_V[k*nbFilter+f] = 0.;
			ClusterLabel_tmp[k*nbFilter+f] = ClusterLabel[k*nbFilter+f];
			ClusterLabel[k*nbFilter+f] = 0.;
			for(unsigned int d1=0;d1<dim;d1++){
				mu_tmp[k*dim*nbFilter+d1*nbFilter+f] = mu[k*dim*nbFilter+d1*nbFilter+f];
				mu[k*dim*nbFilter+d1*nbFilter+f] = 0.;
				for(unsigned int d2=0;d2<dim;d2++){
					V_inv_tmp[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+f] = V_inv[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+f];
					V_inv[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+f] = 0.;
				}
			}
		}
	}
	for(unsigned int f=0;f<nbFilter;f++){
		nbClust[f] = nbClust_tmp[FilterIndexToModify[f]];
		FilterValue[f] = FilterValue_tmp[FilterIndexToModify[f]];
		for(unsigned int k=0;k<nbClust[f];k++){
			weights[k*nbFilter+f] = weights_tmp[k*nbFilter+FilterIndexToModify[f]];
			det_V[k*nbFilter+f] = det_V_tmp[k*nbFilter+FilterIndexToModify[f]];
			ClusterLabel[k*nbFilter+f] = ClusterLabel_tmp[k*nbFilter+FilterIndexToModify[f]];
			for(unsigned int d1=0;d1<dim;d1++){
				mu[k*dim*nbFilter+d1*nbFilter+f] = mu_tmp[k*dim*nbFilter+d1*nbFilter+FilterIndexToModify[f]];
				for(unsigned int d2=0;d2<dim;d2++) V_inv[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+f] = V_inv_tmp[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+FilterIndexToModify[f]];
			}
		}
	}

	delete[] nbClust_tmp;
	delete[] FilterValue_tmp;
	delete[] weights_tmp;
	delete[] det_V_tmp;
	delete[] ClusterLabel_tmp;
	delete[] mu_tmp;
	delete[] V_inv_tmp;

}

// Readers and printers

void GaussianMixtureModel::PrintModelParams(string filename){
	ofstream writefile(filename);
	for(unsigned int l=0;l<nbLabel;l++){
	writefile << "For label " << _MyDescriptors->getLabels(l) << endl;
		for(unsigned int f=0;f<nbFilter;f++){
			if( IsRead ) writefile << "FILTER_VALUE " << FilterValue[f] << endl;
			else writefile << "FILTER_VALUE " << _MyDescriptors->getFilterValue(f) << endl;
			unsigned int nb = 0;
			for(unsigned int k=0;k<nbClust[f];k++) if( ClusterLabel[k*nbFilter+f] == l ) nb++;
			writefile << "NUMBER_OF_CLUSTER " << nb << endl;
			for(unsigned int k=0;k<nbClust[f];k++){
				if( ClusterLabel[k*nbFilter+f] == l ){
					writefile << "WEIGHT " << weights[k*nbFilter+f] << endl;
					writefile << "DETERMINANT_OF_COV_MATRIX " << det_V[k*nbFilter+f] << endl;
					writefile << "ESPERANCE";
					for(unsigned int d=0;d<dim;d++) writefile << " " << mu[k*dim*nbFilter+d*nbFilter+f];
				 	writefile << endl << "COVARIANCE_MATRIX" << endl;
					for(unsigned int d1=0;d1<dim;d1++){
						for(unsigned int d2=0;d2<dim;d2++) writefile << V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+f] << " ";
						writefile << endl;
					}
				}
			}
		}
	}
	writefile.close();
}

void GaussianMixtureModel::PrintToDatabase(const string &name_of_database){
	if( !IsLabelled && !IsRead ){
		cerr << "The GMM is not labelled, we then cannot print it to the database, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	if( !IsDescriptor ){
		if( !IsRead ){
			cerr << "The GMM does not have descriptors, we then cannot print to database the model, aborting" << endl;
			exit(EXIT_FAILURE);
		}else nbLabel = Labels.size();
	}else nbLabel = _MyDescriptors->getNbLabels();

	string path2base = getDatabasePath(name_of_database);	

	if( !IsRead ){
		ofstream writefile_train(path2base+"labelling.out");
		for(unsigned int f=0;f<nbFilter;f++){
			writefile_train << "For descriptor filter \"" << _MyDescriptors->getFilterValue(f) << "\"" << endl;
			for(unsigned int l=0;l<nbLabel;l++){
				writefile_train << "Average cluster probabilities for label : \"" << _MyDescriptors->getLabels(l) << "\":\t\t";
				for(unsigned int k=0;k<nbClust[f];k++) writefile_train << AveLabelProb[k*nbFilter*nbLabel+l*nbFilter+f] << "\t";
				writefile_train << endl;
			}	       
			writefile_train << endl;
		}
		writefile_train.close();
	}

	string ext=".ath";
	for(unsigned int l=0;l<nbLabel;l++){
		string full_filename;
		if( IsRead ) full_filename=Labels[l]+ext;
		else full_filename=_MyDescriptors->getLabels(l)+ext;
		ofstream writefile(path2base+full_filename);
		if( IsRead ){
			writefile << "DESCRIPTOR_NAME " << DescriptorName << endl;
			writefile << "FILTER_TYPE " << FilteringType << endl;
			for(unsigned int p=0;p<DescriptorProperties.size();p++) writefile << DescriptorProperties[p] << endl;
		}else _MyDescriptors->printDescriptorsPropToDatabase(writefile);

		for(unsigned int f=0;f<nbFilter;f++){
			if( IsRead ) writefile << "FILTER_VALUE " << FilterValue[f] << endl;
			else writefile << "FILTER_VALUE " << _MyDescriptors->getFilterValue(f) << endl;
			unsigned int nb = 0;
			for(unsigned int k=0;k<nbClust[f];k++) if( ClusterLabel[k*nbFilter+f] == l ) nb++;
			writefile << "NUMBER_OF_CLUSTER " << nb << endl;
			for(unsigned int k=0;k<nbClust[f];k++){
				if( ClusterLabel[k*nbFilter+f] == l ){
					writefile << "WEIGHT " << weights[k*nbFilter+f] << endl;
					writefile << "DETERMINANT_OF_COV_MATRIX " << det_V[k*nbFilter+f] << endl;
					writefile << "ESPERANCE";
					for(unsigned int d=0;d<dim;d++) writefile << " " << mu[k*dim*nbFilter+d*nbFilter+f];
				 	writefile << endl << "INVERSE_OF_COVARIANCE_MATRIX" << endl;
					for(unsigned int d1=0;d1<dim;d1++){
						for(unsigned int d2=0;d2<dim;d2++) writefile << V_inv[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+f] << " ";
						writefile << endl;
					}
				}
			}
		}
		writefile.close();
	}
}

void GaussianMixtureModel::ReadModelParamFromDatabase(const string &name_of_database){
	string path2base = getMLDatabasePath()+name+"/"+name_of_database+"/";	
	string ext=".ath";
	string end, buffer_s, line;
	struct dirent *diread;
	const char *env = path2base.c_str();
	DIR *dir;
	size_t pos_filter_type, pos_filter_val, pos_des_name, pos_dim, pos_name;
	unsigned int buffer_i;
	nbFilter = 0;
	dim = 0;
	if( (dir = opendir(env) ) != nullptr ){
		cout << "Reading \"" << name_of_database << "\" GMM database" << endl;
		while( (diread = readdir(dir)) != nullptr ){
			buffer_s = diread->d_name;
			if( buffer_s.size() > 4 ){
				end = buffer_s.substr(buffer_s.size()-4,buffer_s.size());
				if( end == ext ) Labels.push_back(buffer_s.substr(0,buffer_s.size()-4));
			}
		}
		nbLabel = Labels.size();
		cout << nbLabel << " different labels : ";
		for(unsigned int l=0;l<nbLabel;l++) cout << Labels[l] << " ";
		cout << endl;
		closedir(dir);
		vector<unsigned int> NClust_temp;
		for(unsigned int l=0;l<nbLabel;l++){
			unsigned int line_fval(1000), line_ftype(1000);
			unsigned int count(0), current_filter_val;
			string full_file_path = path2base+Labels[l]+ext;
			ifstream file(full_file_path.c_str(), ifstream::in);
			if( file ){
				while(getline(file,line)){
					pos_name=line.find("DESCRIPTOR_NAME ");
					if(pos_name!=string::npos){
						istringstream text(line);
						text >> buffer_s;
						text >> DescriptorName;
					}
					pos_filter_type=line.find("FILTER_TYPE ");
					if(pos_filter_type!=string::npos){
						line_ftype = count;
						istringstream text(line);
						text >> buffer_s;
						text >> buffer_s;
						if( nbFilter == 0 ) FilteringType = buffer_s;
						else if( buffer_s != FilteringType ){
							cerr << "The filtering type is not the same in the different labels, aborting" << endl;
							exit(EXIT_FAILURE);
						}
					}
					pos_dim=line.find("NUMBER_OF_DIMENSION ");
					if(pos_dim!=string::npos){
						istringstream text(line);
						text >> buffer_s;
						text >> dim;
						dim2 = dim*dim;
					}
					pos_filter_val=line.find("FILTER_VALUE ");
					if(pos_filter_val!=string::npos){
						line_fval = count;
						istringstream text(line);
						text >> buffer_s;
						text >> buffer_s;
						bool already = false;
						for(unsigned int f=0;f<nbFilter;f++){
							if( buffer_s == FilterValue[f] ){
								current_filter_val = f;
								already = true;
								break;
							}
						}
						if( !already ){
							nbFilter++;
							NClust_temp.push_back(0);
							FilterValue.push_back(buffer_s);
							current_filter_val = nbFilter-1;
						}
					}
					if( count == line_fval+1 ){
						istringstream text(line);
						text >> buffer_s;
						text >> buffer_i;
						NClust_temp[current_filter_val] += buffer_i;
					}
					count++;
				}
			}else{
				cerr << "The file " << Labels[l] << ".ath cannot be openned" << endl;
				exit(EXIT_FAILURE);
			}
		}
		if( dim == 0 ){
			cerr << "The number of dimension cannot be read from the database, aborting" << endl;
			exit(EXIT_FAILURE);
		}
		nbClust = new unsigned int[nbFilter];
		for(unsigned int f=0;f<nbFilter;f++){
			nbClust[f] = NClust_temp[f];
			NClust_temp[f] = 0;
		}
		weights = new double[nbMaxClusters*nbFilter];
		mu = new double[dim*nbMaxClusters*nbFilter];
		V_inv = new long double[dim2*nbMaxClusters*nbFilter];
		det_V = new long double[nbMaxClusters*nbFilter];
		ClusterLabel = new unsigned int[nbFilter*nbMaxClusters];
		IsRead = true;
		vector<string> Prop_temp;
		for(unsigned int l=0;l<nbLabel;l++){
			Prop_temp.clear();
			unsigned int line_ftype(1000), line_fval(1000);
			unsigned int count(0), current_filter_val;
			string filter_val;
			bool stop_read_prop = false;
			string full_file_path = path2base+Labels[l]+ext;
			ifstream file(full_file_path.c_str(), ifstream::in);
			if( file ){
				while(getline(file,line)){
					pos_filter_type=line.find("FILTER_TYPE ");
					if(pos_filter_type!=string::npos) line_ftype = count;
					pos_filter_val=line.find("FILTER_VALUE ");
					if(pos_filter_val!=string::npos){
						istringstream text(line);
						text >> buffer_s;
						text >> filter_val;
						line_fval = count;
						stop_read_prop = true;
					}
					if( !stop_read_prop ){
						if( l == 0 ){
							DescriptorProperties.push_back(line);
							Prop_temp.push_back(line);
						}else Prop_temp.push_back(line);
					}
					if( count > line_fval ){
						for(unsigned int f=0;f<nbFilter;f++){
							if( FilterValue[f] == filter_val ){
								current_filter_val = f;
								break;
							}
						}
						istringstream text(line);
						text >> buffer_s;
						unsigned int currentK;
						text >> currentK;
						for(unsigned int k=0;k<currentK;k++){
							ClusterLabel[NClust_temp[current_filter_val]*nbFilter+current_filter_val] = l;
							getline(file,line);
							istringstream text2(line);
							text2 >> buffer_s;
							if( buffer_s == "WEIGHT" ) text2 >> weights[NClust_temp[current_filter_val]*nbFilter+current_filter_val];
							else{
								cerr << "Issue when reading weight of GMM" << endl;
								exit(EXIT_FAILURE);
							}
							getline(file,line);
							istringstream text3(line);
							text3 >> buffer_s;
							if( buffer_s == "DETERMINANT_OF_COV_MATRIX" ) text3 >> det_V[NClust_temp[current_filter_val]*nbFilter+current_filter_val];
							else{
								cerr << "Issue when reading determinant of cov mat of GMM" << endl;
								exit(EXIT_FAILURE);
							}
							getline(file,line);
							istringstream text4(line);
							text4 >> buffer_s;
							if( buffer_s == "ESPERANCE" ) for(unsigned int d=0;d<dim;d++) text4 >> mu[NClust_temp[current_filter_val]*nbFilter*dim+d*nbFilter+current_filter_val];
							else{
								cerr << "Issue when reading esperance of GMM" << endl;
								exit(EXIT_FAILURE);
							}
							getline(file,line);
							istringstream text5(line);
							text5 >> buffer_s;
							if( buffer_s == "INVERSE_OF_COVARIANCE_MATRIX" ){
								for(unsigned int d1=0;d1<dim;d1++){
									getline(file,line);
									istringstream text6(line);
									for(unsigned int d2=0;d2<dim;d2++) text6 >> V_inv[NClust_temp[current_filter_val]*nbFilter*dim2+d1*nbFilter*dim+nbFilter*d2+current_filter_val];
								}
							}else{
								cerr << "Issue when reading inverse of cov mat of GMM" << endl;
								exit(EXIT_FAILURE);
							}
							NClust_temp[current_filter_val]++;
						}
					}
					// read clusters
					count++;
				}
				if( DescriptorProperties.size() == Prop_temp.size() ){
					bool same = true;
					for(unsigned int s=0;s<DescriptorProperties.size();s++){
						if( DescriptorProperties[s] != Prop_temp[s] ){
							same = false;
							break;
						}
					}
					if( !same ){
						cerr << "The properties of the descriptors are not the same in the different labels of the database, aborting" << endl;
						exit(EXIT_FAILURE);
					}
				}else{
						cerr << "The properties of the descriptors are not the same in the different labels of the database, aborting" << endl;
						exit(EXIT_FAILURE);
				}
			}else{
				cerr << "The file " << Labels[l] << ".ath cannot be openned" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}else{
		cerr << "The database environment \"" << path2base << "\" cannot be openned, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	cout << path2base << " GMM model successfully read !" << endl;
}

void GaussianMixtureModel::readFixedParams(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_nbCMax, pos_tol, pos_maxIter, pos_elfac, pos_nb_bic_inc, pos_after_el, pos_nbInit, pos_InitMethod;
	string buffer_s, line;
	unsigned int ReadOk(0);
	if(file){
		while(file){
			getline(file,line);
			pos_nbCMax=line.find("GMM_NB_MAX_CLUSTER");
			if(pos_nbCMax!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nbMaxClusters;
				ReadOk++;
			}
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
			pos_elfac=line.find("GMM_ELBOW_FACTOR");
			if(pos_elfac!=string::npos){
				istringstream text(line);
				text >> buffer_s >> fac_elbow;
				ReadOk++;
			}
			pos_nb_bic_inc=line.find("GMM_NB_BIC_INCREASE_FOR_MIN");
			if(pos_nb_bic_inc!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nb_bic_increase;
			}
			pos_after_el=line.find("GMM_NO_BIC_MIN_NO_ELBOW_CHOICE");
			if(pos_after_el!=string::npos){
				istringstream text(line);
				text >> buffer_s >> after_elbow_choice;
			}
			pos_nbInit=line.find("GMM_NB_INIT");
			if(pos_nbInit!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nbInit;
			}
			pos_InitMethod=line.find("GMM_INIT_METHOD");
			if(pos_InitMethod!=string::npos){
				istringstream text(line);
				text >> buffer_s >> InitMethod;
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

void GaussianMixtureModel::ReadProperties(vector<string> Properties){
	size_t pos_nbCMax, pos_tol, pos_maxIter, pos_elfac, pos_nb_bic_inc, pos_after_el, pos_nbInit, pos_InitMethod;
	string buffer_s;
	for(unsigned int i=0;i<Properties.size();i++){
		pos_nbCMax=Properties[i].find("GMM_NB_MAX_CLUSTER");
		if(pos_nbCMax!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> nbMaxClusters;
		}
		pos_tol=Properties[i].find("GMM_TOL_LKH_EM");
		if(pos_tol!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> tol_Lkh_EM;
		}
		pos_maxIter=Properties[i].find("GMM_MAX_ITER_EM");
		if(pos_maxIter!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> MaxIter_EM;
		}
		pos_elfac=Properties[i].find("GMM_ELBOW_FACTOR");
		if(pos_elfac!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> fac_elbow;
		}
		pos_nb_bic_inc=Properties[i].find("GMM_NB_BIC_INCREASE_FOR_MIN");
		if(pos_nb_bic_inc!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> nb_bic_increase;
		}
		pos_after_el=Properties[i].find("GMM_NO_BIC_MIN_NO_ELBOW_CHOICE");
		if(pos_after_el!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> after_elbow_choice;
		}
		pos_nbInit=Properties[i].find("GMM_NB_INIT");
		if(pos_nbInit!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> nbInit;
		}
		pos_InitMethod=Properties[i].find("GMM_INIT_METHOD");
		if(pos_InitMethod!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> InitMethod;
		}

	}
}

void GaussianMixtureModel::SetKMeansProperties(vector<string> Properties){
	KMeansProperties = Properties;
	IsKMeansProperties = true;
}

vector<string> GaussianMixtureModel::getAvailableDatabases(){
	string path2base = getMLDatabasePath()+name+"/";
	vector<string> baseAlreadySaved;
	string buffer_s;
	struct dirent *diread;
	const char *env = path2base.c_str();
	DIR *dir;
	if( (dir = opendir(env) ) != nullptr ){
		while( (diread = readdir(dir)) != nullptr ){
			buffer_s = diread->d_name;
			if( buffer_s != "." && buffer_s != ".." ) baseAlreadySaved.push_back(buffer_s);
		}
	}
	return baseAlreadySaved;
}

GaussianMixtureModel::~GaussianMixtureModel(){
	if( IsDescriptor ){
		delete[] weights;
		delete[] weights_old;
		delete[] mu;
		delete[] mu_old;
		delete[] V;
		delete[] V_inv;
		delete[] V_old;
		delete[] V_inv_old;
		delete[] det_V;
		delete[] det_V_old;
		delete[] D_i;
		delete[] C_di;
		delete[] buffer_di;
		delete[] E_d;
		delete[] BIC;
		delete[] LogLikelihood;
		delete[] nbClust;
	}
	if( SavedVariablesInitialized ){
		delete[] saved_bic;
		delete[] saved_nbClust;
		delete[] saved_weights;
		delete[] saved_mu;
		delete[] saved_V;
		delete[] saved_V_inv;
		delete[] saved_det_V;
	}
	if( TrainSavedVariablesInitialized ){
		delete[] train_saved_bic;
		delete[] train_saved_weights;
		delete[] train_saved_mu;
		delete[] train_saved_V;
		delete[] train_saved_V_inv;
		delete[] train_saved_det_V;
	}

	if( IsKMeans ){
		delete MyKM;
	}
	if( IsLabelled ){
		delete[] AveLabelProb;
		delete[] ClusterLabel;
	}
	if( IsRead ){
		if( !IsDescriptor ){
			delete[] nbClust;
			delete[] weights;
			delete[] mu;
			delete[] V_inv;
			delete[] det_V;
		}
		delete[] ClusterLabel;
	}
	if( IsClassified ){
		delete[] Classificator;
	}
}
