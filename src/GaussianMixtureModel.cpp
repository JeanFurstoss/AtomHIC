#include "GaussianMixtureModel.h"
#include "AtomHicConfig.h"
#include <cmath>

using namespace std; 

GaussianMixtureModel::GaussianMixtureModel(){
	this->name = "GaussianMixtureModel"; // TODO the user could affect name to the model or it could be read from database
}

void GaussianMixtureModel::setDescriptors(Descriptors *D){
	MachineLearningModel::setDescriptors(D);

	weights = new double[nbMaxClusters*nbFilter];
	mu = new double[dim*nbMaxClusters*nbFilter];
	V = new long double[dim2*nbMaxClusters*nbFilter];
	V_inv = new long double[dim2*nbMaxClusters*nbFilter];
	det_V = new long double[nbMaxClusters*nbFilter];
	
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
	nbClust = new unsigned int[nbFilter];
	nbDat = new unsigned int[nbFilter];
	for(unsigned int f=0;f<nbFilter;f++) nbDat[f] = _MyDescriptors->getNbDat(f);
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
	if( nbRun < 3 ){
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
			if( saved_bic[current_save_index*nbFilter+f] > saved_bic[(current_save_index-1)*nbFilter+f] ){
			       cout << "We find a minimum in BIC, we then stop the computation and consider that the optimal GMM is for a numer of cluster = " << k-1 << endl;
			       OptModelFound = true;
			       unsigned int opt = current_save_index-1;
			       SetOptimalModel(opt,f);
			       break;
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
				cout << "We then selected N = number_of_cluster_max" << endl;
				cout << "Consider increasing number_of_cluster_max for GMM this training" << endl;
				opt = nbRun-1;
			}else{
				unsigned int opt_N = round( ( saved_bic[(nbRun-1)*nbFilter+f] - saved_bic[f] + ( (double) saved_nbClust[f] * diff_beg ) - ( (double) saved_nbClust[(nbRun-1)*nbFilter+f] * diff_end ) ) / ( diff_beg - diff_end ) );
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
	
	nbClust[filter_value] = _nbClust;
	
	InitFromKMeans(nbClust[filter_value],filter_value);
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
		iter++;
	}while( eps > tol_Lkh_EM && iter < MaxIter_EM );

	ComputeBIC(filter_value);

	if( eps > tol_Lkh_EM ){
		cout << "The EM algorythm does not converge until tolerance for a number of cluster = " << nbClust << ", final residual = " << eps << endl;
		cout << "Maybe increase number of iteration or the tolerance for EM in /data/FixedParameters/FixedParameters.dat" << endl;
	}
}

void GaussianMixtureModel::SetOptimalModel(unsigned int &opt_index, unsigned int &filter_value){
	BIC[filter_value] = saved_bic[opt_index*nbFilter+filter_value];
	nbClust[filter_value] = saved_nbClust[opt_index*nbFilter+filter_value];
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		weights[k*nbFilter+filter_value] = saved_weights[opt_index*nbMaxClusters*nbFilter+k*nbFilter+filter_value];
		det_V[k*nbFilter+filter_value] = saved_det_V[opt_index*nbMaxClusters*nbFilter+k*nbFilter+filter_value];
		for(unsigned int d1=0;d1<dim;d1++){
			mu[k*dim*nbFilter+d1*nbFilter+filter_value] = saved_mu[opt_index*nbMaxClusters*dim*nbFilter+k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = saved_V[opt_index*nbMaxClusters*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
				V_inv[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = saved_V_inv[opt_index*nbMaxClusters*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
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
		saved_weights[current*nbMaxClusters*nbFilter+k*nbFilter+filter_value] = weights[k*nbFilter+filter_value];
		saved_det_V[current*nbMaxClusters*nbFilter+k*nbFilter+filter_value] = det_V[k*nbFilter+filter_value];
		for(unsigned int d1=0;d1<dim;d1++){
			saved_mu[current*nbMaxClusters*dim*nbFilter+k*dim*nbFilter+d1*nbFilter+filter_value] = mu[k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				saved_V[current*nbMaxClusters*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
				saved_V_inv[current*nbMaxClusters*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = V_inv[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
			}
		}
	}
}

void GaussianMixtureModel::InitFromKMeans(unsigned int &_nbClust, unsigned int &filter_value){
	if( !IsKMeans ){
		MyKM = new KMeans();
		MyKM->setDescriptors(_MyDescriptors);
	}

	MyKM->TrainModel(nbClust[filter_value],filter_value);
	for(unsigned int k=0;k<nbClust[filter_value];k++){
		weights[k*nbFilter+filter_value] = (double) MyKM->getNbDat2Cluster()[k*nbFilter+filter_value] / (double) nbDat[filter_value];
		det_V[k*nbFilter+filter_value] = MyKM->getDet_V()[k*nbFilter+filter_value];
		for(unsigned int d1=0;d1<dim;d1++){
			mu[k*dim*nbFilter+d1*nbFilter+filter_value] = MyKM->getCentroids()[k*dim*nbFilter+d1*nbFilter+filter_value];
			for(unsigned int d2=0;d2<dim;d2++){
				V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = MyKM->getV()[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
				V_inv[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value] = MyKM->getV_inv()[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+filter_value];
			}
		}
	}
}

void GaussianMixtureModel::UpdateParams(unsigned int &filter_value){
	for(unsigned int d=0;d<nbClust[filter_value];d++){
		weights_old[d*nbFilter+filter_value] = weights[d*nbFilter+filter_value];
		det_V_old[d*nbFilter+filter_value] = det_V[d*nbFilter+filter_value];
		for(unsigned int x1=0;x1<dim;x1++){
			mu_old[d*dim*nbFilter+x1*nbFilter+filter_value] = mu[d*dim*nbFilter+x1*nbFilter+filter_value];
			for(unsigned int x2=0;x2<dim;x2++){
				V_old[d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value] = V[d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value];
				V_inv_old[d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value] = V_inv[d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value];
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
	for(unsigned int d=0;d<nbClust[filter_value];d++){
		// new weights
		weights[d*nbFilter+filter_value] = E_d[d]/nbDat[filter_value];
		// new means
		for(unsigned int x=0;x<dim;x++){
			mu[d*dim*nbFilter+x*nbFilter+filter_value] = 0.;
			for(unsigned int i=0;i<nbDat[filter_value];i++) mu[d*dim*nbFilter+x*nbFilter+filter_value] += C_di[(d*nbDat[filter_value])+i]*_MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+i*dim+x];
			mu[d*dim*nbFilter+x*nbFilter+filter_value] /= E_d[d];
		}
		// new variances
		for(unsigned int x1=0;x1<dim;x1++){
			for(unsigned int x2=x1;x2<dim;x2++){
				V[d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value] = 0.;
				for(unsigned int i=0;i<nbDat[filter_value];i++) V[d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value] += C_di[(d*nbDat[filter_value])+i]*( (_MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+i*dim+x1] - mu[d*dim*nbFilter+x1*nbFilter+filter_value]) * (_MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+i*dim+x2] - mu[d*dim*nbFilter+x2*nbFilter+filter_value]) );
				V[d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value] /= E_d[d];
			}
			// symmetric part
			for(unsigned int x2=x1+1;x2<dim;x2++){
				V[d*dim2*nbFilter+x2*dim*nbFilter+x1*nbFilter+filter_value] = V[d*dim2*nbFilter+x1*dim*nbFilter+x2*nbFilter+filter_value];
			}
		}
		// invert variances
		MT->invMat_LU(V,V_inv,dim,d,nbFilter,filter_value,det_V[d*nbFilter+filter_value]);
	}
}

long double GaussianMixtureModel::Prob_Cluster(unsigned int index_cluster, unsigned int DescriptorIndex, unsigned int &filter_value){
	double sp = 0.;
	for(unsigned int j=0;j<dim;j++) buffer_vec_2_dim[j] = (_MyDescriptors->getDescriptors()[filter_value*nbDatMax*dim+DescriptorIndex*dim+j]-mu[index_cluster*dim*nbFilter+j*nbFilter+filter_value]);
	for(unsigned int i=0;i<dim;i++){
		buffer_vec_1_dim[i] = 0.;
		for(unsigned int j=0;j<dim;j++) buffer_vec_1_dim[i] += V_inv[index_cluster*dim2*nbFilter+i*dim*nbFilter+j*nbFilter+filter_value]*buffer_vec_2_dim[j];
		sp += buffer_vec_1_dim[i]*buffer_vec_2_dim[i];
	}
	return ( weights[index_cluster*nbFilter+filter_value] / ( pow(2.*M_PI, (double) dim/2.) * sqrt(det_V[index_cluster*nbFilter+filter_value]) ) ) * exp( -.5*sp );
}

void GaussianMixtureModel::ComputeLogLikelihood(unsigned int &filter_value){
	LogLikelihood[filter_value] = 0.;
	for(unsigned int d=0;d<nbClust[filter_value];d++) LogLikelihood[filter_value] += Prob_Cluster(d,0,filter_value);
	LogLikelihood[filter_value] = log(LogLikelihood[filter_value]);
	for(unsigned int i=1;i<nbDat[filter_value];i++){
		long double sum = 0.;
		for(unsigned int d=0;d<nbClust[filter_value];d++) sum += Prob_Cluster(d,i,filter_value);
		LogLikelihood[filter_value] += log(sum);
	}
}

void GaussianMixtureModel::ComputeBIC(unsigned int &filter_value){
	unsigned int NbIndepParams = ( ( nbClust[filter_value] * ( ( dim*( dim + 1)/2) + dim + 1. ) ) - 1. );
	BIC[filter_value] = -( 2.*LogLikelihood[filter_value] ) + ( log((double) nbDat[filter_value]) * (double) NbIndepParams );
}

void GaussianMixtureModel::PrintModelParams(string filename, unsigned int &filter_value){
	ofstream writefile(filename);
	writefile << "The Gaussian Mixture Model is composed by " << nbClust << " clusters" << endl;
        for(unsigned int k=0;k<nbClust[filter_value];k++){
                writefile << "Cluster " << k+1 << endl;
                writefile << "weight = " << weights[k*nbFilter+filter_value] << endl;
                writefile << "mean = ";
                for(unsigned int i=0;i<dim;i++) writefile << mu[k*dim*nbFilter+i*nbFilter+filter_value] << " ";
                writefile << endl;
                writefile << "variance = " << endl;
                for(unsigned int i=0;i<dim;i++){
                        for(unsigned int j=0;j<dim;j++) writefile << V[k*dim2*nbFilter+i*dim*nbFilter+j*nbFilter+filter_value] << " ";
                        writefile << endl;
                }
        }
	writefile << endl;
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
	if( IsKMeans ){
		delete MyKM;
	}
}
