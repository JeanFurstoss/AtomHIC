#include "GaussianMixtureModel.h"
#include "AtomHicConfig.h"
#include <cmath>

using namespace std; 

GaussianMixtureModel::GaussianMixtureModel(){
	this->name = "GaussianMixtureModel"; // TODO the user could affect name to the model or it could be read from database
}

void GaussianMixtureModel::setDescriptors(Descriptors *D){
	MachineLearningModel::setDescriptors(D);
	weights = new double[nbMaxClusters];
	mu = new double[dim*nbMaxClusters];
	V = new long double[dim2*nbMaxClusters];
	V_inv = new long double[dim2*nbMaxClusters];
	det_V = new long double[nbMaxClusters];
	
	D_i = new double[nbDat];
	C_di = new double[nbMaxClusters*nbDat];
	buffer_di = new double[nbMaxClusters*nbDat];
	E_d = new double[nbMaxClusters];
	
	weights_old = new double[nbMaxClusters];
	mu_old = new double[dim*nbMaxClusters];
	V_old = new long double[dim2*nbMaxClusters];
	V_inv_old = new long double[dim2*nbMaxClusters];
	det_V_old = new long double[nbMaxClusters];
}

void GaussianMixtureModel::TrainModel(unsigned int &_nbClust_min, unsigned int &_nbClust_max){
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
		saved_bic = new double[nbRun];
		saved_nbClust = new unsigned int[nbRun];
		saved_weights = new double[nbMaxClusters*nbRun];
		saved_mu = new double[dim*nbMaxClusters*nbRun];
		saved_V = new long double[dim2*nbMaxClusters*nbRun];
		saved_V_inv = new long double[dim2*nbMaxClusters*nbRun];
		saved_det_V = new long double[nbMaxClusters*nbRun];
		SavedVariablesInitialized = true;
	}

	cout << "Training for = " << endl;
	cout << _nbClust_min << " clusters" << endl;

	TrainModel(_nbClust_min);

	cout << "Training done !" << endl;
	unsigned int current_save_index = 0;
	SaveVariables(current_save_index);
	
	bool OptModelFound = false;
	for(unsigned int k=_nbClust_min+1;k<_nbClust_max+1;k++){
		cout << k << " clusters" << endl;
		current_save_index++;
		TrainModel(k);
		cout << "Training done !" << endl;
		SaveVariables(current_save_index);
		if( saved_bic[current_save_index] > saved_bic[current_save_index-1] ){
		       cout << "We find a minimum in BIC, we then stop the computation and consider that the optimal GMM is for a numer of cluster = " << k-1 << endl;
		       OptModelFound = true;
		       unsigned int opt = current_save_index-1;
		       SetOptimalModel(opt);
		       break;
		}
	}
	if( !OptModelFound ){
		cout << "We didn't find a local minimum of the BIC reflecting an optimal GMM" << endl;
		// Elbow method, we search the intersection between the two lines formed by the two first and the two last moints (maybe refine later)
		double diff_beg = saved_bic[1]-saved_bic[0];
		double diff_end = saved_bic[nbRun-1]-saved_bic[nbRun-2];
		unsigned int opt;
		if( diff_end > diff_beg*fac_elbow ){
			cout << "The BIC vs number of cluster curve does not really present an elbow" << endl;
			cout << "We then selected N = number_of_cluster_max" << endl;
			cout << "Consider increasing number_of_cluster_max for GMM this training" << endl;
			opt = nbRun-1;
		}else{
			unsigned int opt_N = round( ( saved_bic[nbRun-1] - saved_bic[0] + ( (double) saved_nbClust[0] * diff_beg ) - ( (double) saved_nbClust[nbRun-1] * diff_end ) ) / ( diff_beg - diff_end ) );
			for(unsigned int i=0;i<nbRun;i++){
				if( saved_nbClust[i] == opt_N ){
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
		SetOptimalModel(opt);
		ofstream writefile("BIC_vs_N.dat");
		writefile << "Number_of_Cluster BIC" << endl;
		for(unsigned int i=0;i<nbRun;i++) writefile << saved_nbClust[i] << " " << saved_bic[i] << endl;
		writefile.close();
		cout << "Yo may give a look on the BIC_vs_N.dat file to see if the optimal number of cluster suites to you" << endl;
	}
}

void GaussianMixtureModel::TrainModel(unsigned int &_nbClust){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, training aborted" << endl;
		exit(EXIT_FAILURE);
	}
	if( _nbClust > nbMaxClusters ){
		cerr << "The number of cluster for GMM training is higher than the maximum number of cluster (see /data/FixedParameters/FixedParameters.dat to modify), aborting" << endl;
		exit(EXIT_FAILURE);
	}
	
	nbClust = _nbClust;
	
	InitFromKMeans(nbClust);
	UpdateParams();
	ComputeLogLikelihood();
	
	long double Lkh_old = LogLikelihood;
	double eps;
	unsigned int iter = 0;
	do{
		EM();
		UpdateParams();
		ComputeLogLikelihood();
		eps = LogLikelihood - Lkh_old;
		Lkh_old = LogLikelihood;
		iter++;
	}while( eps > tol_Lkh_EM && iter < MaxIter_EM );

	ComputeBIC();

	if( eps > tol_Lkh_EM ){
		cout << "The EM algorythm does not converge until tolerance for a number of cluster = " << nbClust << ", final residual = " << eps << endl;
		cout << "Maybe increase number of iteration or the tolerance for EM in /data/FixedParameters/FixedParameters.dat" << endl;
	}
}

void GaussianMixtureModel::SetOptimalModel(unsigned int &opt_index){
	BIC = saved_bic[opt_index];
	nbClust = saved_nbClust[opt_index];
	for(unsigned int k=0;k<nbClust;k++){
		weights[k] = saved_weights[opt_index*nbMaxClusters+k];
		det_V[k] = saved_det_V[opt_index*nbMaxClusters+k];
		for(unsigned int d1=0;d1<dim;d1++){
			mu[k*dim+d1] = saved_mu[opt_index*nbMaxClusters*dim+k*dim+d1];
			for(unsigned int d2=0;d2<dim;d2++){
				V[k*dim2+d1*dim+d2] = saved_V[opt_index*nbMaxClusters*dim2+k*dim2+d1*dim+d2];
				V_inv[k*dim2+d1*dim+d2] = saved_V_inv[opt_index*nbMaxClusters*dim2+k*dim2+d1*dim+d2];
			}
		}
	}
}

void GaussianMixtureModel::SaveVariables(unsigned int &current){
	if( !SavedVariablesInitialized ){
		cerr << "Saved variables for GMM optimization have not been initialized, something went wrong, aborting" << endl;
		exit(EXIT_FAILURE);
	}

	saved_bic[current] = BIC;
	saved_nbClust[current] = nbClust;
	for(unsigned int k=0;k<nbClust;k++){
		saved_weights[current*nbMaxClusters+k] = weights[k];
		saved_det_V[current*nbMaxClusters+k] = det_V[k];
		for(unsigned int d1=0;d1<dim;d1++){
			saved_mu[current*nbMaxClusters*dim+k*dim+d1] = mu[k*dim+d1];
			for(unsigned int d2=0;d2<dim;d2++){
				saved_V[current*nbMaxClusters*dim2+k*dim2+d1*dim+d2] = V[k*dim2+d1*dim+d2];
				saved_V_inv[current*nbMaxClusters*dim2+k*dim2+d1*dim+d2] = V_inv[k*dim2+d1*dim+d2];
			}
		}
	}
}

void GaussianMixtureModel::InitFromKMeans(unsigned int &_nbClust){
	if( !IsKMeans ){
		MyKM = new KMeans();
		MyKM->setDescriptors(_MyDescriptors);
	}

	MyKM->TrainModel(nbClust);
	for(unsigned int k=0;k<nbClust;k++){
		weights[k] = (double) MyKM->getNbDat2Cluster()[k] / (double) nbDat;
		det_V[k] = MyKM->getDet_V()[k];
		for(unsigned int d1=0;d1<dim;d1++){
			mu[k*dim+d1] = MyKM->getCentroids()[k*dim+d1];
			for(unsigned int d2=0;d2<dim;d2++){
				V[k*dim2+d1*dim+d2] = MyKM->getV()[k*dim2+d1*dim+d2];
				V_inv[k*dim2+d1*dim+d2] = MyKM->getV_inv()[k*dim2+d1*dim+d2];
			}
		}
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

void GaussianMixtureModel::EM(){
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

long double GaussianMixtureModel::Prob_Cluster(unsigned int index_cluster, unsigned int DescriptorIndex){
	double sp = 0.;
	for(unsigned int j=0;j<dim;j++) buffer_vec_2_dim[j] = (_MyDescriptors->getDescriptors()[DescriptorIndex*dim+j]-mu[index_cluster*dim+j]);
	for(unsigned int i=0;i<dim;i++){
		buffer_vec_1_dim[i] = 0.;
		for(unsigned int j=0;j<dim;j++) buffer_vec_1_dim[i] += V_inv[index_cluster*dim2+i*dim+j]*buffer_vec_2_dim[j];
		sp += buffer_vec_1_dim[i]*buffer_vec_2_dim[i];
	}
	return ( weights[index_cluster] / ( pow(2.*M_PI, (double) dim/2.) * sqrt(det_V[index_cluster]) ) ) * exp( -.5*sp );
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

void GaussianMixtureModel::PrintModelParams(string filename){
	ofstream writefile(filename);
	writefile << "The Gaussian Mixture Model is composed by " << nbClust << " clusters" << endl;
        for(unsigned int k=0;k<nbClust;k++){
                writefile << "Cluster " << k+1 << endl;
                writefile << "weight = " << weights[k] << endl;
                writefile << "mean = ";
                for(unsigned int i=0;i<dim;i++) writefile << mu[k*dim+i] << " ";
                writefile << endl;
                writefile << "variance = " << endl;
                for(unsigned int i=0;i<dim;i++){
                        for(unsigned int j=0;j<dim;j++) writefile << V[k*dim2+i*dim+j] << " ";
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
