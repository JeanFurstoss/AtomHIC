#include "GaussianMixtureModel.h"
#include "AtomHicConfig.h"
#include <cmath>
#include <filesystem>
#include <dirent.h>

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
	cout << "KMeans initialization done" << endl;
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

void GaussianMixtureModel::InitFromKMeans(unsigned int &_nbClust, unsigned int &filter_value){
	if( !IsKMeans ){
		MyKM = new KMeans();
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

long double GaussianMixtureModel::Prob_Cluster(unsigned int &index_cluster, unsigned int &DescriptorIndex, unsigned int &filter_value){
	double sp = 0.;
	for(unsigned int j=0;j<dim;j++) buffer_vec_2_dim[j] = (_MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(filter_value*nbDatMax+DescriptorIndex)*dim+j]-mu[index_cluster*dim*nbFilter+j*nbFilter+filter_value]); //not sure
	for(unsigned int i=0;i<dim;i++){
		buffer_vec_1_dim[i] = 0.;
		for(unsigned int j=0;j<dim;j++) buffer_vec_1_dim[i] += V_inv[index_cluster*dim2*nbFilter+i*dim*nbFilter+j*nbFilter+filter_value]*buffer_vec_2_dim[j];
		sp += buffer_vec_1_dim[i]*buffer_vec_2_dim[i];
	}
	return ( weights[index_cluster*nbFilter+filter_value] / ( pow(2.*M_PI, (double) dim/2.) * sqrt(det_V[index_cluster*nbFilter+filter_value]) ) ) * exp( -.5*sp );
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

// Label the GMM by affecting to each cluster the label which has the highest probability in its and return the second highest value to see if there is overlapping between labels 
void GaussianMixtureModel::Labelling(){
	unsigned int nbLabel = _MyDescriptors->getNbLabels();
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
			for(unsigned int k=0;k<nbClust[f];k++) AveLabelProb[k*nbFilter*nbLabel+_MyDescriptors->getLabels_uint()[f*nbDatMax+i]*nbFilter+f] += MaximumLikelihoodClassifier(k,i,f);
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

void GaussianMixtureModel::PrintToDatabase(const string &name_of_database){
	if( !IsDescriptor ){
		cerr << "The GMM does not have descriptors, we then cannot print to database the model, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	if( !IsLabelled ){
		cerr << "The GMM is not labelled, we then cannot print it to the database, aborting" << endl;
		exit(EXIT_FAILURE);
	}

	string path2base = getDatabasePath(name_of_database);	

	unsigned int nbLabel = _MyDescriptors->getNbLabels();

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

	string ext=".ath";
	for(unsigned int l=0;l<nbLabel;l++){
		string full_filename=_MyDescriptors->getLabels(l)+ext;
		ofstream writefile(path2base+full_filename);
		_MyDescriptors->printDescriptorsPropToDatabase(writefile);
		for(unsigned int f=0;f<nbFilter;f++){
			writefile << "FILTER_VALUE " << _MyDescriptors->getFilterValue(f) << endl;
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

double GaussianMixtureModel::MaximumLikelihoodClassifier(unsigned int &index_cluster, unsigned int &DescriptorIndex, unsigned int &filter_value){
	long double sumProbs = 0.;
	for(unsigned int k=0;k<nbClust[filter_value];k++) sumProbs += Prob_Cluster(k,DescriptorIndex,filter_value);
	return Prob_Cluster(index_cluster,DescriptorIndex,filter_value) / sumProbs;
}

void GaussianMixtureModel::ReadModelParamFromDatabase(const std::string &name_of_database){
	string path2base = getMLDatabasePath()+name+"/"+name_of_database+"/";	
	string ext=".ath";
	string end, buffer_s;
	struct dirent *diread;
	const char *env = path2base.c_str();
	DIR *dir;
	if( (dir = opendir(env) ) != nullptr ){
		cout << "Reading \"" << name_of_database << "\" GMM database" << endl;
		while( (diread = readdir(dir)) != nullptr ){
			buffer_s = diread->d_name;
			if( buffer_s.size() > 4 ){
				end = buffer_s.substr(buffer_s.size()-4,buffer_s.size());
				if( end == ext ) Labels.push_back(buffer_s.substr(0,buffer_s.size()-4));
			}
		}
		cout << Labels.size() << " different labels : ";
		for(unsigned int l=0;l<Labels.size();l++) cout << Labels[l] << " ";
		cout << endl;
		closedir(dir);
		//for(unsigned int l=0;l<Labels.size();l++){	
	}else{
		cerr << "The database environment \"" << path2base << "\" cannot be openned, aborting" << endl;
		exit(EXIT_FAILURE);
	}
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
	if( IsLabelled ){
		delete[] AveLabelProb;
		delete[] ClusterLabel;
	}
}
