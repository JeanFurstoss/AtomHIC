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
#include <filesystem>
#include <dirent.h>

using namespace std; 

KMeans::KMeans(){
	this->name = "KMeans"; // TODO the user could affect name to the model or it could be read from database
	readFixedParams();
}

void KMeans::setDescriptors(Descriptors *D){
	MachineLearningModel::setDescriptors(D);
	
	if( IsFilterIndexModified ) ChangeFilterIndex();
	if( !IsRead ){
		for(unsigned int n=0;n<nbFilter;n++){
			nbClust.push_back(0);
			centroids.push_back(vector<vector<double>>());
			ClusterLabel.push_back(vector<unsigned int>());
			IsLabelled.push_back(false);
			MyKMeansTools.push_back(nullptr);
			AveLabelProb.push_back(nullptr);
			IsKMeansTools.push_back(false);
		}
	}

}

void KMeans::TrainModel(unsigned int &nbClust_min, unsigned int &nbClust_max, const std::string &filter_value, const bool &softly_labeled){
	nbLabel = _MyDescriptors->getNbLabels();
	bool noLabel = false;
	if( nbLabel < 2 ){
		cout << "The descriptors are not labelled, or the label has only one value, the fitted KMeans wont then be labelled" << endl;
		noLabel = true;
	}
	if( noLabel ){
		_MyDescriptors->constructSubarrays(true,false);
		_dataMat = _MyDescriptors->getSubarray(current_nbDat,filter_value);
		fitOptimalKMeans(nbClust_min,nbClust_max);
	}else{
		if( !softly_labeled ){
			unsigned int current_f = getCurrentFIndex(filter_value);
			_MyDescriptors->constructSubarrays(RatioTestTrain,true,true);
			for(unsigned int l=0;l<nbLabel;l++){
				string current_label = _MyDescriptors->getLabels(l);
				cout << "Training KMeans model for descriptors : " << filter_value << ", and label : " << current_label << endl << endl;
				_dataMat = _MyDescriptors->getSubarray(current_nbDat,filter_value,current_label);
				string forfilename = filter_value+"_"+current_label;
				fitOptimalKMeans(nbClust_min,nbClust_max,forfilename);
				appendOptimalKMeans(filter_value);
				for(unsigned int n=0;n<optimal_nbClust;n++) ClusterLabel[current_f].push_back(l);
			}
			setKMeansTools(filter_value);
			ComputeAveLabelProb(filter_value);
			PrintLabelling(filter_value);
			IsLabelled[current_f] = true;
		}else{
			cout << "Training KMeans model for descriptors : " << filter_value << endl << endl;
			_MyDescriptors->constructSubarrays(RatioTestTrain,true,false);
			_dataMat = _MyDescriptors->getSubarray(current_nbDat,filter_value);
			fitOptimalKMeans(nbClust_min,nbClust_max,filter_value);
			appendOptimalKMeans(filter_value);
			setKMeansTools(filter_value);
			Labelling(filter_value);
		}
	}

}

void KMeans::TrainModel(unsigned int &nbClust_min, unsigned int &nbClust_max, const bool &softly_labeled){
	for(unsigned int f=0;f<nbFilter;f++) TrainModel(nbClust_min,nbClust_max,FilterValue[f],softly_labeled);
}

void KMeans::TrainModel(vector<unsigned int> &nbClust_min, vector<unsigned int> &nbClust_max, vector<std::string> &filter_value){
	for(unsigned int f=0;f<filter_value.size();f++) TrainModel(nbClust_min[f],nbClust_max[f],filter_value[f]);
}

void KMeans::TrainModel(unsigned int &_nbClust, std::string filter_value){
	unsigned int current_f = getCurrentFIndex(filter_value);
	nbClust[current_f] = _nbClust;
	if( IsKMeansTools[current_f] ) delete MyKMeansTools[current_f];
	MyKMeansTools[current_f] = new KMeansTools(nbClust[current_f],_dataMat,current_nbDat,dim);
	IsKMeansTools[current_f] = true;
	MyKMeansTools[current_f]->ReadProperties(current_Properties);
	MyKMeansTools[current_f]->fit();
}

void KMeans::fitOptimalKMeans(unsigned int &nbClust_min, unsigned int &nbClust_max, string namefile){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, fitting aborted" << endl;
		exit(EXIT_FAILURE);
	}
	unsigned int nbRuns = nbClust_max - nbClust_min + 1;
	vector<KMeansTools*> MyKMeans;
	vector<double> SILs(nbRuns);
	for(unsigned int n=0;n<nbRuns;n++){
		unsigned int current_nbclust = nbClust_min+n;
		MyKMeans.push_back(new KMeansTools(current_nbclust,_dataMat,current_nbDat,dim));
		MyKMeans[n]->ReadProperties(current_Properties);
	}
	
	//pragma
	for(unsigned int n=0;n<nbRuns;n++){
		cout << "Number of cluster = " << n+nbClust_min << endl;
		MyKMeans[n]->fit();
		SILs[n] = MyKMeans[n]->ComputeSilhouette();
	}
	cout << endl;
	// search the index giving the best KMeans 
	unsigned int optimal_index = 0;
	
	// 1. See if there is a minimum in the silhouette curve (the silhouette has to increase nb_silhouette_increase times after a minimum to be considered as the optimal -> defined in FixedParameters, default 1)
	unsigned int nb_increase = 0;
	double current_sil_low = SILs[0];
	bool minFound = false;
	for(unsigned int n=1;n<nbRuns;n++){
		if( current_sil_low < SILs[n] ){
			nb_increase++;
			if( nb_increase >= nb_sil_increase ){
				minFound = true;
				break;
			}
		}else{
			current_sil_low = SILs[n];
			optimal_index = n;
		}
	}

	// 2. if no minimum, search if there is an Elbow (i.e. a change in the slope of the SIL=f(N) curve lower than fac_elbow -> defined in FixedParameters, default 0.1)
	bool ElbowFound = false;
	if( !minFound ){
		optimal_index = 1;
		double current_slope = fabs(SILs[1]-SILs[0]);
		for(unsigned int n=1;n<nbRuns-1;n++){
			double new_slope = fabs(SILs[n+1]-SILs[n]);
			if( fac_elbow*current_slope > new_slope ){
				optimal_index = n;
				ElbowFound = true;
				break;
			}else current_slope = new_slope;
		}
		// 3. if no Elbow use the choice of the user (i.e. max, min or a given number of cluster)
		if( !ElbowFound ){
			if( after_elbow_choice == "Max" ){
				optimal_index = nbRuns-1;
			}else if( after_elbow_choice == "Min" ){
				optimal_index = 0.;
			}else{
				unsigned int nbclust_to_choose;
				istringstream iss_N(after_elbow_choice);
				iss_N >> nbclust_to_choose;
				if( nbclust_to_choose < nbClust_min ){
					cout << "The provided number of cluster by KMEANS_NO_SIL_MIN_NO_ELBOW_CHOICE in /data/FixedParameters/FixedParameters.dat is lower than nbClust_min, we then selected N = nbClust_max" << endl;
					optimal_index = nbRuns-1;
				}else if( nbclust_to_choose > nbClust_max ){
					cout << "The provided number of cluster by KMEANS_NO_SIL_MIN_NO_ELBOW_CHOICE in /data/FixedParameters/FixedParameters.dat is higher than nbClust_max, we then selected N = nbClust_max" << endl;
					optimal_index = nbRuns-1;
				}else optimal_index = nbclust_to_choose-nbClust_min;
			}
			cout << "As no minimum in the SIL=f(N) curve and no Elbow have been detected, the optimal number of cluster for this KMEANS (" << optimal_index+nbClust_min << ") has been chosen following KMEANS_NO_SIL_MIN_NO_ELBOW_CHOICE parameter in FixedParameters" << endl;
		}else cout << "The optimal number of cluster for this KMeans (" << optimal_index+nbClust_min << ") has been chosen with the Elbow method" << endl;
	}else cout << "The optimal number of cluster for this KMeans (" << optimal_index+nbClust_min << ") has been chosen using a minimum in the SIL=f(N) curve" << endl;

	if( namefile != "" ) namefile = "_" + namefile;
	ofstream writefile("SIL_vs_N"+namefile+".dat");
	writefile << "Number_of_Cluster SIL" << endl;
	for(unsigned int i=0;i<nbRuns;i++) writefile << nbClust_min+i << " " << SILs[i] << endl;
	writefile.close();
	cout << "SIL_vs_N" << namefile << ".dat file successfuly writted" << endl;
	cout << endl;


	// clean variables of previous optimal KMeans
	for(unsigned int n=0;n<optimal_centroids.size();n++) optimal_centroids[n].clear();
	optimal_centroids.clear();
	
	// store new optimal KMeans
	optimal_nbClust = optimal_index+nbClust_min;
	for(unsigned int n=0;n<optimal_nbClust;n++){
		optimal_centroids.push_back(vector<double>());
		for(unsigned int d1=0;d1<dim;d1++) optimal_centroids[n].push_back(MyKMeans[optimal_index]->getCentroids()(n,d1));	
	}

	for(unsigned int n=0;n<nbRuns;n++){
		delete MyKMeans[0];
		MyKMeans.erase(MyKMeans.begin());
	}
}

void KMeans::appendOptimalKMeans(const string &filter_value){
	unsigned int current_f = getCurrentFIndex(filter_value);
	nbClust[current_f] += optimal_nbClust;
	for(unsigned int n=0;n<optimal_nbClust;n++){
		centroids[current_f].push_back(vector<double>());
		unsigned int current_size = centroids[current_f].size()-1;
		for(unsigned int d1=0;d1<dim;d1++) centroids[current_f][current_size].push_back(optimal_centroids[n][d1]);
	}
}

void KMeans::setKMeansTools(string filter_value){
	unsigned int current_f = getCurrentFIndex(filter_value);
	if( IsKMeansTools[current_f] ) delete MyKMeansTools[current_f];
	MyKMeansTools[current_f] = new KMeansTools(nbClust[current_f],_dataMat,current_nbDat,dim,centroids[current_f]);
	MyKMeansTools[current_f]->ReadProperties(current_Properties);
	IsKMeansTools[current_f] = true;
}

void KMeans::ComputeAveLabelProb(string filter_value){
	unsigned int current_f = getCurrentFIndex(filter_value);
	if( AveLabelProb[current_f] ) delete[] AveLabelProb[current_f];

	AveLabelProb[current_f] = new long double[nbClust[current_f]*nbLabel];

	double zero = 0.;
	_MyDescriptors->constructSubarrays(zero,true,false);
	_dataMat = _MyDescriptors->getSubarray(current_nbDat,filter_value);
	MyKMeansTools[current_f]->setDataMat(_dataMat,current_nbDat);
	CorresIndexDescriptors = _MyDescriptors->getCorresIndexSubarray(current_nbDat,filter_value);
	
	MyKMeansTools[current_f]->AffectData2Cluster();
	unsigned int *data2Cluster = MyKMeansTools[current_f]->getData2Cluster();
	
	unsigned int *NbDesLabel = new unsigned int[nbLabel];
	
	for(unsigned int l=0;l<nbLabel;l++){
		NbDesLabel[l] = 0;
		for(unsigned int k=0;k<nbClust[current_f];k++) AveLabelProb[current_f][k*nbLabel+l] = 0.;
	}

	for(unsigned int i=0;i<current_nbDat;i++){
		NbDesLabel[_MyDescriptors->getLabels_uint(CorresIndexDescriptors[i])] += 1;
		AveLabelProb[current_f][data2Cluster[i]*nbLabel+_MyDescriptors->getLabels_uint(CorresIndexDescriptors[i])] += 1;
	}
	
	for(unsigned int k=0;k<nbClust[current_f];k++)
		for(unsigned int l=0;l<nbLabel;l++) AveLabelProb[current_f][k*nbLabel+l] /= NbDesLabel[l];
	delete[] NbDesLabel;
}

void KMeans::PrintLabelling(string filter_value){
	unsigned int current_f = getCurrentFIndex(filter_value);
	bool printWarning = false;
	cout << "For descriptor filter \"" << _MyDescriptors->getFilterValue(current_f) << "\"" << endl;
	vector<vector<string>> element_arr;
	element_arr.push_back(vector<string>());
	element_arr[0].push_back("");
	for(unsigned int l1=0;l1<nbLabel;l1++) element_arr[0].push_back(_MyDescriptors->getLabels(l1));
	for(unsigned int l1=0;l1<nbLabel;l1++){
		element_arr.push_back(vector<string>());
		element_arr[l1+1].push_back(_MyDescriptors->getLabels(l1));
		for(unsigned int l2=0;l2<nbLabel;l2++){
			double sum = 0.;
			for(unsigned int k=0;k<nbClust[current_f];k++){
				if( ClusterLabel[current_f][k] == l2 ) sum += AveLabelProb[current_f][k*nbLabel+l1];
			}
			if( ( l1 != l2 ) && ( sum > tol2ndProb ) ) printWarning = true;
			element_arr[l1+1].push_back(to_string(sum));
		}
	}
	vector<vector<unsigned int>> fusion_arr(nbLabel+1, vector<unsigned int>(nbLabel+1,1));
	Dis.DisplayArray(element_arr,fusion_arr);
	
	if( printWarning )
		cout << "Warning, the labels seem to be not well separated in the descriptor space !" << endl;
}

void KMeans::Labelling(){
	for(unsigned int f=0;f<nbFilter;f++) Labelling(FilterValue[f]);
}

void KMeans::Labelling(string filter_value){
	unsigned int current_f = getCurrentFIndex(filter_value);
	unsigned int zero_lab = 0;
	unsigned int max = _MyDescriptors->getLabelsSize(current_f,zero_lab);
	unsigned int min = _MyDescriptors->getLabelsSize(current_f,zero_lab);
	for(unsigned int l=1;l<nbLabel;l++){
		 if( _MyDescriptors->getLabelsSize(current_f,l) > max ) max = _MyDescriptors->getLabelsSize(current_f,l);
		 else if( _MyDescriptors->getLabelsSize(current_f,l) < min ) min = _MyDescriptors->getLabelsSize(current_f,l);
	}
	if( max*tolLabelSize > min ) cout << "Warning the \"" << filter_value << "\" descriptors are not homogeneously distributed within the labels (high difference between number of descriptors in each labels), which could biased the results obtained with the fitted KMeans" << endl;
	
	if( IsLabelled[current_f] ) ClusterLabel[current_f].clear();
	ComputeAveLabelProb(filter_value);

	// Search which label has the highest probability in each cluster
	for(unsigned int k=0;k<nbClust[current_f];k++){
		long double max = AveLabelProb[current_f][k*nbLabel];
		unsigned int ind = 0;
		for(unsigned int l=1;l<nbLabel;l++){
			if( AveLabelProb[current_f][k*nbLabel+l] > max ){
				max = AveLabelProb[current_f][k*nbLabel+l];
				ind = l;
			}
		}
		ClusterLabel[current_f].push_back(ind);
	}

	IsLabelled[current_f] = true;
	PrintLabelling(filter_value);
}

void KMeans::Classify(){
	bool isFullyLabelled = true;
	for(unsigned int f=0;f<nbFilter;f++){
		Classify(FilterValue[f]);
		isFullyLabelled *= IsLabelled[f];
	}
	if( isFullyLabelled || IsRead ){
        	// write the StructureIndex.txt file
        	ofstream writefile("StructureIndex.txt");
        	writefile << "Structures index used for the KMeans structural analysis" << endl;
        	for(unsigned int l=0;l<nbLabel;l++) writefile << l << " " << Labels[l] << endl;
        	writefile << nbLabel << " Not identified" << endl;
        	writefile.close();
		cout << "The StructureIndex.txt file containing the names of the labels has been printed" << endl;
	}
}

void KMeans::Classify(string filter_value){
	unsigned int current_f = getCurrentFIndex(filter_value);
	if( !IsDescriptor ){
		cerr << "We dont have descriptor to classify, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Classifying the descriptors " << filter_value << endl;
	
	_MyDescriptors->constructSubarrays(true,false);
	CorresIndexDescriptors = _MyDescriptors->getCorresIndexSubarray(current_nbDat,filter_value);
	_dataMat = _MyDescriptors->getSubarray(current_nbDat,filter_value);
	MyKMeansTools[current_f]->setDataMat(_dataMat,current_nbDat);
	
	MyKMeansTools[current_f]->AffectData2Cluster();
	unsigned int *data2Cluster = MyKMeansTools[current_f]->getData2Cluster();
	
	unsigned int nbDatTot = 0;
	for(unsigned int f=0;f<nbFilter_descriptors;f++) nbDatTot += nbDat[f];
	if( !IsClassified ){
		IsClassified = true;
		unsigned int nbclass = 2*nbDatTot;
		Classificator = new double[nbclass];
		for(unsigned int i=0;i<nbclass;i++) Classificator[i] = -1.;
	}

	// Label classification	
	if( IsLabelled[current_f] || IsRead ){
		for(unsigned int j=0;j<current_nbDat;j++) Classificator[CorresIndexDescriptors[j]*2] = ClusterLabel[current_f][data2Cluster[j]]; // TODO search how to estimate error (i.e. second column of Classificator)
		cout << "Done" << endl;
	}else{ // Cluster classification
		for(unsigned int j=0;j<current_nbDat;j++) Classificator[CorresIndexDescriptors[j]*2] = data2Cluster[j]; // TODO search how to estimate error (i.e. second column of Classificator)
		cout << "Done" << endl;
	}
}

void KMeans::ChangeFilterIndex(){
	vector<unsigned int> nbClust_tmp(nbFilter);
	vector<string> FilterValue_tmp(nbFilter);
	vector<vector<unsigned int>> ClusterLabel_tmp(nbFilter);
	vector<vector<vector<double>>> centroids_tmp(nbFilter);

	// Copy data and reinitialize them
	for(unsigned int f=0;f<nbFilter;f++){
		nbClust_tmp[f] = nbClust[f];
		nbClust[f] = 0;
		FilterValue_tmp[f] = FilterValue[f];
		FilterValue[f] = 0.;
		for(unsigned int k=0;k<nbClust_tmp[f];k++){
			ClusterLabel_tmp[f].push_back(ClusterLabel[f][k]);
			ClusterLabel[f][k] = 0.;
			centroids_tmp[f].push_back(vector<double>());
			for(unsigned int d1=0;d1<dim;d1++){
				centroids_tmp[f][k].push_back(centroids[f][k][d1]);
				centroids[f][k][d1] = 0.;
			}
		}
	}
	for(unsigned int f=0;f<nbFilter;f++){
		nbClust[f] = nbClust_tmp[FilterIndexToModify[f]];
		FilterValue[f] = FilterValue_tmp[FilterIndexToModify[f]];
		for(unsigned int k=0;k<nbClust[f];k++){
			ClusterLabel[f][k] = ClusterLabel_tmp[FilterIndexToModify[f]][k];
			for(unsigned int d1=0;d1<dim;d1++) centroids[f][k][d1] = centroids_tmp[FilterIndexToModify[f]][k][d1];
		}
	}
}

// Readers and printers

void KMeans::PrintModelParams(string filename){
	ofstream writefile(filename);
	for(unsigned int l=0;l<nbLabel;l++){
	writefile << "For label " << _MyDescriptors->getLabels(l) << endl;
		for(unsigned int f=0;f<nbFilter;f++){
			if( IsRead ) writefile << "FILTER_VALUE " << FilterValue[f] << endl;
			else writefile << "FILTER_VALUE " << _MyDescriptors->getFilterValue(f) << endl;
			unsigned int nb = 0;
			for(unsigned int k=0;k<nbClust[f];k++) if( ClusterLabel[f][k] == l ) nb++;
			writefile << "NUMBER_OF_CLUSTER " << nb << endl;
			for(unsigned int k=0;k<nbClust[f];k++){
				if( ClusterLabel[f][k] == l ){
					writefile << "CENTROID";
					for(unsigned int d=0;d<dim;d++) writefile << " " << centroids[f][k][d];
					writefile << endl;
				}
			}
		}
	}
	writefile.close();
}

//void KMeans::PrintModelParams(string filename, vector<string> label_order){
//	ofstream writefile(filename);
//	for(unsigned int lo=0;lo<label_order.size();lo++){
//		for(unsigned int l=0;l<nbLabel;l++){
//			if( label_order[lo] == _MyDescriptors->getLabels(l) ){
//				writefile << "For label " << _MyDescriptors->getLabels(l) << endl;
//				for(unsigned int f=0;f<nbFilter;f++){
//					if( IsRead ) writefile << "FILTER_VALUE " << FilterValue[f] << endl;
//					else writefile << "FILTER_VALUE " << _MyDescriptors->getFilterValue(f) << endl;
//					unsigned int nb = 0;
//					for(unsigned int k=0;k<nbClust[f];k++) if( ClusterLabel[f][k] == l ) nb++;
//					writefile << "NUMBER_OF_CLUSTER " << nb << endl;
//					// print clusters in order with increasing weight
//					vector<double> clust_id_to_print;
//					for(unsigned int k=0;k<nbClust[f];k++){
//						if( ClusterLabel[f][k] == l ){
//							clust_id_to_print.push_back(k);
//							clust_id_to_print.push_back(weights[f][k]);
//						}
//					}
//					MT->sort(clust_id_to_print,1,2,clust_id_to_print);
//					unsigned int current_nb_clust = clust_id_to_print.size()/2;
//					for(unsigned int k=0;k<current_nb_clust;k++){
//						writefile << "CENTROID";
//						for(unsigned int d=0;d<dim;d++) writefile << " " << centroids[f][((unsigned int) clust_id_to_print[k*current_nb_clust])][d];
//						writefile << endl;
//					}
//				}
//				break;
//			}
//		}
//	}
//	writefile.close();
//}

void KMeans::PrintToDatabase(const string &name_of_database){
	bool isFullyLabelled = true;
	for(unsigned int f=0;f<nbFilter;f++) isFullyLabelled *= IsLabelled[f];
	if( !isFullyLabelled && !IsRead ){
		cerr << "The KMeans is not labelled, we then cannot print it to the database, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	if( !IsDescriptor ){
		if( !IsRead ){
			cerr << "The KMeans does not have descriptors, we then cannot print to database the model, aborting" << endl;
			exit(EXIT_FAILURE);
		}else nbLabel = Labels.size();
	}else nbLabel = _MyDescriptors->getNbLabels();

	string path2base = getDatabasePath(name_of_database);	

	if( !IsRead ){
		ofstream writefile_train(path2base+"labelling.out");
		writefile_train << "Confusion matrix after labelling of the KMeans:" << endl;
		for(unsigned int f=0;f<nbFilter;f++){
			writefile_train << "\tFor descriptor filter \"" << _MyDescriptors->getFilterValue(f) << "\"" << endl;
			vector<vector<string>> element_arr;
			element_arr.push_back(vector<string>());
			element_arr[0].push_back("");
			for(unsigned int l1=0;l1<nbLabel;l1++) element_arr[0].push_back(_MyDescriptors->getLabels(l1));
			for(unsigned int l1=0;l1<nbLabel;l1++){
				element_arr.push_back(vector<string>());
				element_arr[l1+1].push_back(_MyDescriptors->getLabels(l1));
				for(unsigned int l2=0;l2<nbLabel;l2++){
					double sum = 0.;
					for(unsigned int k=0;k<nbClust[f];k++){
						if( ClusterLabel[f][k] == l2 ) sum += AveLabelProb[f][k*nbLabel+l1];
					}
					element_arr[l1+1].push_back(to_string(sum));
				}
			}
			vector<vector<unsigned int>> fusion_arr(nbLabel+1, vector<unsigned int>(nbLabel+1,1));
			Dis.DisplayArray(element_arr,fusion_arr,writefile_train);
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
			for(unsigned int k=0;k<nbClust[f];k++) if( ClusterLabel[f][k] == l ) nb++;
			writefile << "NUMBER_OF_CLUSTER " << nb << endl;
			for(unsigned int k=0;k<nbClust[f];k++){
				if( ClusterLabel[f][k] == l ){
					writefile << "CENTROID";
					for(unsigned int d=0;d<dim;d++) writefile << " " << centroids[f][k][d];
					writefile << endl;
				}
			}
		}
		writefile.close();
	}
}

void KMeans::ReadModelParamFromDatabase(const string &name_of_database){
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
		cout << "Reading \"" << name_of_database << "\" KMeans database" << endl;
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
		for(unsigned int f=0;f<nbFilter;f++){
			nbClust.push_back(NClust_temp[f]);
			centroids.push_back(vector<vector<double>>(NClust_temp[f]));
			ClusterLabel.push_back(vector<unsigned int>(NClust_temp[f]));
			for(unsigned int k=0;k<NClust_temp[f];k++)
				for(unsigned int d1=0;d1<dim;d1++) centroids[f][k].push_back(0.);
			IsLabelled.push_back(true);
			MyKMeansTools.push_back(nullptr);
			IsKMeansTools.push_back(true);
			NClust_temp[f] = 0;
		}
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
							ClusterLabel[current_filter_val][NClust_temp[current_filter_val]] = l;
							getline(file,line);
							istringstream text3(line);
							text3 >> buffer_s;
							if( buffer_s == "CENTROID" ) for(unsigned int d=0;d<dim;d++) text3 >> centroids[current_filter_val][NClust_temp[current_filter_val]][d];
							else{
								cerr << "Issue when reading centroids of KMeans" << endl;
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
		// initialize KMeansTools
		for(unsigned int f=0;f<nbFilter;f++) MyKMeansTools[f] = new KMeansTools(nbClust[f],dim,centroids[f]);
	}else{
		cerr << "The database environment \"" << path2base << "\" cannot be openned, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	cout << path2base << " KMeans model successfully read !" << endl;
}

void KMeans::readFixedParams(){
	MachineLearningModel::readFixedParams();
	//string fp;
	//#ifdef FIXEDPARAMETERS
	//fp = FIXEDPARAMETERS;
	//#endif
	//string backslash="/";
	//string filename=fp+backslash+FixedParam_Filename;
	//ifstream file(filename, ios::in);
	ifstream file(FixedParam_Filename, ios::in);
	size_t pos_elfac, pos_nb_sil_inc, pos_after_el;
	string buffer_s, line;
	if(file){
		while(file){
			getline(file,line);
			pos_elfac=line.find("KMEANS_ELBOW_FACTOR");
			if(pos_elfac!=string::npos){
				istringstream text(line);
				text >> buffer_s >> fac_elbow;
				current_Properties.push_back(line);
			}
			pos_nb_sil_inc=line.find("KMEANS_NB_SIL_INCREASE_FOR_MIN");
			if(pos_nb_sil_inc!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nb_sil_increase;
				current_Properties.push_back(line);
			}
			pos_after_el=line.find("KMEANS_NO_BIC_MIN_NO_ELBOW_CHOICE");
			if(pos_after_el!=string::npos){
				istringstream text(line);
				text >> buffer_s >> after_elbow_choice;
				current_Properties.push_back(line);
			}
		}
		file.close();
	}
	//else{
	//	cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
	//	exit(EXIT_FAILURE);
	//}
}

void KMeans::ReadProperties(vector<string> Properties){
	MachineLearningModel::ReadProperties(Properties);
	size_t pos_nbCMax, pos_tol, pos_maxIter, pos_elfac, pos_nb_sil_inc, pos_after_el, pos_nbInit, pos_InitMethod;
	size_t pos_kmtol, pos_kmmaxIter, pos_kmnbInit;
	string buffer_s;
	for(unsigned int i=0;i<Properties.size();i++){

		pos_elfac=Properties[i].find("KMEANS_ELBOW_FACTOR");
		if(pos_elfac!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> fac_elbow;
			// search if this properties is already stored in current_Properties and remove it to store the new value
			size_t thispos;
			for(unsigned int p=0;p<current_Properties.size();p++){
				thispos = current_Properties[p].find("KMEANS_ELBOW_FACTOR");
				if( thispos!=string::npos ) current_Properties.erase(current_Properties.begin()+p);
			}
			current_Properties.push_back(Properties[i]);
		}
		pos_nb_sil_inc=Properties[i].find("KMEANS_NB_SIL_INCREASE_FOR_MIN");
		if(pos_nb_sil_inc!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> nb_sil_increase;
			size_t thispos;
			for(unsigned int p=0;p<current_Properties.size();p++){
				thispos = current_Properties[p].find("KMEANS_NB_SIL_INCREASE_FOR_MIN");
				if( thispos!=string::npos ) current_Properties.erase(current_Properties.begin()+p);
			}
			current_Properties.push_back(Properties[i]);
		}
		pos_after_el=Properties[i].find("KMEANS_NO_SIL_MIN_NO_ELBOW_CHOICE");
		if(pos_after_el!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> after_elbow_choice;
			size_t thispos;
			for(unsigned int p=0;p<current_Properties.size();p++){
				thispos = current_Properties[p].find("KMEANS_NO_SIL_MIN_NO_ELBOW_CHOICE");
				if( thispos!=string::npos ) current_Properties.erase(current_Properties.begin()+p);
			}
			current_Properties.push_back(Properties[i]);
		}
		// Read properties of KMeansTools to propagate them to MyKMeanss if this function is used
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
}

KMeans::~KMeans(){
}
