//**********************************************************************************
//*   DBScan.cpp                                                                   *
//**********************************************************************************
//* This file contains the implementation of the DBScan class                      *
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


#include "DBScan.h"
#include "AtomHicConfig.h"
#include <cmath>

using namespace std; 

DBScan::DBScan(){
	this->name = "DBScan"; 
	readFixedParams();
}

void DBScan::setDescriptors(Descriptors *D){
	MachineLearningModel::setDescriptors(D);
	this->IsDescriptor = true;
	nbClust = new unsigned int[nbFilter];
}

void DBScan::TrainModel(string filter_name){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, training aborted" << endl;
		exit(EXIT_FAILURE);
	}

	if( !_MyDescriptors->getIsNeighbours(filter_name) || _MyDescriptors->get_current_rc(filter_name) != eps ){
		_MyDescriptors->searchNeighbours(eps,filter_name);
	}
	unsigned int current_filter;
	bool found = false;
	for(unsigned int f=0;f<nbFilter;f++){
		if( FilterValue[f] == filter_name ){
			current_filter = f;
			found = true;
			break;
		}
	}
	if( !found ){
		cerr << "The provided filter value (" << filter_name << ") does not correspond to a filter value of the ML model, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	unsigned int current_filter_neigh;
	unsigned int nbNMax = _MyDescriptors->getNbMaxNAndFilter(filter_name,current_filter_neigh);
	unsigned int nbDes = nbDat[current_filter];
	double mean_minPts = 0.;
	for(unsigned int i=0;i<nbDes;i++) mean_minPts += _MyDescriptors->getNeighbours(current_filter_neigh)[i*(nbNMax+1)];
	mean_minPts /= (double) nbDes;
	minPts = (unsigned int) mean_minPts;
	minPts--;
	cout << "Average number of neighbours : " << minPts << endl;

	cout << "Performing density based clustering.." << endl;
	unsigned int nbDatTot = 0;
        for(unsigned int f=0;f<nbFilter;f++) nbDatTot += nbDat[f];
        if( !IsClassified ){
                IsClassified = true;
                Classificator = new double[2*nbDatTot]; // first column, id of cluster ([n*2]) : -1 => undefined, 0 => noisy, i => id of cluster, second column status of point ([n*2+1]) : core (1), outlier (0) point -10 not treated (filter)
		for(unsigned int i=0;i<nbDatTot*2;i++) Classificator[i] = -10;
        }
	// DBScan algorythm
	// initialize status of point
	for(unsigned int i=0;i<nbDes;i++){
		unsigned int ind = _MyDescriptors->getFilterIndex(current_filter*nbDatMax+i);
		Classificator[ind*2] = -1;
		if( _MyDescriptors->getNeighbours(current_filter_neigh)[i*(nbNMax+1)] > minPts ) Classificator[ind*2+1] = 1;
		else Classificator[ind*2+1] = 0;
	}
	vector<unsigned int> stack;
	unsigned int clust_count = 1;
	unsigned int *temp_arr_neigh = new unsigned int[nbNMax];
	unsigned int current_nbneigh;
	for(unsigned int i=0;i<nbDes;i++){
		unsigned int ind = _MyDescriptors->getFilterIndex(current_filter*nbDatMax+i);
		if( Classificator[ind*2] != -1 || Classificator[ind*2+1] == 0 ) continue;
		unsigned int this_i = i;
		while( true ){
			if( Classificator[ind*2] == -1 ){
				Classificator[ind*2] = clust_count;
				if( Classificator[ind*2+1] == 1 ){
					current_nbneigh = _MyDescriptors->getNeighbours(current_filter_neigh)[i*(nbNMax+1)];
					for(unsigned int n=0;n<current_nbneigh;n++) temp_arr_neigh[n] = _MyDescriptors->getNeighbours(current_filter_neigh)[i*(nbNMax+1)+1+n];
					for(unsigned int i2=0;i2<current_nbneigh;i2++){
						unsigned int ind2 = _MyDescriptors->getFilterIndex(current_filter*nbDatMax+temp_arr_neigh[i2]);
						if( Classificator[ind2*2] == -1 ) stack.push_back(temp_arr_neigh[i2]);
					}
				}
			}
			if( stack.size() == 0 ) break;
			i = stack.back();
			ind = _MyDescriptors->getFilterIndex(current_filter*nbDatMax+i);
			stack.pop_back();
		}
		clust_count++;
		i = this_i;
	}
	delete[] temp_arr_neigh;	
	cout << "Done !" << endl;
        cout << "Number of cluster found : " << clust_count-1 << endl;
	nbClust[current_filter] = clust_count-1;
}

unsigned int DBScan::getNbClust(string filter_name){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, we cannot return the number of cluster" << endl;
		exit(EXIT_FAILURE);
	}

	unsigned int current_filter;
	bool found = false;
	for(unsigned int f=0;f<nbFilter;f++){
		if( FilterValue[f] == filter_name ){
			current_filter = f;
			found = true;
			break;
		}
	}
	if( !found ){
		cerr << "The provided filter value (" << filter_name << ") does not correspond to a filter value of the ML model, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	return nbClust[current_filter];
}

unsigned int DBScan::getFilterValue(string filter_name){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, we cannot return the filter value" << endl;
		exit(EXIT_FAILURE);
	}

	unsigned int current_filter;
	bool found = false;
	for(unsigned int f=0;f<nbFilter;f++){
		if( FilterValue[f] == filter_name ){
			current_filter = f;
			found = true;
			break;
		}
	}
	if( !found ){
		cerr << "The provided filter name (" << filter_name << ") does not correspond to a filter name of the ML model, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	return current_filter;
}
void DBScan::ComputeMuAndV(string filter_name){
	if( !IsDescriptor ){
		cerr << "The ML model does not have descriptors, we cannot compute means and variances of clusters" << endl;
		exit(EXIT_FAILURE);
	}

	unsigned int current_filter;
	bool found = false;
	for(unsigned int f=0;f<nbFilter;f++){
		if( FilterValue[f] == filter_name ){
			current_filter = f;
			found = true;
			break;
		}
	}
	
	if( !found ){
		cerr << "The provided filter value (" << filter_name << ") does not correspond to a filter value of the ML model, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	
	if( !IsMuAndV ){
		centroids = new double[nbClustMax*dim*nbFilter];
		V = new long double[nbClustMax*dim*dim*nbFilter];
		IsMuAndV = true;
	}

	// iniitalize variances and centroids to zero
	unsigned int nbDes = nbDat[current_filter];
	unsigned int *nbDat2Clust = new unsigned int[nbClust[current_filter]]; // WARNING!
	for(unsigned int k=0;k<nbClust[current_filter];k++){
		nbDat2Clust[k] = 0;
		for(unsigned int d1=0;d1<dim;d1++){
			centroids[k*dim*nbFilter+d1*nbFilter+current_filter] = 0.;
			for(unsigned int d2=d1;d2<dim;d2++) V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+current_filter] = 0.;
		}
	}

	// compute centroids
	for(unsigned int i=0;i<nbDes;i++){
		unsigned int ind = _MyDescriptors->getFilterIndex(current_filter*nbDatMax+i);
		unsigned int current_clust_index;
		if( Classificator[ind*2] > 0 ) current_clust_index = ((unsigned int) Classificator[ind*2]) - 1;
		else continue;
		nbDat2Clust[current_clust_index] += 1;
		for(unsigned int d=0;d<dim;d++) centroids[current_clust_index*dim*nbFilter+d*nbFilter+current_filter] += _MyDescriptors->getDescriptors()[ind*dim+d];
	}	
	for(unsigned int k=0;k<nbClust[current_filter];k++){
		for(unsigned int d=0;d<dim;d++) centroids[k*dim*nbFilter+d*nbFilter+current_filter] /= nbDat2Clust[k];
	}

	// compute variances
	for(unsigned int i=0;i<nbDes;i++){
		unsigned int ind = _MyDescriptors->getFilterIndex(current_filter*nbDatMax+i);
		unsigned int current_clust_index;
		if( Classificator[ind*2] > 0 ) current_clust_index = ((unsigned int) Classificator[ind*2]) - 1;
		else continue;
		for(unsigned int d1=0;d1<dim;d1++){
			for(unsigned int d2=0;d2<dim;d2++){
				V[current_clust_index*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+current_filter] += ( _MyDescriptors->getDescriptors()[ind*dim+d1] - centroids[current_clust_index*dim*nbFilter+d1*nbFilter+current_filter] ) * ( _MyDescriptors->getDescriptors()[ind*dim+d2] - centroids[current_clust_index*dim*nbFilter+d2*nbFilter+current_filter] );
			}
		}
	}
	for(unsigned int k=0;k<nbClust[current_filter];k++){
		for(unsigned int d1=0;d1<dim;d1++){
			for(unsigned int d2=d1;d2<dim;d2++){
				unsigned int ind = k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+current_filter;
				if( nbDat2Clust[k] != 0 ) V[ind] /= nbDat2Clust[k];
				else V[ind] = 0.;
				V[k*dim2*nbFilter+d2*dim*nbFilter+d1*nbFilter+current_filter] = V[ind];
			}
		}
	}

	delete[] nbDat2Clust;

}

void DBScan::readFixedParams(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_eps, pos_minPts, pos_nbClustMax;
	string buffer_s, line;
	unsigned int ReadOk(0);
	if(file){
		while(file){
			getline(file,line);
			pos_nbClustMax=line.find("DBSCAN_NB_MAX_CLUSTER");
			if(pos_nbClustMax!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nbClustMax;
				ReadOk++;
			}
			pos_eps=line.find("DBSCAN_EPS");
			if(pos_eps!=string::npos){
				istringstream text(line);
				text >> buffer_s >> eps;
				ReadOk++;
			}
			pos_minPts=line.find("DBSCAN_MINPTS");
			if(pos_minPts!=string::npos){
				istringstream text(line);
				text >> buffer_s >> minPts;
				ReadOk++;
			}
		}
	}else{
		cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
		exit(EXIT_FAILURE);
	}
	file.close();
	if( ReadOk != 3 ){
		cerr << "Error during reading of FixedParameters.dat for DBScan, aborting" << endl;
		exit(EXIT_FAILURE);
	}
}

void DBScan::ReadProperties(vector<string> Properties){
	size_t pos_eps, pos_minPts;
	string buffer_s;
	for(unsigned int i=0;i<Properties.size();i++){
		pos_eps=Properties[i].find("DBSCAN_EPS");
		if(pos_eps!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> eps;
		}
		pos_minPts=Properties[i].find("DBSCAN_MINPTS");
		if(pos_minPts!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> minPts;
		}
	}
}

DBScan::~DBScan(){
	if( this->IsDescriptor ){
		delete[] nbClust;
	}
	if( IsMuAndV ){
		delete[] centroids;
		delete[] V;
	}
}
