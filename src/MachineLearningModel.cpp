//**********************************************************************************
//*   MachineLearningModel.cpp                                                     *
//**********************************************************************************
//* This file contains the implementation of the MachineLearningModel class        *
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


#include "MachineLearningModel.h"
#include "AtomHicConfig.h"
#include <filesystem>
#include <dirent.h>

using namespace std;

MachineLearningModel::MachineLearningModel(){
	MT = new MathTools;
}

void MachineLearningModel::setDescriptors(Descriptors *D){ 
	_MyDescriptors = D;
	if( IsRead ){ // Case where the ML model is read from database, check that descriptors are consistent with base informations
		if( FilteringType != _MyDescriptors->getFilteringType() ){
			cerr << "The filtering type is different between the descriptors and the read ML database, aborting" << endl;
			exit(EXIT_FAILURE);
		}
		//if( nbFilter != _MyDescriptors->getNbFilter() ){
		//	cerr << "The number of filter is different between the descriptors and the read ML database, aborting" << endl;
		//	exit(EXIT_FAILURE);
		//}
		nbFilter_descriptors = _MyDescriptors->getNbFilter();
		if( dim != _MyDescriptors->getDim() ){
			cerr << "The number of dimension is different between the descriptors and the read ML database, aborting" << endl;
			exit(EXIT_FAILURE);
		}
		if( this->IsDescriptor ) delete[] FilterIndexToModify;
		FilterIndexToModify = new unsigned int[nbFilter];
		for(unsigned int f1=0;f1<nbFilter;f1++){
			bool found = false;
			for(unsigned int f2=0;f2<nbFilter_descriptors;f2++){
				if( FilterValue[f1] == _MyDescriptors->getFilterValue(f2) ){
					found = true;       
					if( f1 != f2 ) IsFilterIndexModified = true;
					FilterIndexToModify[f1] = f2;
					break;
				}
			}
			if( !found ){
				cerr << "The filter value are different between the decriptors and the read ML database, aborting" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}else{
		dim = _MyDescriptors->getDim();
		dim2 = dim*dim;
		nbFilter = _MyDescriptors->getNbFilter();
		nbFilter_descriptors = _MyDescriptors->getNbFilter();
		FilteringType = _MyDescriptors->getFilteringType();
		for(unsigned int f=0;f<nbFilter;f++) FilterValue.push_back(_MyDescriptors->getFilterValue(f));
	}

	if( this->IsDescriptor ){
		delete[] nbDat;	
		delete[] buffer_vec_1_dim;
		delete[] buffer_vec_2_dim;
	}
	nbDat = new unsigned int[nbFilter_descriptors];
	for(unsigned int f=0;f<nbFilter_descriptors;f++) nbDat[f] = _MyDescriptors->getNbDat(f);
	
	nbDatMax = _MyDescriptors->getNbDatMax();
	buffer_vec_1_dim = new double[dim];
	buffer_vec_2_dim = new double[dim];
	this->IsDescriptor = true;
}

string MachineLearningModel::getMLDatabasePath(){
	string database;	
	string backslash="/";
	#ifdef MACHINELEARNINGMODELS_DATABASE
	database = MACHINELEARNINGMODELS_DATABASE;
	#endif
	if( database.empty() ){
		cerr << "Warning machine learning database environment is empty" << endl;
		exit(EXIT_FAILURE);
	} else {
		return database.c_str()+backslash;
	}
}

string MachineLearningModel::getDatabasePath(const string &name_of_database){
	string ML_database_path = getMLDatabasePath();
	string path2base = ML_database_path+name;
	string buffer_s;
	struct dirent *diread;
	const char *env = path2base.c_str();
	DIR *dir;
	vector<string> baseAlreadySaved;
	if( (dir = opendir(env) ) != nullptr ){
		while( (diread = readdir(dir)) != nullptr ){
			buffer_s = diread->d_name;
			if( buffer_s != "." && buffer_s != ".." ) baseAlreadySaved.push_back(buffer_s);
		}
		bool alreadyStored = false;
		for(unsigned int b=0;b<baseAlreadySaved.size();b++){
			if( baseAlreadySaved[b] == name_of_database ){
				alreadyStored = true;
				break;
			}
		}
		if( !alreadyStored ){
			closedir(dir);
			path2base += "/"+name_of_database+"/";
			if( !filesystem::create_directory(path2base) ){
				cout << "Error when creating directory of the database, we will print the base in the current directory" << endl;
				path2base = "";
			}else cout << "We will store the base in " << path2base << endl;
		}else{ // search for directory with the same name and an _ and digits at the end to set the name with higher digit
			cout << "The base \"" << name_of_database << "\" already exist in " << path2base << endl;
			unsigned int index_of_base(0);
			for(unsigned int b=0;b<baseAlreadySaved.size();b++){
				size_t pos_ = baseAlreadySaved[b].find_last_of("_");
				if( pos_!=string::npos){
					string beg = baseAlreadySaved[b].substr(0,pos_);
					if( beg == name_of_database ){
						string end = baseAlreadySaved[b].substr(pos_+1,baseAlreadySaved[b].size());
						bool isNumber = true;
						for(size_t t=0;t<end.size();t++){
							if( !isdigit(end[t]) ){
								isNumber = false;
								break;
							}
						}
						if( isNumber ){
							unsigned int temp_index;
							istringstream iss(end);
							iss >> temp_index;
							if( temp_index > index_of_base ) index_of_base = temp_index;
						}
					}
				}
			}
			index_of_base++;
			closedir(dir);
			path2base += "/"+name_of_database+"_"+to_string(index_of_base)+"/";
			if( !filesystem::create_directory(path2base) ){
				cout << "Error when creating directory of the database, we will print the base in the current directory" << endl;
				path2base = "";
			}else cout << "We will store the base in " << path2base << endl;
		}
	}else{
		cout << "The base environment cannot be openned" << endl;
	}
	return path2base;
}

void MachineLearningModel::PrintClassifiedData(string filename){
	if( !IsClassified ){
		cerr << "The data are not classified, we thus cannot print them, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	ofstream writefile(filename);
	unsigned int nbDatTot = 0;
	for(unsigned int f=0;f<nbFilter_descriptors;f++){
		for(unsigned int j=0;j<nbDat[f];j++){
			for(unsigned int d=0;d<dim;d++) writefile << _MyDescriptors->getDescriptors()[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*dim+d] << " ";
			writefile << Classificator[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*2] << " " << Classificator[_MyDescriptors->getFilterIndex(f*nbDatMax+j)*2+1] << endl;
		}
	}
	writefile.close();
	cout << "Classified data successfully printed in \"" << filename << "\" under format : DescriptorValues LabelIndex MaximumLikelihoodClassifier" << endl;
}

unsigned int MachineLearningModel::getCurrentFIndex(string filter_value){
	unsigned int current_f;
	bool found = false;
	for(unsigned int f=0;f<nbFilter;f++){
		if( FilterValue[f] == filter_value ){
			current_f = f;
			found = true;
			break;
		}
	}
	if( !found ){
		cerr << "The provided filter value does not correspond to a filter of the ML model, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	return current_f;
}

void MachineLearningModel::readFixedParams(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	size_t pos_rattrain;
	ifstream file(filename, ios::in);
	string buffer_s, line;
	if(file){
		while(file){
			getline(file,line);
			pos_rattrain=line.find("ML_RATIO_TEST_TRAIN");
			if(pos_rattrain!=string::npos){
				istringstream text(line);
				text >> buffer_s >> RatioTestTrain;
				current_Properties.push_back(line);
			}
		}
	}
}

void MachineLearningModel::ReadProperties(vector<string> &Properties){
	size_t pos_rattrain;
	string buffer_s;
	for(unsigned int i=0;i<Properties.size();i++){
		current_Properties.push_back(Properties[i]);
		pos_rattrain=Properties[i].find("ML_RATIO_TEST_TRAIN");
		if(pos_rattrain!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> RatioTestTrain;
		}
	}
}

vector<string> MachineLearningModel::getAvailableDatabases(){
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

MachineLearningModel::~MachineLearningModel(){
	delete MT;
	if( this->IsDescriptor ){
		delete[] nbDat;
		delete[] buffer_vec_1_dim;
		delete[] buffer_vec_2_dim;
	}
	if( FilterIndexToModify ) delete[] FilterIndexToModify;
}
