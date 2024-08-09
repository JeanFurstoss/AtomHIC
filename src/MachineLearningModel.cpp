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
	this->IsDescriptor = true;
	dim = _MyDescriptors->getDim();
	dim2 = dim*dim;
	nbDatMax = _MyDescriptors->getNbDatMax();
	buffer_vec_1_dim = new double[dim];
	buffer_vec_2_dim = new double[dim];
	nbFilter = _MyDescriptors->getNbFilter();
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

MachineLearningModel::~MachineLearningModel(){
	delete MT;
	delete[] buffer_vec_1_dim;
	delete[] buffer_vec_2_dim;
}
