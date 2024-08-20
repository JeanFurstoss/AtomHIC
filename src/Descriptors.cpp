#include "Descriptors.h"
#include "AtomHicConfig.h"
#include <iostream>
#include <dirent.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <filesystem>

using namespace std;


Descriptors::Descriptors(AtomicSystem *_MySystem):_MySystem(_MySystem){
	readFixedParams();
	MT = new MathTools();
	nbDatTot = _MySystem->getNbAtom();
	ConstructFilterIndexArray(_MySystem);
}

Descriptors::Descriptors(AtomicSystem *_MySystem, string& DescriptorName, string& _FilteringType):_MySystem(_MySystem){
	FilteringType = _FilteringType;
	MT = new MathTools();
	nbDatTot = _MySystem->getNbAtom();
	ConstructFilterIndexArray(_MySystem);
	unsigned int aux_size, aux_id;
	aux_id = _MySystem->getAuxIdAndSize(DescriptorName,aux_size);
	if( aux_size == 0 ){
		cerr << "The provided descriptor name cannot be found in the auxiliary properties of the atomic system, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	dim = aux_size;
	_Descriptors = _MySystem->getAux(aux_id);
	AreDescriptorsMine = false;
}

Descriptors::Descriptors(AtomicSystem *_MySystem, vector<string> _Properties):_MySystem(_MySystem){
	readFixedParams(); // read fixed parameters before in order to have all values initialized even if some are missing in _Properties
	readProperties(_Properties);
	MT = new MathTools();
	nbDatTot = _MySystem->getNbAtom();
	ConstructFilterIndexArray(_MySystem);
}

void Descriptors::ConstructFilterIndexArray(AtomicSystem *_MySystem){
	FilterValue.clear();
	if( FilteringType == "none" ){
		nbFilter = 1;
		FilterValue.push_back("none");
		nbDat.push_back(_MySystem->getNbAtom());
		FilterIndex = new unsigned int[nbDat[0]];
		for(unsigned int i=0;i<nbDat[0];i++) FilterIndex[i] = i;
	}else if( FilteringType == "element" || FilteringType == "type" ){
		nbFilter = _MySystem->getNbAtomType();
		for(unsigned int f=0;f<nbFilter;f++) FilterValue.push_back("");
		if( nbFilter > nbMaxFilter ) cout << "Warning there is " << nbFilter << " different filters, which is a lot" << endl;
		if( _MySystem->getAtomType(0) == "" ){
			vector<unsigned int> type_t2e;
			vector<string> element_t2e;
			unsigned int buffer_i;
			string buffer_s, line;
			cout << "The element names are not provided in the atomic system" << endl;
			ifstream file_e2t("Type2Element.ath", ifstream::in);
			if( file_e2t ){
				cout << "Reading the correspondance between type and element in Type2Element.ath" << endl;
				while(getline(file_e2t,line)){
					istringstream text(line);
					text >> buffer_i >> buffer_s;
					type_t2e.push_back(buffer_i);
					element_t2e.push_back(buffer_s);
				}
				for(unsigned int f=0;f<nbFilter;f++){
					nbDat.push_back(0);
					FilterValue[type_t2e[f]-1] = element_t2e[f];
				}
				FilteringType = "element";
			}else{
				cout << "You can provide a Type2Element.ath file giving the correspondance between type and element in the working directory (an example is present in /data/ExampleFiles/)" << endl;
				cout << "As this file has not been found, we will simply filter the descriptors by type which could lead to dramatic confusion depending on what you are doing" << endl;
				FilteringType = "type";
			}
		}else{
			for(unsigned int f=0;f<nbFilter;f++){
				nbDat.push_back(0);
				FilterValue[f] = _MySystem->getAtomType(f);
			}
			FilteringType = "element";
		}
		for(unsigned int i=0;i<nbDatTot;i++) nbDat[_MySystem->getAtom(i).type_uint-1]++;
		nbDatMax = MT->max_vec(nbDat);
		FilterIndex = new unsigned int[nbDatMax*nbFilter];
		for(unsigned int f=0;f<nbFilter;f++) nbDat[f] = 0;
		unsigned int current_f;
		for(unsigned int i=0;i<nbDatTot;i++){
			current_f = _MySystem->getAtom(i).type_uint-1; // TODO modify here when changing the thing with type_uint
			FilterIndex[current_f*nbDatMax+nbDat[current_f]] = i;
			nbDat[current_f]++;
		}
	}else{
		unsigned int aux_size_filter, aux_id_filter;
		aux_id_filter = _MySystem->getAuxIdAndSize(FilteringType,aux_size_filter);
		if( aux_size_filter == 0 ){
			cerr << "The fitlering type \"" << FilteringType << "\" cannot be found in the auxiliary properties of the atomic system, aborting" << endl;
			exit(EXIT_FAILURE);
		}else if( aux_size_filter > 1 ){
			cerr << "The filtering type \"" << FilteringType << "\" is based on a vector, which is not implemented yet, aborting" << endl;
			exit(EXIT_FAILURE);
		}
		for(unsigned int i=0;i<nbDatTot;i++){
			bool already = false;
			for(unsigned int f=0;f<FilterValue.size();f++){
				if( to_string(_MySystem->getAux(aux_id_filter)[i]) == FilterValue[f] ){
					nbDat[f]++;
					already = true;
					break;
				}
			}
			if( !already ){
				nbDat.push_back(1);
				FilterValue.push_back(to_string(_MySystem->getAux(aux_id_filter)[i]));
			}
		}
		nbFilter = FilterValue.size();
		if( nbFilter > nbMaxFilter ) cout << "Warning there is " << nbFilter << " different filters, which is a lot" << endl;
		nbDatMax = MT->max_vec(nbDat);
		FilterIndex = new unsigned int[nbDatMax*nbFilter];
		for(unsigned int f=0;f<nbFilter;f++) nbDat[f] = 0;
		unsigned int current_f;
		for(unsigned int i=0;i<nbDatTot;i++){
			for(unsigned int f=0;f<FilterValue.size();f++){
				if( to_string(_MySystem->getAux(aux_id_filter)[i]) == FilterValue[f] ){
					FilterIndex[f*nbDatMax+nbDat[f]] = i;
					nbDat[f]++;
					break;
				}
			}
		}
	}
	// TODO verbose
	//if( nbFilter < nbMaxFilter ){
	//	cout << "The computed descriptors will be filtered by " << FilteringType << ", with :" << endl;
	//	for(unsigned int f=0;f<nbFilter;f++) cout << nbDat[f] << " descriptor having filter value: " << FilterValue[f] << endl;
	//}
}

Descriptors::Descriptors(const string &FilenameOrDir){ // This constructor read unlabelled descriptors from file or labelled descriptors from a directory containing subdirectories having the name of the label and containing files with the descriptors. This constructor considere that the data are not filtered
	// Test if the provided string is a directory or not
	string wd = filesystem::current_path();
	string full_path = wd+'/'+FilenameOrDir;
	string buffer_s, beg;
	struct dirent *diread;
	const char *env = full_path.c_str();
	DIR *dir;
	if( (dir = opendir(env) ) != nullptr ){
		cout << "Reading labelled descriptors from directory " << FilenameOrDir << endl;
		ifstream prop_file(full_path+"/DescriptorProperties.ath");
		if(prop_file){
			string line;
			cout << "Reading descriptor properties from " << FilenameOrDir+"/DescriptorProperties.ath" << endl;
			size_t pos;
			while(getline(prop_file,line)){
				pos = line.find("DESCRIPTOR_NAME");
				if(pos!=string::npos){
					istringstream text(line);
					text >> buffer_s;
					text >> name;
				}else Properties.push_back(line);
			}
			if( name == "none" ) cout << "Warning, the descriptor name has not been read from the DescriptorProperties.ath file" << endl;
			prop_file.close();
		}else{
			cout << "No file describing the descriptor properties is present in " << full_path << endl;
			cout << "The descriptor properties will be set to null" << endl;
			cout << "If you are fitting and labelling a machine learning model to put in the AtomHIC database, it is strongly encouraged to describe the descriptor properties by generating a DescriptorProperties.ath file in " << full_path << " directory" << endl;
			cout << "An example of such file can be found in /data/ExampleFiles/DescriptorProperties.ath" << endl;
		}

		while( (diread = readdir(dir)) != nullptr ){
			buffer_s = diread->d_name;
			size_t pos = buffer_s.find(".");
			if( pos==string::npos ) Labels.push_back(buffer_s);
			//if( buffer_s != "." && buffer_s != ".." ) Labels.push_back(buffer_s);
		}
		cout << Labels.size() << " different labels : ";
		for(unsigned int l=0;l<Labels.size();l++) cout << Labels[l] << " ";
		cout << endl;
		closedir(dir);
		vector<vector<string>> full_filenames;
		// list all files contained in the directories
		for(unsigned int l=0;l<Labels.size();l++){
			full_filenames.push_back(vector<string>());
			string full_path_dir = full_path+'/'+Labels[l];
			env = full_path_dir.c_str();
			if( (dir = opendir(env) ) != nullptr ){
				while( (diread = readdir(dir)) != nullptr ){
					buffer_s = diread->d_name;
					beg = buffer_s.substr(0,1); 
					if( beg != "." && buffer_s != "." && buffer_s != ".." ) full_filenames[l].push_back(full_path_dir+'/'+buffer_s);
				}
				closedir(dir);
			}
		}
		// first read of file to know the number of data point and dimension
		unsigned int count_dat(0), count_dim(0);
		string line;
		for(unsigned int l1=0;l1<Labels.size();l1++){
			for(unsigned int l2=0;l2<full_filenames[l1].size();l2++){
				unsigned int subcount_dat(0), subcount_dim(0);
				ifstream file(full_filenames[l1][l2]);
				if(file){
					string line;
					while(getline(file,line)){
						if( subcount_dat == 0 ){
							istringstream text(line);
							do{
								string sub;
								text >> sub;
								if( sub.length() ) subcount_dim++;
							}while( text );
						}
						subcount_dat++;
					}
					file.close();
				}else{
					cerr << "The file cannot be openned" << endl;
					exit(EXIT_FAILURE);
				}
				if( l1 == 0 && l2 == 0 ) count_dim = subcount_dim;
				if( subcount_dim != count_dim ){
					cerr << "The dimension of the provided descriptor files are different, especially for file " << full_filenames[l1][l2] << ", aborting" << endl;
					exit(EXIT_FAILURE);
				}
				count_dat += subcount_dat;
			}
		}
		dim = count_dim;
		nbDat.push_back(count_dat);
		nbDatMax = nbDat[0];
		_Descriptors = new double[dim*nbDat[0]];
		Labels_uint = new unsigned int[nbDat[0]];
		LabelsSize = new unsigned int[Labels.size()];
		nbFilter = 1;
		FilterIndex = new unsigned int[nbDatMax];
		for(unsigned int i=0;i<nbDatMax;i++) FilterIndex[i] = i;
		FilterValue.push_back("none");


		count_dat = 0;
		for(unsigned int l1=0;l1<Labels.size();l1++){
			LabelsSize[l1] = 0;
			for(unsigned int l2=0;l2<full_filenames[l1].size();l2++){
				unsigned int subcount_dat(0), subcount_dim(0);
				ifstream file(full_filenames[l1][l2]);
				if(file){
					string line;
					while(getline(file,line)){
						istringstream text(line);
						for(unsigned int i=0;i<dim;i++) text >> _Descriptors[count_dat*dim+i];
						Labels_uint[count_dat] = l1;
						LabelsSize[l1]++;
						count_dat++;
					}
					file.close();
				}else{
					cerr << "The file cannot be openned" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}
		cout << nbDat[0] << " descriptors of dimension " << dim << ", successfully read from directory " << FilenameOrDir << endl;
	}else{
		cout << "Reading descriptors from file : " << FilenameOrDir << endl;
		ifstream file(FilenameOrDir, ios::in);
		unsigned int count_line(0), count_dim(0);
		if(file){
			string line;
			while(getline(file,line)){
				if( count_line == 0 ){
					istringstream text(line);
					do{
						string sub;
						text >> sub;
						if( sub.length() ) count_dim++;
					}while( text );
				}
				count_line++;
			}
			file.close();
		}else{
			cerr << "The file cannot be openned" << endl;
			exit(EXIT_FAILURE);
		}

		dim = count_dim;
		nbDat.push_back(count_line);
		nbDatMax = nbDat[0];
		_Descriptors = new double[dim*nbDat[0]];
		FilterIndex = new unsigned int[nbDatMax];
		for(unsigned int i=0;i<nbDatMax;i++) FilterIndex[i] = i;
		FilterValue.push_back("none");
		nbFilter = 1;

		ifstream file2(FilenameOrDir, ios::in);
		if(file2){
			string line;
			unsigned int count(0);
			while(getline(file2,line)){
				istringstream text(line);
				for(unsigned int i=0;i<dim;i++) text >> _Descriptors[count*dim+i];
				count++;
			}
			file2.close();
		}else{
			cerr << "The file cannot be openned" << endl;
			exit(EXIT_FAILURE);
		}
		cout << nbDat[0] << " descriptors of dimension " << dim << ", successfully read from file " << FilenameOrDir << endl;
	}
	this->MT = new MathTools();
}

Descriptors::Descriptors(const string &FilenameOrDir, const string &DescriptorName){ // This constructor read labelled descriptors from atomic file or labelled descriptors from a directory containing subdirectories having the name of the label and containing files with the descriptors
	// Test if the provided string is a directory or not
	readFixedParams();
	string wd = filesystem::current_path();
	string full_path = wd+'/'+FilenameOrDir;
	string buffer_s, beg;
	struct dirent *diread;
	const char *env = full_path.c_str();
	DIR *dir;
	if( (dir = opendir(env) ) != nullptr ){
		cout << "Reading labelled descriptors from directory " << FilenameOrDir << endl;
		ifstream prop_file(full_path+"/DescriptorProperties.ath");
		if(prop_file){
			string line;
			cout << "Reading descriptor properties from " << FilenameOrDir+"/DescriptorProperties.ath" << endl;
			size_t pos;
			while(getline(prop_file,line)){
				pos = line.find("DESCRIPTOR_NAME");
				if(pos!=string::npos){
					istringstream text(line);
					text >> buffer_s;
					text >> name;
				}else Properties.push_back(line);
			}
			if( name == "none" ) cout << "Warning, the descriptor name has not been read from the DescriptorProperties.ath file" << endl;
			prop_file.close();
		}else{
			cout << "No file describing the descriptor properties is present in " << full_path << endl;
			cout << "The descriptor properties will be set to null" << endl;
			cout << "If you are fitting and labelling a machine learning model to put in the AtomHIC database, it is strongly encouraged to describe the descriptor properties by generating a DescriptorProperties.ath file in " << full_path << " directory" << endl;
			cout << "An example of such file can be found in /data/ExampleFiles/DescriptorProperties.ath" << endl;
		}
		while( (diread = readdir(dir)) != nullptr ){
			buffer_s = diread->d_name;
			size_t pos = buffer_s.find(".");
			if( pos==string::npos ) Labels.push_back(buffer_s);
		}
		cout << Labels.size() << " different labels : ";
		for(unsigned int l=0;l<Labels.size();l++) cout << Labels[l] << " ";
		cout << endl;
		closedir(dir);
		vector<vector<string>> full_filenames;
		// list all files contained in the directories
		for(unsigned int l=0;l<Labels.size();l++){
			full_filenames.push_back(vector<string>());
			string full_path_dir = full_path+'/'+Labels[l];
			env = full_path_dir.c_str();
			if( (dir = opendir(env) ) != nullptr ){
				while( (diread = readdir(dir)) != nullptr ){
					buffer_s = diread->d_name;
					beg = buffer_s.substr(0,1); 
					if( beg != "." && buffer_s != "." && buffer_s != ".." ) full_filenames[l].push_back(full_path_dir+'/'+buffer_s);
				}
				closedir(dir);
			}
		}
		// first read of file to know the number of data point and dimension
		bool IsFiltered = false;
		bool Type2Element = false;
		if( FilteringType != "none" ) IsFiltered = true;
		unsigned int count_dat(0), count_dim(0), line_At(1000), line_aux(1000), ReadOk, index_des, index_fil(100000), nbAux;
		size_t pos_At, pos_Aux, pos_aux_vec;
		string line, aux_name;
		for(unsigned int l1=0;l1<Labels.size();l1++){
			for(unsigned int l2=0;l2<full_filenames[l1].size();l2++){
				unsigned int subcount_dat(0), subcount_dim(0), subcount_filter(0);
				ReadOk = 0;
				ifstream file(full_filenames[l1][l2]);
				unsigned int count(0), subcount(0);
				if(file){
					string line;
					while(getline(file,line)){
						pos_At=line.find("NUMBER OF ATOMS");
						if( pos_At!=string::npos ) line_At = count;
						if( count == line_At+1 ){
							istringstream text(line);
							text >> subcount_dat;
							ReadOk++;
						}
						pos_Aux=line.find("ITEM: ATOMS");
						if( pos_Aux!=string::npos ){
							line_aux = count;
							istringstream text(line);
							while( text >> buffer_s ){
								pos_aux_vec = buffer_s.find("[");
								if( pos_aux_vec!=string::npos )	aux_name=buffer_s.substr(0,pos_aux_vec);
								else aux_name=buffer_s;
								if( aux_name == DescriptorName ){
									if( subcount_dim == 0 ) index_des = subcount-2;
									subcount_dim++;
								}else if( aux_name == FilteringType ){
									if( subcount_filter == 0 ) index_fil = subcount-2;
									subcount_filter++;
								}
								subcount++;
							}
							nbAux = subcount - 1;
							if( IsFiltered ){
								if( subcount_filter > 1 ){
									cerr << "The filtering type is based on a vector, which is not implemented yet, aborting" << endl;
									exit(EXIT_FAILURE);
								}else if( subcount_filter == 0 ){
									if( FilteringType == "element" ){
										istringstream text(line);
										unsigned int subsub = 0;
										Type2Element = false;
										while( text >> buffer_s ){
											if( buffer_s == "type" ){
												index_fil = subsub-2;
												Type2Element = true;
												break;
											}
											subsub++;
										}
										if( !Type2Element ){
											cerr << "Neither element nor type can be read from atomic file for filtering data, aborting" << endl;
											cerr << "If you don't want to filter the data you can set the DESCRIPTORS_FILTERING_TYPE parameter to none in /data/FixedParameters/FixedParameters.dat" << endl;
											exit(EXIT_FAILURE);
										}
									}else{
										cerr << "The filtering value \"" << FilteringType << "\" has not been found in the file \"" << full_filenames[l1][l2] << "\", aborting" << endl;
										exit(EXIT_FAILURE);
									}
								}
							}
							ReadOk++;
						}
						if( IsFiltered && count > line_aux ){
							istringstream text(line);
							for(unsigned int i=0;i<index_fil+1;i++) text >> buffer_s;
							bool already = false;
							for(unsigned int f=0;f<nbFilter;f++){
								if( buffer_s == FilterValue[f] ){
									nbDat[f]++;
									already = true;
									break;
								}
							}
							if( !already ){
								FilterValue.push_back(buffer_s);
								nbDat.push_back(1);
								nbFilter++;
								if( nbFilter >= nbMaxFilter ) cout << "Warning there is " << nbFilter << " different filters, which is a lot" << endl;
							}
						}
						count++;
					}
					file.close();
				}else{
					cerr << "The file cannot be openned" << endl;
					exit(EXIT_FAILURE);
				}
				if( ReadOk < 2 ){
					cerr << "The file format is not supported, aborting" << endl;
					cerr << "Try to generate and write your data using AtomHIC" << endl;
					exit(EXIT_FAILURE);
				}
				if( l1 == 0 && l2 == 0 ) count_dim = subcount_dim;
				if( count_dim != subcount_dim ){
					cerr << "The dimension of the provided descriptor files are different, especially for file " << full_filenames[l1][l2] << ", aborting" << endl;
					exit(EXIT_FAILURE);
				}
				count_dat += subcount_dat;
			}
		}
		

		dim = count_dim;
		if( !IsFiltered ){
			nbDat.push_back(count_dat);
			nbFilter = 1;
		}
		MT = new MathTools();
		nbDatMax = MT->max_vec(nbDat);
		nbDatTot = count_dat;
		FilterIndex = new unsigned int[nbFilter*nbDatMax];
		_Descriptors = new double[dim*nbDatTot];
		Labels_uint = new unsigned int[nbFilter*nbDatMax];
		LabelsSize = new unsigned int[Labels.size()*nbFilter];

		unsigned int *count_fil = new unsigned int[nbFilter];
		double *temp_des = new double[dim];
		for(unsigned int f=0;f<nbFilter;f++){
			count_fil[f] = 0;
			for(unsigned int l=0;l<Labels.size();l++) LabelsSize[f*Labels.size()+l] = 0;
		}
	
		count_dat = 0;
		for(unsigned int l1=0;l1<Labels.size();l1++){
			for(unsigned int l2=0;l2<full_filenames[l1].size();l2++){
				size_t pos_aux;
				unsigned int count(0);
				line_aux = 1000;
				ifstream file(full_filenames[l1][l2].c_str(), ifstream::in);
				if(file){
					string line;
					while(getline(file,line)){
						pos_Aux=line.find("ITEM: ATOMS");
						if( pos_Aux!=string::npos ) line_aux = count;
						if( count > line_aux ){
							istringstream text(line);
							unsigned int current_fil = 0;
							for(unsigned int ind=0;ind<nbAux;ind++){
								if( ind >= index_des && ind<index_des+dim ) text >> _Descriptors[count_dat*dim+ind-index_des];
								else if( ind == index_fil ){
									text >> buffer_s;
									for(unsigned int f=0;f<nbFilter;f++){
										if( buffer_s == FilterValue[f] ){
											FilterIndex[f*nbDatMax+count_fil[f]] = count_dat;
											current_fil = f;
											count_fil[f]++;
											break;
										}
									}
								}else text >> buffer_s;
							}
							Labels_uint[current_fil*nbDatMax+count_fil[current_fil]-1] = l1;
							LabelsSize[current_fil*Labels.size()+l1]++;
							count_dat++;
						}
						count++;
					}
				}else{
					cerr << "The file cannot be openned" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}

		vector<unsigned int> type_t2e;
		vector<string> element_t2e;
		unsigned int buffer_i;
		if( FilteringType == "element" && Type2Element ){
			cout << "The elements are not provided in the atomic systems" << endl;
			ifstream file_e2t("Type2Element.ath", ifstream::in);
			if( file_e2t ){
				cout << "Reading the correspondance between type and element in Type2Element.ath" << endl;
				while(getline(file_e2t,line)){
					istringstream text(line);
					text >> buffer_i >> buffer_s;
					type_t2e.push_back(buffer_i);
					element_t2e.push_back(buffer_s);
				}
				for(unsigned int f1=0;f1<nbFilter;f1++){
					for(unsigned int f2=0;f2<type_t2e.size();f2++){
						if( FilterValue[f1] == to_string(type_t2e[f2]) ){
							FilterValue[f1] = element_t2e[f2];
							break;
						}
					}
				}
			}else{
				cout << "You can provide a Type2Element.ath file giving the correspondance between type and element in the working directory (an example is present in /data/ExampleFiles/)" << endl;
				cout << "As this file has not been found, we will simply filter the descriptors by type" << endl;
			}
		}
		cout << nbDatTot << " descriptors of dimension " << this->dim << ", successfully read from directory " << FilenameOrDir << endl;
		if( nbFilter > 1 ){
			cout << "The descriptors are filtered by " << FilteringType << ", with :" << endl;
			for(unsigned int f=0;f<nbFilter;f++){
				cout << nbDat[f] << " descriptor having filter value: " << FilterValue[f] << " among which:" << endl;
				for(unsigned int l=0;l<Labels.size();l++) cout << LabelsSize[f*Labels.size()+l] << " are labelled by " << Labels[l] << endl;
			}
		}else{
			cout << "The descriptor are not filtered" << endl;
			for(unsigned int l=0;l<Labels.size();l++) cout << LabelsSize[l] << " are labelled by " << Labels[l] << endl;
		}
		delete[] count_fil;
		delete[] temp_des;	
	}else{ 
		cerr << "Reading a single labelled atomic file is not currently implemented" << endl;
		cerr << "The best is to use the AtomicSystem constructor to build the descriptors" << endl;
		cerr << "Or to place your single file in directories to make labelled descriptors (the unique label wont then be used)" << endl;
		cerr << "For instance do: mkdir TempDir; mkdir TempDir/Label; mv " << FilenameOrDir << " TempDir/Label/.; and finally recall AtomHIC executable giving TempDir directory as argument" << endl;
		exit(EXIT_FAILURE);
	}
}
	
void Descriptors::printDescriptorsPropToDatabase(ofstream &writefile){
	setProperties();
	string buffer_s;
	bool IsDesName = false;
	bool Isftype = false;
	bool IsNbDim = false;
	size_t pos_d, pos_t, pos_ndim;
	for(unsigned int s=0;s<Properties.size();s++){
		pos_d = Properties[s].find("DESCRIPTOR_NAME");
		if( pos_d!=string::npos ) IsDesName = true;
		pos_t = Properties[s].find("FILTER_TYPE");
		if( pos_t!=string::npos ) Isftype = true;
		pos_ndim = Properties[s].find("NUMBER_OF_DIMENSION");
		if( pos_ndim!=string::npos ) IsNbDim = true;
	}
	if( !IsDesName ) writefile << "DESCRIPTOR_NAME " << name << endl;
	if( !Isftype ) writefile << "FILTER_TYPE " << FilteringType << endl;
	if( !IsNbDim ) writefile << "NUMBER_OF_DIMENSION " << dim << endl;
	for(unsigned int s=0;s<Properties.size();s++) writefile << Properties[s] << endl;
}

void Descriptors::readProperties(vector<string> _Properties){
	Properties.clear();
	for(unsigned int s=0;s<_Properties.size();s++) Properties.push_back(_Properties[s]);
	size_t pos_ftype;
	string buffer_s;
	for(unsigned int s=0;s<_Properties.size();s++){
		pos_ftype = _Properties[s].find("FILTER_TYPE");
		if( pos_ftype!=string::npos ){
			istringstream text(_Properties[s]);
			text >> buffer_s >> FilteringType;
			break;
		}
	}
}

void Descriptors::readFixedParams(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_filter;
	string buffer_s, line;
	if(file){
		while(file){
			getline(file,line);
			pos_filter=line.find("DESCRIPTORS_FILTERING_TYPE ");
			if(pos_filter!=string::npos){
				istringstream text(line);
				text >> buffer_s >> FilteringType;
			}
		}
	}else{
		cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
		exit(EXIT_FAILURE);
	}
	file.close();
	//cout << "From /data/FixedParameters/FixedParameters.dat the descriptors will be filtered by \"" << FilteringType << "\"" << endl;
}

Descriptors::~Descriptors(){
	delete[] FilterIndex;
	if( AreDescriptorsMine ) delete[] _Descriptors;
	delete MT;
}
