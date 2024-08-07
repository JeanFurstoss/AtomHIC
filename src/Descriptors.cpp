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
		while( (diread = readdir(dir)) != nullptr ){
			buffer_s = diread->d_name;
			if( buffer_s != "." && buffer_s != ".." ) Labels.push_back(buffer_s);
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
		nbFilter = 1;
		FilterValue.push_back("none");


		count_dat = 0;
		for(unsigned int l1=0;l1<Labels.size();l1++){
			for(unsigned int l2=0;l2<full_filenames[l1].size();l2++){
				unsigned int subcount_dat(0), subcount_dim(0);
				ifstream file(full_filenames[l1][l2]);
				if(file){
					string line;
					while(getline(file,line)){
						istringstream text(line);
						for(unsigned int i=0;i<dim;i++) text >> _Descriptors[count_dat*dim+i];
						Labels_uint[count_dat] = l1;
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

Descriptors::Descriptors(const string &FilenameOrDir, const string &DescriptorName, const string &_FilteringType){ // This constructor read labelled descriptors from atomic file or labelled descriptors from a directory containing subdirectories having the name of the label and containing files with the descriptors
	// Test if the provided string is a directory or not
	FilteringType = _FilteringType;
	string wd = filesystem::current_path();
	string full_path = wd+'/'+FilenameOrDir;
	string buffer_s, beg;
	struct dirent *diread;
	const char *env = full_path.c_str();
	DIR *dir;
	if( (dir = opendir(env) ) != nullptr ){
		cout << "Reading labelled descriptors from directory " << FilenameOrDir << endl;
		while( (diread = readdir(dir)) != nullptr ){
			buffer_s = diread->d_name;
			if( buffer_s != "." && buffer_s != ".." ) Labels.push_back(buffer_s);
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
									cerr << "The filtering value \"" << FilteringType << "\" has not been found in the file \"" << full_filenames[l1][l2] << "\", aborting" << endl;
									exit(EXIT_FAILURE);
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
								if( nbFilter >= nbMaxFilter ){
									cerr << "The number of filter is higher than maximum allowed, if you know what you are doing you can change its value in /data/FixedParameters/FixedParameters.dat" << endl;
									exit(EXIT_FAILURE);
								}
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
		_Descriptors = new double[dim*nbFilter*nbDatMax];
		Labels_uint = new unsigned int[nbFilter*nbDatMax];
		LabelsSize = new unsigned int[Labels.size()*nbFilter];

		unsigned int *count_fil = new unsigned int[nbFilter];
		double *temp_des = new double[dim];
		for(unsigned int f=0;f<nbFilter;f++){
			count_fil[f] = 0;
			for(unsigned int l=0;l<Labels.size();l++) LabelsSize[f*Labels.size()+l] = 0;
		}
	
		for(unsigned int l1=0;l1<Labels.size();l1++){
			for(unsigned int l2=0;l2<full_filenames[l1].size();l2++){
				size_t pos_aux;
				unsigned int count(0);
				line_aux = 1000;
				ifstream file(full_filenames[l1][l2]);
				if(file){
					string line;
					while(getline(file,line)){
						pos_Aux=line.find("ITEM: ATOMS");
						if( pos_Aux!=string::npos ) line_aux = count;
						if( count > line_aux ){
							istringstream text(line);
							unsigned int current_fil = 0;
							for(unsigned int ind=0;ind<nbAux;ind++){
								if( ind >= index_des && ind<=index_des+dim ) text >> temp_des[ind-index_des];
								else if( ind == index_fil ){
									text >> buffer_s;
									for(unsigned int f=0;f<nbFilter;f++){
										if( buffer_s == FilterValue[f] ){
											current_fil = f;
											break;
										}
									}
								}else text >> buffer_s;
							}
							for(unsigned int d=0;d<dim;d++) _Descriptors[current_fil*nbDatMax*dim+count_fil[current_fil]*dim+d] = temp_des[d];
							Labels_uint[current_fil*nbDatMax+count_fil[current_fil]] = l1;
							LabelsSize[current_fil*Labels.size()+l1]++;
							count_fil[current_fil]++;
						}
						count++;
					}
					file.close();
				}else{
					cerr << "The file cannot be openned" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}
		unsigned int nbDatTot = 0;
		for(unsigned int f=0;f<nbFilter;f++) nbDatTot += nbDat[f];
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
Descriptors::~Descriptors(){
	delete[] _Descriptors;
	delete MT;
}
