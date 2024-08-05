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

Descriptors::Descriptors(const string &FilenameOrDir){ // This constructor read unlabelled descriptors from file or labelled descriptors from a directory containing subdirectories having the name of the label and containing files with the descriptors
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

		this->dim = count_dim;
		this->nbDat = count_dat;
		this->_Descriptors = new double[this->dim*this->nbDat];
		Labels_uint = new unsigned int[nbDat];

		count_dat = 0;
		for(unsigned int l1=0;l1<Labels.size();l1++){
			for(unsigned int l2=0;l2<full_filenames[l1].size();l2++){
				unsigned int subcount_dat(0), subcount_dim(0);
				ifstream file(full_filenames[l1][l2]);
				if(file){
					string line;
					while(getline(file,line)){
						istringstream text(line);
						for(unsigned int i=0;i<this->dim;i++) text >> _Descriptors[count_dat*dim+i];
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
		cout << this->nbDat << " descriptors of dimension " << this->dim << ", successfully read from directory " << FilenameOrDir << endl;
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
		this->dim = count_dim;
		this->nbDat = count_line;
		this->_Descriptors = new double[this->dim*this->nbDat];

		ifstream file2(FilenameOrDir, ios::in);
		if(file2){
			string line;
			unsigned int count(0);
			while(getline(file2,line)){
				istringstream text(line);
				for(unsigned int i=0;i<this->dim;i++) text >> _Descriptors[count*dim+i];
				count++;
			}
			file2.close();
		}else{
			cerr << "The file cannot be openned" << endl;
			exit(EXIT_FAILURE);
		}
		cout << this->nbDat << " descriptors of dimension " << this->dim << ", successfully read from file " << FilenameOrDir << endl;
	}
	this->MT = new MathTools();
}

Descriptors::Descriptors(const string &FilenameOrDir, const string &DescriptorName){ // This constructor read labelled descriptors from atomic file or labelled descriptors from a directory containing subdirectories having the name of the label and containing files with the descriptors
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
		unsigned int count_dat(0), count_dim(0), line_At(1000), ReadOk, index_des;
		size_t pos_At, pos_Aux, pos_aux_vec;
		string line, aux_name;
		for(unsigned int l1=0;l1<Labels.size();l1++){
			for(unsigned int l2=0;l2<full_filenames[l1].size();l2++){
				unsigned int subcount_dat(0), subcount_dim(0);
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
							istringstream text(line);
							while( text >> buffer_s ){
								pos_aux_vec = buffer_s.find("[");
								if( pos_aux_vec!=string::npos )	aux_name=buffer_s.substr(0,pos_aux_vec);
								else aux_name=buffer_s;
								if( aux_name == DescriptorName ){
									if( subcount_dim == 0 ) index_des = subcount;
									subcount_dim++;
								}
								subcount++;
							}
							ReadOk++;
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
		
		index_des -= 2; //the ITEM: ATOM

		this->dim = count_dim;
		this->nbDat = count_dat;
		this->_Descriptors = new double[this->dim*this->nbDat];
		Labels_uint = new unsigned int[nbDat];
		count_dat = 0;
		for(unsigned int l1=0;l1<Labels.size();l1++){
			for(unsigned int l2=0;l2<full_filenames[l1].size();l2++){
				size_t pos_aux;
				unsigned int count(0);
				unsigned int line_aux(1000);
				ifstream file(full_filenames[l1][l2]);
				if(file){
					string line;
					while(getline(file,line)){
						pos_Aux=line.find("ITEM: ATOMS");
						if( pos_Aux!=string::npos ) line_aux = count;
						if( count > line_aux ){
							istringstream text(line);
							for(unsigned int ind=0;ind<index_des;ind++) text >> buffer_s;
							for(unsigned int ind=index_des;ind<index_des+dim+1;ind++) text >> _Descriptors[count_dat*dim+ind-index_des];
							Labels_uint[count_dat] = l1;
							count_dat++;
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
		cout << this->nbDat << " descriptors of dimension " << this->dim << ", successfully read from directory " << FilenameOrDir << endl;
	}else{ // TODO
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
		this->dim = count_dim;
		this->nbDat = count_line;
		this->_Descriptors = new double[this->dim*this->nbDat];

		ifstream file2(FilenameOrDir, ios::in);
		if(file2){
			string line;
			unsigned int count(0);
			while(getline(file2,line)){
				istringstream text(line);
				for(unsigned int i=0;i<this->dim;i++) text >> _Descriptors[count*dim+i];
				count++;
			}
			file2.close();
		}else{
			cerr << "The file cannot be openned" << endl;
			exit(EXIT_FAILURE);
		}
		cout << this->nbDat << " descriptors of dimension " << this->dim << ", successfully read from file " << FilenameOrDir << endl;
	}
	this->MT = new MathTools();
}
Descriptors::~Descriptors(){
	delete[] _Descriptors;
	delete MT;
}
