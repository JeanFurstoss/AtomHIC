//**********************************************************************************
//*   Descriptors.cpp                                                              *
//**********************************************************************************
//* This file contains the implementation of the Descriptors class                 *
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


#include "Descriptors.h"
#include "AtomHicConfig.h"
#include <iostream>
#include <dirent.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <filesystem>
#include <random>

using namespace std;


Descriptors::Descriptors(AtomicSystem *_MySystem):_MySystem(_MySystem){
	readFixedParams();
	MT = new MathTools();
	nbDatTot = _MySystem->getNbAtom();
	ConstructFilterIndexArray(_MySystem);
}

Descriptors::Descriptors(AtomicSystem *_MySystem, string DescriptorName, string _FilteringType):_MySystem(_MySystem){
	FilteringType = _FilteringType;
	MT = new MathTools();
	nbDatTot = _MySystem->getNbAtom();
	ConstructFilterIndexArray(_MySystem);
	if( DescriptorName == "Position" ){
		dim = 3;
		ConstructDescriptorFromAtomicPosition();
	}else{
		unsigned int aux_size, aux_id;
		aux_id = _MySystem->getAuxIdAndSize(DescriptorName,aux_size);
		if( aux_size == 0 ){
			cerr << "The provided descriptor name cannot be found in the auxiliary properties of the atomic system (it could also be the keyword \"Position\"), aborting" << endl;
			exit(EXIT_FAILURE);
		}
		dim = aux_size;
		_Descriptors = _MySystem->getAux(aux_id);
		AreDescriptorsMine = false;
	}
}

Descriptors::Descriptors(AtomicSystem *_MySystem, vector<string> _Properties):_MySystem(_MySystem){
	readFixedParams(); // read fixed parameters before in order to have all values initialized even if some are missing in _Properties
	readProperties(_Properties);
	MT = new MathTools();
	nbDatTot = _MySystem->getNbAtom();
	ConstructFilterIndexArray(_MySystem);
}

Descriptors::Descriptors(double *Descriptors, unsigned int nbdat, unsigned int dim){
	_Descriptors = Descriptors;
	MT = new MathTools();
	nbDat.push_back(nbdat);
	nbDatMax = nbDat[0];
	nbDatTot = nbDat[0];
	nbFilter = 1;
	FilterValue.push_back("none");
	FilterIndex = new unsigned int[nbDat[0]];
	for(unsigned int i=0;i<nbDat[0];i++) FilterIndex[i] = i;
	this->dim = dim;
	AreDescriptorsMine = false;

}

void Descriptors::ConstructLabelIndexArray(){
	if( nbLabels == 0 ){
		cerr << "The descriptors are not labelled, cannot construct LabelIndex array" << endl;
		exit(EXIT_FAILURE);
	}
	nbDatMaxLabel = MT->max_p(LabelsSize,nbLabels*nbFilter);
	LabelIndex = new unsigned int[nbDatMaxLabel*nbLabels*nbFilter];
	DescriptorLabel = new unsigned int[nbDatTot];
	unsigned int *count_dat = new unsigned int[nbLabels*nbFilter];
	for(unsigned int i=0;i<nbLabels;i++){
		for(unsigned int j=0;j<nbFilter;j++) count_dat[i*nbFilter+j] = 0;
	}
	for(unsigned int f=0;f<nbFilter;f++){
		for(unsigned int i=0;i<nbDat[f];i++){
			unsigned int current_label = Labels_uint[f*nbDatMax+i];
			DescriptorLabel[FilterIndex[f*nbDatMax+i]] = current_label;
			LabelIndex[f*nbLabels*nbDatMaxLabel+current_label*nbDatMaxLabel+count_dat[current_label*nbFilter+f]] = FilterIndex[f*nbDatMax+i];
			count_dat[current_label*nbFilter+f] += 1;
		}
	}
}

void Descriptors::ConstructFilterIndexArray(AtomicSystem *_MySystem){
	FilterValue.clear();
	if( FilteringType == "none" ){
		nbFilter = 1;
		FilterValue.push_back("none");
		nbDat.push_back(_MySystem->getNbAtom());
		nbDatMax = nbDat[0];
		nbDatTot = nbDat[0];
		FilterIndex = new unsigned int[nbDat[0]];
		DescriptorFilter = new unsigned int[nbDat[0]];
		for(unsigned int i=0;i<nbDat[0];i++){
			FilterIndex[i] = i;
			DescriptorFilter[i] = 0;
		}
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
				for(unsigned int f=0;f<nbFilter;f++){
					nbDat.push_back(0);
					FilterValue[f] = f+1;
				}
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
		DescriptorFilter = new unsigned int[nbDatTot];
		for(unsigned int f=0;f<nbFilter;f++) nbDat[f] = 0;
		unsigned int current_f;
		for(unsigned int i=0;i<nbDatTot;i++){
			current_f = _MySystem->getAtom(i).type_uint-1; // TODO modify here when changing the thing with type_uint
			FilterIndex[current_f*nbDatMax+nbDat[current_f]] = i;
			DescriptorFilter[i] = current_f;
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
		DescriptorFilter = new unsigned int[nbDatTot];
		for(unsigned int f=0;f<nbFilter;f++) nbDat[f] = 0;
		unsigned int current_f;
		for(unsigned int i=0;i<nbDatTot;i++){
			for(unsigned int f=0;f<FilterValue.size();f++){
				if( to_string(_MySystem->getAux(aux_id_filter)[i]) == FilterValue[f] ){
					FilterIndex[f*nbDatMax+nbDat[f]] = i;
					DescriptorFilter[i] = f;
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
	this->MT = new MathTools();
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
		nbDatTot = nbDat[0];
		_Descriptors = new double[dim*nbDat[0]];
		Labels_uint = new unsigned int[nbDat[0]];
		LabelsSize = new unsigned int[Labels.size()];
		nbLabels = Labels.size();
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
		ConstructLabelIndexArray();
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
		nbDatTot = nbDat[0];
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
		nbLabels = Labels.size();
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
		ConstructLabelIndexArray();
	}else{ 
		cerr << "Reading a single labelled atomic file is not currently implemented" << endl;
		cerr << "The best is to use the AtomicSystem constructor to build the descriptors" << endl;
		cerr << "Or to place your single file in directories to make labelled descriptors (the unique label wont then be used)" << endl;
		cerr << "For instance do: mkdir TempDir; mkdir TempDir/Label; mv " << FilenameOrDir << " TempDir/Label/.; and finally recall AtomHIC executable giving TempDir directory as argument" << endl;
		exit(EXIT_FAILURE);
	}
}

void Descriptors::ConstructDescriptorFromAtomicPosition(){
	// TODO better thing when struct Atom has changed (i.e. not a copy)
	_Descriptors = new double[dim*nbDatTot];
	for(unsigned int i=0;i<nbDatTot;i++){
		_Descriptors[i*3] = _MySystem->getAtom(i).pos.x;
		_Descriptors[i*3+1] = _MySystem->getAtom(i).pos.y;
		_Descriptors[i*3+2] = _MySystem->getAtom(i).pos.z;
	}
}

void Descriptors::ComputeExtremums(){
	min_val = new double[nbFilter*dim];
	max_val = new double[nbFilter*dim];
	for(unsigned int f=0;f<nbFilter;f++){
		for(unsigned int d=0;d<dim;d++){
			unsigned int min_ind = 0;
			unsigned int max_ind = 0;
			for(unsigned int i=1;i<nbDat[f];i++){
				if( _Descriptors[FilterIndex[f*nbDatMax+i]*dim+d] > _Descriptors[FilterIndex[f*nbDatMax+max_ind]*dim+d] ) max_ind = i;
				if( _Descriptors[FilterIndex[f*nbDatMax+i]*dim+d] < _Descriptors[FilterIndex[f*nbDatMax+min_ind]*dim+d] ) min_ind = i;
			}
			min_val[f*dim+d] = _Descriptors[FilterIndex[f*nbDatMax+min_ind]*dim+d];
			max_val[f*dim+d] = _Descriptors[FilterIndex[f*nbDatMax+max_ind]*dim+d];
		}
	}
	AreExtremums = true;
}

unsigned int Descriptors::getNbMaxNAndFilter(string filter_name, unsigned int &f){
	bool already = false;
	for(unsigned int ft=0;ft<Neigh_FilterValue.size();ft++){
		if( filter_name == Neigh_FilterValue[ft] ){
			already = true;
			f = ft;
			break;
		}
	}
	if( !already ){
		cerr << "Neighbour list has not been computed for filter \"" << filter_name << "\", aborting" << endl;
		exit(EXIT_FAILURE);
		return 0;
	}

	return this->nbMaxN[f];
}

bool Descriptors::getIsNeighbours(string filter_name){ 
	bool already = false;
	unsigned int f;
	for(unsigned int ft=0;ft<Neigh_FilterValue.size();ft++){
		if( filter_name == Neigh_FilterValue[ft] ){
			already = true;
			break;
		}
	}
	return already;
}

double Descriptors::get_current_rc(std::string filter_name){
	bool already = false;
	unsigned int f;
	for(unsigned int ft=0;ft<Neigh_FilterValue.size();ft++){
		if( filter_name == Neigh_FilterValue[ft] ){
			already = true;
			f = ft;
			break;
		}
	}
	if( !already ) return 0.;
	else return rc_neighbours[f];
}

double Descriptors::ComputeAverageMinDistance(string FilterVal, string distFunc){
	unsigned int current_f = current_filter(FilterVal);
	// search extremums for computing the apparent volume of the system 
	if( !AreExtremums ) ComputeExtremums();
	double AppVol = 1.;
	for(unsigned int d=0;d<dim;d++){
		double length = max_val[current_f*dim+d]-min_val[current_f*dim+d];
		if( length == 0 ){
			length = 1;
			cout << "Warning, the dimension " << d+1 << " of descriptors is flat" << endl;
			//cerr << "The dimension " << d+1 << " of descriptors is flat, cannot search neighbours, aborting" << endl;
			//exit(EXIT_FAILURE);
		}
		AppVol *= length;
	}
	double Density = AppVol / nbDat[current_f];
	double secFac = 2.; // TODO test this value
	double rc = secFac*pow(Density,1./dim);
	searchNeighbours(rc,FilterVal);
	unsigned int current_f_neigh;
        unsigned int current_nbMaxN = getNbMaxNAndFilter(FilterVal,current_f_neigh);
	double ave_dist = 0.;
	for(unsigned int i=0;i<nbDat[current_f];i++){
		vector<double> dist_q;
		unsigned int curid = FilterIndex[current_f*nbDatMax+i];
		for(unsigned int j=0;j<Neighbours[current_f_neigh][i*(current_nbMaxN+1)];j++){
			unsigned int id = FilterIndex[current_f*nbDatMax+Neighbours[current_f_neigh][i*(current_nbMaxN+1)+j+1]];
			double curdist = 0.;
			for(unsigned int d=0;d<dim;d++) curdist += (_Descriptors[curid*dim+d] - _Descriptors[id*dim+d]) * (_Descriptors[curid*dim+d] - _Descriptors[id*dim+d]);
			dist_q.push_back(curdist);
		}
		ave_dist += sqrt(MT->min_vec(dist_q));
	}
	return ave_dist / nbDat[current_f];
}

void Descriptors::searchNeighbours(double rc, string FilterVal, string distFunc){
	unsigned int MaxNbCellPerDim = 10;
	unsigned int current_f = current_filter(FilterVal);
	bool already = false;
	bool already_filter = false;
	unsigned int current_ind_neigh = 0;
	for(unsigned int f=0;f<Neigh_FilterValue.size();f++){
		if( FilterVal == Neigh_FilterValue[f] ){
			already_filter = true;
			current_ind_neigh = f;
		        if( rc_neighbours[f] == rc ){
				already = true;
				break;
			}
		}
	}
	if( already ){
		cout << "Neighbour list has already been computed for this filter and this cutoff" << endl;
		return;
	}
	if( distFunc != "Euclidian" ){ // TODO maybe this ceil list algo also work for other type of distances (should have a look on the maths of distances)
		cerr << "For the moment other distance function than Euclidian are not implemented, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	if( !already_filter ){
		Neigh_FilterValue.push_back(FilterVal);
		current_ind_neigh = Neigh_FilterValue.size()-1;
	}
	double rc_squared = rc*rc;
	// search extremums for computing the ceils 
	if( !AreExtremums ) ComputeExtremums();
	// ceil list algo
	unsigned int *nbCells = new unsigned int[dim];
	double *CellSizes = new double[dim];
	unsigned int nbCellTot = 0;
	double length = max_val[current_f*dim]-min_val[current_f*dim];
	if( length == 0 ){
		length = 1;
		cout << "Warning the dimension 1 of descriptors is flat" << endl;
		//cerr << "The dimension 1 of descriptors is flat, cannot search neighbours, aborting" << endl;
		//exit(EXIT_FAILURE);
	}
	double CellVolume;
	nbCells[0] = floor(length/rc);
	if( nbCells[0] > MaxNbCellPerDim ) nbCells[0] = MaxNbCellPerDim;
	if( nbCells[0] == 0 ){
		nbCells[0] = 1;
		CellSizes[0] = length;
	}else CellSizes[0] = length / (double) nbCells[0];
	nbCellTot = nbCells[0];
	CellVolume = CellSizes[0];
	for(unsigned int d=1;d<dim;d++){
		double length = max_val[current_f*dim+d]-min_val[current_f*dim+d];
		if( length == 0 ){
			length = 1;
			cout << "Warning, the dimension " << d+1 << " of descriptors is flat" << endl;
			//cerr << "The dimension " << d+1 << " of descriptors is flat, cannot search neighbours, aborting" << endl;
			//exit(EXIT_FAILURE);
		}
		nbCells[d] = floor(length/rc);
		if( nbCells[d] > MaxNbCellPerDim ) nbCells[d] = MaxNbCellPerDim;
		if( nbCells[d] == 0 ){
			nbCells[d] = 1;
			CellSizes[d] = length;
		}else CellSizes[d] = length / (double) nbCells[d];
		nbCellTot *= nbCells[d];
		CellVolume *= CellSizes[d];
	}
	// indexes for accessing cells
	unsigned int *CellIndexes = new unsigned int[dim];
	CellIndexes[dim-1] = 1;
	for(int d1=dim-2;d1>=0;d1--){
		CellIndexes[d1] = nbCells[dim-1];
		for(int d2=dim-2;d2>=d1+1;d2--){
			CellIndexes[d1] *= nbCells[d2];
		}
	}
	// Fill cells
	vector<vector<unsigned int>> Cells(nbCellTot); // Cells[ i1*(prod nbCell(d!=d1)) + i2*(prod nbCell(d!=d1 && d!=d2)) + .. + iN ][i] = id of ith point belonging to Cell (i1,i2,i3,..,iN)
	//#pragma omp parallel for
	for(unsigned int i=0;i<nbDat[current_f];i++){
		unsigned int current_cell_index = 0;
		unsigned int current_index = FilterIndex[current_f*nbDatMax+i];
		for(unsigned int d=0;d<dim;d++){ // apply homotethy to retrieve the cell index
			unsigned int tmp = floor( ( _Descriptors[current_index*dim+d] - min_val[current_f*dim+d] ) / CellSizes[d]);
			if( tmp == nbCells[d] )	tmp--;
			current_cell_index += tmp*CellIndexes[d];
		}
		Cells[current_cell_index].push_back(i);
	}
	// compute nbMaxN
	unsigned int maxN, nbtotverif;
	maxN = Cells[0].size();
	nbtotverif = maxN;
	for(unsigned int n=1;n<nbCellTot;n++){
		unsigned int current_size = Cells[n].size();
		nbtotverif += current_size;
		if( current_size > maxN ) maxN = current_size;
	}
	if( nbtotverif != nbDat[current_f] ){
	        cerr << "Some descriptor points were lost during construction of ceils for neighbour research, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	// Compute hypersphere volume
	double HypSphVol, tmp(1.);
	if( dim%2 == 0 ){
		for(unsigned int d=1;d<=dim/2;d++) tmp *= d;
		HypSphVol = 1. / tmp;
		HypSphVol *= pow(M_PI,(double) dim/2.);
	}else{
		for(unsigned int d=1;d<=dim;d+=2) tmp *= d;
		HypSphVol = 1. / tmp;
		HypSphVol *= pow(M_PI,((double) dim-1.)/2.);
		HypSphVol *= pow(2.,((double) dim+1.)/2.);
	}
	HypSphVol *= pow(rc,dim);
	double secFac = 2.;
	if( (unsigned int) (maxN*secFac*HypSphVol/CellVolume) == 0 ) maxN = (unsigned int) (maxN*secFac);
	else maxN = (unsigned int) (maxN*secFac*HypSphVol/CellVolume); 
	if( !already_filter ){
		Neighbours.push_back(new unsigned int[(maxN+1)*nbDat[current_f]]);
		nbMaxN.push_back(maxN);
		rc_neighbours.push_back(rc);
	}else{
		delete[] Neighbours[current_ind_neigh];
		Neighbours[current_ind_neigh] = new unsigned int[(maxN+1)*nbDat[current_f]];
		nbMaxN[current_ind_neigh] = maxN;
		rc_neighbours[current_ind_neigh] = rc;
	}
	// search neighbours
	vector<vector<int>> combinations = MT->GenerateNDCombinations(dim,-1,1);
	bool allok = false;
	while( !allok ){
		allok = true;
		//#pragma omp parallel for // cannot parallelized for the moment due to the break statement..
		for(unsigned int n=0;n<nbCellTot;n++){
			Dis.ProgressBar(nbCellTot,n);
			double squared_dist;
			// get the dimension index of the cell
			vector<unsigned int> current_cell_indexes(dim);
			vector<unsigned int> neigh_cell_indexes(dim);
			current_cell_indexes[0] = floor(n/CellIndexes[0]);
			unsigned int n_new = n - current_cell_indexes[0]*CellIndexes[0];
		       for(unsigned int d=1;d<dim;d++){
			       current_cell_indexes[d] = floor( n_new / CellIndexes[d]);
			       n_new -= current_cell_indexes[d]*CellIndexes[d];
		       }
			// search the neighbouring cells indexes
			vector<unsigned int> full_cell_id;
			for(unsigned int c=0;c<combinations.size();c++){
				bool totreat = true;
				for(unsigned int d=0;d<dim;d++){
					neigh_cell_indexes[d] = current_cell_indexes[d] + combinations[c][d];
					if( neigh_cell_indexes[d] < 0 || neigh_cell_indexes[d] >= nbCells[d] ){ // cell outside of bond => do not treat
						totreat = false;
						break;
					}
				}
				if( !totreat ) continue;
				full_cell_id.push_back(0);
				unsigned int lastind = full_cell_id.size()-1;
				for(unsigned int d=0;d<dim;d++) full_cell_id[lastind] += neigh_cell_indexes[d]*CellIndexes[d];
			}

		       // search neighbours of the current cell
		       for(unsigned int i=0;i<Cells[n].size();i++){
		       		// initialize counter
		       		unsigned int count_neigh = 0;
				for(unsigned int c=0;c<full_cell_id.size();c++){
		       			for(unsigned int j=0;j<Cells[full_cell_id[c]].size();j++){
						if( Cells[n][i] != Cells[full_cell_id[c]][j] ){
							SquaredDistance(Cells[n][i],Cells[full_cell_id[c]][j],current_f,squared_dist);
			       				if( squared_dist < rc_squared ){
								Neighbours[current_ind_neigh][Cells[n][i]*(maxN+1)+1+count_neigh] = Cells[full_cell_id[c]][j];
								count_neigh++;
								if( count_neigh > maxN ){
									cout << "Maximum number of neighbours allowed overpassed, doubling this number and relaunch neighbour research" << endl;
									allok = false;
									break;
									//cerr << "Maximum number of neighbours allowed overpassed, aborting" << endl;
									//exit(EXIT_FAILURE);
								}
							}
						}
					}
					if( !allok ) break;
				}
				if( !allok ) break;
				Neighbours[current_ind_neigh][Cells[n][i]*(maxN+1)] = count_neigh;
		       }
		       if( !allok ) break;
		} // end loop on nbCellTot
		if( !allok ){
			maxN *= 2;
			delete[] Neighbours[current_ind_neigh];
			Neighbours[current_ind_neigh] = new unsigned int[(maxN+1)*nbDat[current_f]];
			nbMaxN[current_ind_neigh] = maxN;
		}else IsNeighbours = true;
	}
	delete[] nbCells;
	delete[] CellSizes;
	delete[] CellIndexes;
}

//double Descriptors::SquaredDistance(unsigned int id_1, unsigned int id_2, unsigned int filter){// For the moment Euclidian distance (to add other distance and rename)
void Descriptors::SquaredDistance(unsigned int id_1, unsigned int id_2, unsigned int filter, double &squared_dist){// For the moment Euclidian distance (to add other distance and rename)
	//double squared_dist = 0.;
	squared_dist = 0.;
	for(unsigned int d=0;d<dim;d++){
		double diff = _Descriptors[FilterIndex[filter*nbDatMax+id_1]*dim+d]-_Descriptors[FilterIndex[filter*nbDatMax+id_2]*dim+d];
		squared_dist += diff*diff;
	}
	//return squared_dist;
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

// TEST subarrays
void Descriptors::constructSubarrays(bool filter, bool label){
	// TODO delete if already allocated
	if( !filter && !label ){
		if( !subarray_defined ){
			DescriptorsSubarray.push_back(new MatrixXd(nbDatTot,dim)); 
			SubarraySize.push_back(nbDatTot);
			CorresIndexSubarray.push_back(new unsigned int[nbDatTot]);
			for(unsigned int i=0;i<nbDatTot;i++){
				CorresIndexSubarray[0][i] = i;
				for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[0]))(i,d) = _Descriptors[i*dim+d];
			}
			subarray_defined = true;
		}else if( subarray_filter || subarray_label ){
			for(unsigned int i=0;i<DescriptorsSubarray.size();i++){
				delete DescriptorsSubarray[0];// TODO test this
				DescriptorsSubarray.erase(DescriptorsSubarray.begin());
				SubarraySize.erase(SubarraySize.begin());
				delete[] CorresIndexSubarray[0];
				CorresIndexSubarray.erase(CorresIndexSubarray.begin());
			}
			DescriptorsSubarray.push_back(new MatrixXd(nbDatTot,dim));
			SubarraySize.push_back(nbDatTot);
			CorresIndexSubarray.push_back(new unsigned int[nbDatTot]);
			for(unsigned int i=0;i<nbDatTot;i++){
				CorresIndexSubarray[0][i] = i;
				for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[0]))(i,d) = _Descriptors[i*dim+d];
			}
		}
		subarray_filter = false;
		subarray_label = false;
	}else if( filter && !label ){
		if( !subarray_defined ){
			//cout << "Constructing descriptor subarrays using fitler" << endl;
			for(unsigned int f=0;f<nbFilter;f++){
				DescriptorsSubarray.push_back(new MatrixXd(nbDat[f],dim));
				CorresIndexSubarray.push_back(new unsigned int[nbDat[f]]);
				for(unsigned int i=0;i<nbDat[f];i++){
					CorresIndexSubarray[f][i] = FilterIndex[f*nbDatMax+i];
					for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[f]))(i,d) = _Descriptors[FilterIndex[f*nbDatMax+i]*dim+d];
				}
				SubarraySize.push_back(nbDat[f]);
			}
			subarray_defined = true;
		}else if( !subarray_filter && !subarray_label ){
			//cout << "Constructing descriptor subarrays using fitler" << endl;
			delete DescriptorsSubarray[0];
			DescriptorsSubarray.erase(DescriptorsSubarray.begin());
			SubarraySize.erase(SubarraySize.begin());
			delete[] CorresIndexSubarray[0];
			CorresIndexSubarray.erase(CorresIndexSubarray.begin());
			for(unsigned int f=0;f<nbFilter;f++){
				DescriptorsSubarray.push_back(new MatrixXd(nbDat[f],dim));
				CorresIndexSubarray.push_back(new unsigned int[nbDat[f]]);
				for(unsigned int i=0;i<nbDat[f];i++){
					CorresIndexSubarray[f][i] = FilterIndex[f*nbDatMax+i];
					for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[f]))(i,d) = _Descriptors[FilterIndex[f*nbDatMax+i]*dim+d];
				}
				SubarraySize.push_back(nbDat[f]);
			}
		}else if( ( !subarray_filter && subarray_label ) || ( !subarray_filter && !subarray_label ) ){	
			//cout << "Constructing descriptor subarrays using fitler" << endl;
			for(unsigned int i=0;i<DescriptorsSubarray.size();i++){
				delete DescriptorsSubarray[0];
				DescriptorsSubarray.erase(DescriptorsSubarray.begin());
				SubarraySize.erase(SubarraySize.begin());
				delete[] CorresIndexSubarray[0];
				CorresIndexSubarray.erase(CorresIndexSubarray.begin());
			}
			for(unsigned int f=0;f<nbFilter;f++){
				DescriptorsSubarray.push_back(new MatrixXd(nbDat[f],dim));
				CorresIndexSubarray.push_back(new unsigned int[nbDat[f]]);
				for(unsigned int i=0;i<nbDat[f];i++){
					CorresIndexSubarray[f][i] = FilterIndex[f*nbDatMax+i];
					for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[f]))(i,d) = _Descriptors[FilterIndex[f*nbDatMax+i]*dim+d];
				}
				SubarraySize.push_back(nbDat[f]);
			}
		}
		subarray_filter = true;
		subarray_label = false;
	}else if( !filter && label ){
		if( !subarray_defined ){
			//cout << "Constructing descriptor subarrays using labels" << endl;
			for(unsigned int l=0;l<nbLabels;l++){
				unsigned int nbdatlab = 0;
				for(unsigned int f=0;f<nbFilter;f++) nbdatlab += LabelsSize[f*nbLabels+l];
				DescriptorsSubarray.push_back(new MatrixXd(nbdatlab,dim));
				CorresIndexSubarray.push_back(new unsigned int[nbdatlab]);
				SubarraySize.push_back(nbdatlab);
				nbdatlab = 0;
				for(unsigned int f=0;f<nbFilter;f++){
					for(unsigned int i=0;i<LabelsSize[f*nbLabels+l];i++){
						CorresIndexSubarray[l][nbdatlab] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
						for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[l]))(nbdatlab,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
						nbdatlab++;
					}
				}
			}
			subarray_defined = true;
		}else if( !subarray_filter && !subarray_label ){
			//cout << "Constructing descriptor subarrays using labels" << endl;
			delete DescriptorsSubarray[0];
			DescriptorsSubarray.erase(DescriptorsSubarray.begin());
			SubarraySize.erase(SubarraySize.begin());
			delete[] CorresIndexSubarray[0];
			CorresIndexSubarray.erase(CorresIndexSubarray.begin());
			for(unsigned int l=0;l<nbLabels;l++){
				unsigned int nbdatlab = 0;
				for(unsigned int f=0;f<nbFilter;f++) nbdatlab += LabelsSize[f*nbLabels+l];
				DescriptorsSubarray.push_back(new MatrixXd(nbdatlab,dim));
				CorresIndexSubarray.push_back(new unsigned int[nbdatlab]);
				SubarraySize.push_back(nbdatlab);
				nbdatlab = 0;
				for(unsigned int f=0;f<nbFilter;f++){
					for(unsigned int i=0;i<LabelsSize[f*nbLabels+l];i++){
						CorresIndexSubarray[l][nbdatlab] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
						for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[l]))(nbdatlab,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
						nbdatlab++;
					}
				}
			}
		}else if( ( subarray_filter && !subarray_label ) || ( !subarray_filter && !subarray_label ) ){	
			//cout << "Constructing descriptor subarrays using labels" << endl;
			for(unsigned int i=0;i<DescriptorsSubarray.size();i++){
				delete DescriptorsSubarray[0];
				DescriptorsSubarray.erase(DescriptorsSubarray.begin());
				SubarraySize.erase(SubarraySize.begin());
				delete[] CorresIndexSubarray[0];
				CorresIndexSubarray.erase(CorresIndexSubarray.begin());
			}
			for(unsigned int l=0;l<nbLabels;l++){
				unsigned int nbdatlab = 0;
				for(unsigned int f=0;f<nbFilter;f++) nbdatlab += LabelsSize[f*nbLabels+l];
				DescriptorsSubarray.push_back(new MatrixXd(nbdatlab,dim));
				CorresIndexSubarray.push_back(new unsigned int[nbdatlab]);
				SubarraySize.push_back(nbdatlab);
				nbdatlab = 0;
				for(unsigned int f=0;f<nbFilter;f++){
					for(unsigned int i=0;i<LabelsSize[f*nbLabels+l];i++){
						CorresIndexSubarray[l][nbdatlab] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
						for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[l]))(nbdatlab,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
						nbdatlab++;
					}
				}
			}
		}
		subarray_filter = false;
		subarray_label = true;
	}else{
		if( !subarray_defined ){
			//cout << "Constructing descriptor subarrays using filter and labels" << endl;
			for(unsigned int f=0;f<nbFilter;f++){
				for(unsigned int l=0;l<nbLabels;l++){
					CorresIndexSubarray.push_back(new unsigned int[LabelsSize[f*nbLabels+l]]);
					DescriptorsSubarray.push_back(new MatrixXd(LabelsSize[f*nbLabels+l],dim));
					SubarraySize.push_back(LabelsSize[f*nbLabels+l]);
					for(unsigned int i=0;i<LabelsSize[f*nbLabels+l];i++){
						CorresIndexSubarray[f*nbLabels+l][i] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
						for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[f*nbLabels+l]))(i,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
					}
				}
			}
			subarray_defined = true;
		}else if( !subarray_filter && !subarray_label ){
			//cout << "Constructing descriptor subarrays using filter and labels" << endl;
			delete DescriptorsSubarray[0];
			DescriptorsSubarray.erase(DescriptorsSubarray.begin());
			SubarraySize.erase(SubarraySize.begin());
			delete[] CorresIndexSubarray[0];
			CorresIndexSubarray.erase(CorresIndexSubarray.begin());
			for(unsigned int f=0;f<nbFilter;f++){
				for(unsigned int l=0;l<nbLabels;l++){
					DescriptorsSubarray.push_back(new MatrixXd(LabelsSize[f*nbLabels+l],dim));
					CorresIndexSubarray.push_back(new unsigned int[LabelsSize[f*nbLabels+l]]);
					SubarraySize.push_back(LabelsSize[f*nbLabels+l]);
					for(unsigned int i=0;i<LabelsSize[f*nbLabels+l];i++){
						CorresIndexSubarray[f*nbLabels+l][i] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
						for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[f*nbLabels+l]))(i,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
					}
				}
			}
		}else if( ( subarray_filter && !subarray_label ) || ( !subarray_filter && subarray_label ) ){	
			//cout << "Constructing descriptor subarrays using filter and labels" << endl;
			for(unsigned int i=0;i<DescriptorsSubarray.size();i++){
				delete[] DescriptorsSubarray[0];
				DescriptorsSubarray.erase(DescriptorsSubarray.begin());
				SubarraySize.erase(SubarraySize.begin());
				delete[] CorresIndexSubarray[0];
				CorresIndexSubarray.erase(CorresIndexSubarray.begin());
			}
			for(unsigned int f=0;f<nbFilter;f++){
				for(unsigned int l=0;l<nbLabels;l++){
					DescriptorsSubarray.push_back(new MatrixXd(LabelsSize[f*nbLabels+l],dim));
					CorresIndexSubarray.push_back(new unsigned int[LabelsSize[f*nbLabels+l]]);
					SubarraySize.push_back(LabelsSize[f*nbLabels+l]);
					for(unsigned int i=0;i<LabelsSize[f*nbLabels+l];i++){
						CorresIndexSubarray[f*nbLabels+l][i] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
						for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[f*nbLabels+l]))(i,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
					}
				}
			}
		}
		subarray_filter = true;
		subarray_label = true;
	}
} // TODO corres index here also

void Descriptors::constructSubarrays(double &RatioTestDataset, bool filter, bool label){
	bool diffratiotest = false;
	if( _RatioTestDataset != RatioTestDataset ){
		diffratiotest = true;
		_RatioTestDataset = RatioTestDataset;
	}
	if( label && Labels.size() == 0 ){
		cerr << "The descriptors are not labelled, we cannot construct a test dataset, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	if( RatioTestDataset > 1. ){
		cerr << "The ratio between the size of testing and training dataset is higher than one, aborting" << endl;
		exit(EXIT_FAILURE);
	}

	random_device rd;
	if( filter && !label ){
		if( !subarray_defined ){
			//cout << "Constructing descriptor subarrays using fitler" << endl;
			for(unsigned int f=0;f<nbFilter;f++){
				unsigned int nbdat_test = round(nbDat[f]*RatioTestDataset);
				unsigned int nbdat_train = nbDat[f]-nbdat_test;
				DescriptorsSubarray.push_back(new MatrixXd(nbdat_train,dim));
				TestDataset.push_back(new MatrixXd(nbdat_test,dim));
				CorresIndexSubarray.push_back(new unsigned int[nbdat_train]);
				CorresIndexTestDataset.push_back(new unsigned int[nbdat_test]);
				SubarraySize.push_back(nbdat_train);
				TestDatasetSize.push_back(nbdat_test);
				vector<unsigned int> test_or_train;
				test_or_train.insert(test_or_train.end(),nbdat_test,0);
				test_or_train.insert(test_or_train.end(),nbdat_train,1);
				shuffle(test_or_train.begin(), test_or_train.end(), default_random_engine(rd()));
				unsigned int nbdatlab = 0;
				nbdat_test = 0;
				nbdat_train = 0;
				for(unsigned int i=0;i<nbDat[f];i++){
					if( test_or_train[i] == 0 ){
						for(unsigned int d=0;d<dim;d++) (*(TestDataset[f]))(nbdat_test,d) = _Descriptors[FilterIndex[f*nbDatMax+i]*dim+d];
						CorresIndexTestDataset[f][nbdat_test] = FilterIndex[f*nbDatMax+i];
						nbdat_test++;
					}else{
						for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[f]))(nbdat_train,d) = _Descriptors[FilterIndex[f*nbDatMax+i]*dim+d];
						CorresIndexSubarray[f][nbdat_train] = FilterIndex[f*nbDatMax+i];
						nbdat_train++;
					}
				}
			}
			subarray_defined = true;
		}else if( ( !subarray_filter && subarray_label ) || ( !subarray_filter && !subarray_label ) || diffratiotest ){	
			//cout << "Constructing descriptor subarrays using fitler" << endl;
			unsigned int init_size = DescriptorsSubarray.size();
			for(unsigned int i=0;i<init_size;i++){
				delete DescriptorsSubarray[0];
				DescriptorsSubarray.erase(DescriptorsSubarray.begin());
				SubarraySize.erase(SubarraySize.begin());
				delete TestDataset[0];
				TestDataset.erase(TestDataset.begin());
				TestDatasetSize.erase(TestDatasetSize.begin());
				delete[] CorresIndexSubarray[0];
				CorresIndexSubarray.erase(CorresIndexSubarray.begin());
				delete[] CorresIndexTestDataset[0];
				CorresIndexTestDataset.erase(CorresIndexTestDataset.begin());
			}
			for(unsigned int f=0;f<nbFilter;f++){
				unsigned int nbdat_test = round(nbDat[f]*RatioTestDataset);
				unsigned int nbdat_train = nbDat[f]-nbdat_test;
				DescriptorsSubarray.push_back(new MatrixXd(nbdat_train,dim));
				TestDataset.push_back(new MatrixXd(nbdat_test,dim));
				CorresIndexSubarray.push_back(new unsigned int[nbdat_train]);
				CorresIndexTestDataset.push_back(new unsigned int[nbdat_test]);
				SubarraySize.push_back(nbdat_train);
				TestDatasetSize.push_back(nbdat_test);
				vector<unsigned int> test_or_train;
				test_or_train.insert(test_or_train.end(),nbdat_test,0);
				test_or_train.insert(test_or_train.end(),nbdat_train,1);
				shuffle(test_or_train.begin(), test_or_train.end(), default_random_engine(rd()));
				unsigned int nbdatlab = 0;
				nbdat_test = 0;
				nbdat_train = 0;
				for(unsigned int i=0;i<nbDat[f];i++){
					if( test_or_train[i] == 0 ){
						for(unsigned int d=0;d<dim;d++) (*(TestDataset[f]))(nbdat_test,d) = _Descriptors[FilterIndex[f*nbDatMax+i]*dim+d];
						CorresIndexTestDataset[f][nbdat_test] = FilterIndex[f*nbDatMax+i];
						nbdat_test++;
					}else{
						for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[f]))(nbdat_train,d) = _Descriptors[FilterIndex[f*nbDatMax+i]*dim+d];
						CorresIndexSubarray[f][nbdat_train] = FilterIndex[f*nbDatMax+i];
						nbdat_train++;
					}
				}
			}
		}
		subarray_filter = true;
		subarray_label = false;
		IsTestDataset = true;
	}else if( !filter && label ){
		if( !subarray_defined ){
			//cout << "Constructing descriptor subarrays using labels" << endl;
			for(unsigned int l=0;l<nbLabels;l++){
				unsigned int nbdatlab = 0;
				for(unsigned int f=0;f<nbFilter;f++) nbdatlab += LabelsSize[f*nbLabels+l];
				unsigned int nbdat_test = round(nbdatlab*RatioTestDataset);
				unsigned int nbdat_train = nbdatlab-nbdat_test;
				DescriptorsSubarray.push_back(new MatrixXd(nbdat_train,dim));
				TestDataset.push_back(new MatrixXd(nbdat_test,dim));
				CorresIndexSubarray.push_back(new unsigned int[nbdat_train]);
				CorresIndexTestDataset.push_back(new unsigned int[nbdat_test]);
				SubarraySize.push_back(nbdat_train);
				TestDatasetSize.push_back(nbdat_test);
				vector<unsigned int> test_or_train;
				test_or_train.insert(test_or_train.end(),nbdat_test,0);
				test_or_train.insert(test_or_train.end(),nbdat_train,1);
				shuffle(test_or_train.begin(), test_or_train.end(), default_random_engine(rd()));
				nbdatlab = 0;
				nbdat_test = 0;
				nbdat_train = 0;
				for(unsigned int f=0;f<nbFilter;f++){
					for(unsigned int i=0;i<LabelsSize[f*nbLabels+l];i++){
						if( test_or_train[nbdatlab] == 0 ){
							for(unsigned int d=0;d<dim;d++) (*(TestDataset[l]))(nbdat_test,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
							CorresIndexTestDataset[l][nbdat_test] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
							nbdat_test++;
						}else{
							for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[l]))(nbdat_train,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
							CorresIndexSubarray[l][nbdat_train] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
							nbdat_train++;
						}
						nbdatlab++;
					}
				}
			}
			subarray_defined = true;
		}else if( ( subarray_filter && !subarray_label ) || diffratiotest ){	
			//cout << "Constructing descriptor subarrays using labels" << endl;
			for(unsigned int i=0;i<DescriptorsSubarray.size();i++){
				delete DescriptorsSubarray[0];
				DescriptorsSubarray.erase(DescriptorsSubarray.begin());
				SubarraySize.erase(SubarraySize.begin());
				delete TestDataset[0];
				TestDataset.erase(TestDataset.begin());
				TestDatasetSize.erase(TestDatasetSize.begin());
				delete[] CorresIndexSubarray[0];
				CorresIndexSubarray.erase(CorresIndexSubarray.begin());
				delete[] CorresIndexTestDataset[0];
				CorresIndexTestDataset.erase(CorresIndexTestDataset.begin());
			}
			for(unsigned int l=0;l<nbLabels;l++){
				unsigned int nbdatlab = 0;
				for(unsigned int f=0;f<nbFilter;f++) nbdatlab += LabelsSize[f*nbLabels+l];
				unsigned int nbdat_test = round(nbdatlab*RatioTestDataset);
				unsigned int nbdat_train = nbdatlab-nbdat_test;
				DescriptorsSubarray.push_back(new MatrixXd(nbdat_train,dim));
				TestDataset.push_back(new MatrixXd(nbdat_test,dim));
				CorresIndexSubarray.push_back(new unsigned int[nbdat_train]);
				CorresIndexTestDataset.push_back(new unsigned int[nbdat_test]);
				SubarraySize.push_back(nbdat_train);
				TestDatasetSize.push_back(nbdat_test);
				vector<unsigned int> test_or_train;
				test_or_train.insert(test_or_train.end(),nbdat_test,0);
				test_or_train.insert(test_or_train.end(),nbdat_train,1);
				shuffle(test_or_train.begin(), test_or_train.end(), default_random_engine(rd()));
				nbdatlab = 0;
				nbdat_test = 0;
				nbdat_train = 0;
				for(unsigned int f=0;f<nbFilter;f++){
					for(unsigned int i=0;i<LabelsSize[f*nbLabels+l];i++){
						if( test_or_train[nbdatlab] == 0 ){
							for(unsigned int d=0;d<dim;d++) (*(TestDataset[l]))(nbdat_test,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
							CorresIndexTestDataset[l][nbdat_test] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
							nbdat_test++;
						}else{
							for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[l]))(nbdat_train,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
							CorresIndexSubarray[l][nbdat_train] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
							nbdat_train++;
						}
						nbdatlab++;
					}
				}
			}
		}
		subarray_filter = false;
		subarray_label = true;
		IsTestDataset = true;
	}else if( filter && label ){
		if( !subarray_defined ){
			//cout << "Constructing descriptor subarrays using filter and labels" << endl;
			for(unsigned int f=0;f<nbFilter;f++){
				for(unsigned int l=0;l<nbLabels;l++){
					unsigned int nbdat_test = round(LabelsSize[f*nbLabels+l]*RatioTestDataset);
					unsigned int nbdat_train = LabelsSize[f*nbLabels+l]-nbdat_test;
					DescriptorsSubarray.push_back(new MatrixXd(nbdat_train,dim));
					TestDataset.push_back(new MatrixXd(nbdat_test,dim));
					CorresIndexSubarray.push_back(new unsigned int[nbdat_train]);
					CorresIndexTestDataset.push_back(new unsigned int[nbdat_test]);
					SubarraySize.push_back(nbdat_train);
					TestDatasetSize.push_back(nbdat_test);
					vector<unsigned int> test_or_train;
					test_or_train.insert(test_or_train.end(),nbdat_test,0);
					test_or_train.insert(test_or_train.end(),nbdat_train,1);
					shuffle(test_or_train.begin(), test_or_train.end(), default_random_engine(rd()));
					nbdat_test = 0;
					nbdat_train = 0;
					for(unsigned int i=0;i<LabelsSize[f*nbLabels+l];i++){
						if( test_or_train[i] == 0 ){
							for(unsigned int d=0;d<dim;d++) (*(TestDataset[f*nbLabels+l]))(nbdat_test,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
							CorresIndexTestDataset[f*nbLabels+l][nbdat_test] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
							nbdat_test++;
						}else{
							for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[f*nbLabels+l]))(nbdat_train,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
							CorresIndexSubarray[f*nbLabels+l][nbdat_train] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
							nbdat_train++;
						}
					}
				}
			}
			subarray_defined = true;
		}else if( ( !subarray_filter && subarray_label ) || diffratiotest ){	
			//cout << "Constructing descriptor subarrays using filter and labels" << endl;
			for(unsigned int i=0;i<DescriptorsSubarray.size();i++){
				delete[] DescriptorsSubarray[0];
				DescriptorsSubarray.erase(DescriptorsSubarray.begin());
				SubarraySize.erase(SubarraySize.begin());
				delete TestDataset[0];
				TestDataset.erase(TestDataset.begin());
				TestDatasetSize.erase(TestDatasetSize.begin());
				delete[] CorresIndexSubarray[0];
				CorresIndexSubarray.erase(CorresIndexSubarray.begin());
				delete[] CorresIndexTestDataset[0];
				CorresIndexTestDataset.erase(CorresIndexTestDataset.begin());
			}
			for(unsigned int f=0;f<nbFilter;f++){
				for(unsigned int l=0;l<nbLabels;l++){
					unsigned int nbdat_test = round(LabelsSize[f*nbLabels+l]*RatioTestDataset);
					unsigned int nbdat_train = LabelsSize[f*nbLabels+l]-nbdat_test;
					DescriptorsSubarray.push_back(new MatrixXd(nbdat_train,dim));
					TestDataset.push_back(new MatrixXd(nbdat_test,dim));
					CorresIndexSubarray.push_back(new unsigned int[nbdat_train]);
					CorresIndexTestDataset.push_back(new unsigned int[nbdat_test]);
					SubarraySize.push_back(nbdat_train);
					TestDatasetSize.push_back(nbdat_test);
					vector<unsigned int> test_or_train;
					test_or_train.insert(test_or_train.end(),nbdat_test,0);
					test_or_train.insert(test_or_train.end(),nbdat_train,1);
					shuffle(test_or_train.begin(), test_or_train.end(), default_random_engine(rd()));
					nbdat_test = 0;
					nbdat_train = 0;
					for(unsigned int i=0;i<LabelsSize[f*nbLabels+l];i++){
						if( test_or_train[i] == 0 ){
							for(unsigned int d=0;d<dim;d++) (*(TestDataset[f*nbLabels+l]))(nbdat_test,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
							CorresIndexTestDataset[f*nbLabels+l][nbdat_test] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
							nbdat_test++;
						}else{
							for(unsigned int d=0;d<dim;d++) (*(DescriptorsSubarray[f*nbLabels+l]))(nbdat_train,d) = _Descriptors[LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i]*dim+d];
							CorresIndexSubarray[f*nbLabels+l][nbdat_train] = LabelIndex[f*nbDatMaxLabel*nbLabels+l*nbDatMaxLabel+i];
							nbdat_train++;
						}
					}
				}
			}
		}
		subarray_filter = true;
		subarray_label = true;
		IsTestDataset = true;
	}
}

MatrixXd* Descriptors::getSubarray(unsigned int &subarray_nbdat, std::string filter_name, std::string label_name){
	unsigned int f = current_filter(filter_name);
	unsigned int l = current_label(label_name);
	subarray_nbdat = SubarraySize[f*nbLabels+l];
	return DescriptorsSubarray[f*nbLabels+l];
}

MatrixXd* Descriptors::getTestDataset(unsigned int &subarray_nbdat, std::string filter_name, std::string label_name){
	unsigned int f = current_filter(filter_name);
	unsigned int l = current_label(label_name);
	subarray_nbdat = TestDatasetSize[f*nbLabels+l];
	return TestDataset[f*nbLabels+l];
}

unsigned int* Descriptors::getCorresIndexSubarray(unsigned int &subarray_nbdat, std::string filter_name, std::string label_name){
	unsigned int f = current_filter(filter_name);
	unsigned int l = current_label(label_name);
	subarray_nbdat = SubarraySize[f*nbLabels+l];
	return CorresIndexSubarray[f*nbLabels+l];
}

unsigned int* Descriptors::getCorresIndexTestDataset(unsigned int &subarray_nbdat, std::string filter_name, std::string label_name){
	unsigned int f = current_filter(filter_name);
	unsigned int l = current_label(label_name);
	subarray_nbdat = TestDatasetSize[f*nbLabels+l];
	return CorresIndexTestDataset[f*nbLabels+l];
}
// END TEST subarray

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

unsigned int Descriptors::current_filter(const std::string &filter_name){
	bool already = false;
	unsigned int f = 0;
	if( filter_name != "none" ){
		for(unsigned int ft=0;ft<FilterValue.size();ft++){
			if( filter_name == FilterValue[ft] ){
				already = true;
				f = ft;
				break;
			}
		}
		if( !already ){
			cerr << "Filter value not found" << endl;
			exit(EXIT_FAILURE);
		}
	}
	return f;
}

unsigned int Descriptors::current_label(const std::string &label_name){
	bool already = false;
	unsigned int l = 0;
	if( label_name != "none" ){
		for(unsigned int lt=0;lt<Labels.size();lt++){
			if( label_name == Labels[lt] ){
				already = true;
				l = lt;
				break;
			}
		}
		if( !already ){
			cerr << "Label value not found" << endl;
			exit(EXIT_FAILURE);
		}
	}
	return l;
}


Descriptors::~Descriptors(){
	if( FilterIndex ) delete[] FilterIndex;
	if( DescriptorLabel) delete[] DescriptorLabel;
	if( DescriptorFilter ) delete[] DescriptorFilter;
	if( Labels_uint ) delete[] Labels_uint;
	if( LabelsSize ) delete[] LabelsSize;
	if( LabelIndex ) delete[] LabelIndex;
	// TEST subarray
	for(unsigned int s=0;s<DescriptorsSubarray.size();s++) if( DescriptorsSubarray[s] ) delete DescriptorsSubarray[s];
	for(unsigned int s=0;s<TestDataset.size();s++) if( TestDataset[s] ) delete TestDataset[s];
	for(unsigned int s=0;s<CorresIndexSubarray.size();s++) if( CorresIndexSubarray[s] ) delete CorresIndexSubarray[s];
	for(unsigned int s=0;s<CorresIndexTestDataset.size();s++) if( CorresIndexTestDataset[s] ) delete CorresIndexTestDataset[s];
	// END TEST subarray
	if( AreDescriptorsMine ) delete[] _Descriptors;
	delete MT;
	if( AreExtremums ){
		delete[] min_val;
		delete[] max_val;
	}
	if( IsNeighbours ){
		for(unsigned int i=0;i<Neighbours.size();i++) delete[] Neighbours[i];
	}
}
