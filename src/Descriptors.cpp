#include "Descriptors.h"
#include "AtomHicConfig.h"
#include <iostream>
#include <dirent.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

Descriptors::Descriptors(const string &Filename){
	cout << "Reading descriptors from file : " << Filename << endl;
	ifstream file(Filename, ios::in);
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
	this->MT = new MathTools();

	ifstream file2(Filename, ios::in);
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
	cout << this->nbDat << " descriptors of dimension " << this->dim << ", successfully read from file " << Filename << endl;

}

Descriptors::~Descriptors(){
	delete[] _Descriptors;
	delete MT;
}
