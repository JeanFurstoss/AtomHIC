// AtomHic library files
#include "AtomHicConfig.h"
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"

using namespace std;

int main(int argc, char *argv[])
{
	if( argc < 4 ){
		cerr << "Usage: SaveGaussianMixture_Steinhardt Filename_L_Q number_of_dim CrystalType DefectName AtomType" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename, CrystalType, DefectName, line, AtomType;	
	unsigned int dim, count;
	int buffer_int;
	double buffer_d;
	InputFilename = argv[1];
	istringstream iss_l(argv[2]);
	iss_l >> dim;
	CrystalType = argv[3];
	DefectName = argv[4];
	AtomType = argv[5];
	vector<vector<double>> data;
	ifstream file(InputFilename, ios::in);
	count = 0;
	if(file){
                do{
                        getline(file,line);
			if( !file ) break;
			if( count == 0 ) data.push_back(vector<double>());
			istringstream text(line);
			text >> buffer_int >> buffer_d;
			data[data.size()-1].push_back(buffer_d);
			count += 1;
			if( count == dim ) count = 0;
		}while(file);
	}
	MathTools MT;
	vector<double> mu;
	vector<vector<double>> C, C_inv;
	long double det;
	MT.MultidimGaussian(data,mu,C);
	MT.invMat_LU(C,C_inv,det);
	cout << "deter : " << det << endl;
	string database;	
	string GMM="/GaussianMixtureModel/"; // TODO create dir if does not exist
	string ext=".dat";
	if( CrystalType == "Forsterite" ){
		#ifdef STEINHARDT_FORSTERITE_DATABASE
		database = STEINHARDT_FORSTERITE_DATABASE;
		#endif
	}else if( CrystalType == "Periclase" ){
		#ifdef STEINHARDT_PERICLASE_DATABASE
		database = STEINHARDT_PERICLASE_DATABASE;
		#endif
	}
	if( database.empty() ){
		cerr << "Warning database environment for crystal is empty" << endl;
		exit(EXIT_FAILURE);
	}
	string fullpathname=database.c_str()+GMM+DefectName+ext;
	ifstream ifile;
	ifile.open(fullpathname);
	string buffer_s, buffer_s_1;
        vector<string> attype_already;
	size_t pos_attype, pos_dim;
	if( ifile ){
		while(ifile){
			getline(ifile,line);
			pos_attype = line.find("ATOM_TYPE");
			if(pos_attype!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1;
				attype_already.push_back(buffer_s_1);
			}
			pos_dim = line.find("NUMBER_OF_DIMENSION");
			if(pos_dim!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_d;
			}

		}
		ifile.close();
		bool already = false;
		for(unsigned int i=0;i<attype_already.size();i++){
			if( AtomType == attype_already[i] ){
				already = true;
				break;
			}
		}
		if( already ){
			cerr << "This atom type has already been stored for this defect, aborting.." << endl;
			exit(EXIT_FAILURE);
		}else if( buffer_d != dim ){
			cerr << "The number of dimension used for the atom type already stored is different to the provided one, aborting.." << endl;
			exit(EXIT_FAILURE);
		}else{
			ofstream writefile(fullpathname, std::ios::app);
			writefile << "ATOM_TYPE " << AtomType << endl;
			writefile << "DETERMINANT " << det << endl;
			writefile << "ESPERANCE";
			for(unsigned int i=0;i<dim;i++) writefile << " " << mu[i];
			writefile << endl;
			writefile << "INVERSE OF COVARIANCE MATRIX" << endl;
			for(unsigned int i=0;i<dim;i++){
				for(unsigned int j=0;j<dim;j++) writefile << C_inv[i][j] << " ";
				writefile << endl;
			}
			writefile.close();
		}
	}else{
		ofstream writefile(fullpathname);
		writefile << "NUMBER_OF_DIMENSION " << dim << endl;
		writefile << "ATOM_TYPE " << AtomType << endl;
		writefile << "DETERMINANT " << det << endl;
		writefile << "ESPERANCE";
		for(unsigned int i=0;i<dim;i++) writefile << " " << mu[i];
		writefile << endl;
		writefile << "INVERSE OF COVARIANCE MATRIX" << endl;
		for(unsigned int i=0;i<dim;i++){
			for(unsigned int j=0;j<dim;j++) writefile << C_inv[i][j] << " ";
			writefile << endl;
		}
		writefile.close();
	}

	return 0;
}
