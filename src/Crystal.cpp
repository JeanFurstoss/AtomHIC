#include "Crystal.h"
#include "MyStructs.h"
#include "AtomHicConfig.h"
#include <iostream>
#include <dirent.h>
#include <vector>
#include <fstream>
#include <sstream>
#include "MathTools.h"
#include <cmath>

using namespace std;

Crystal::Crystal(const string& crystalName){
	this->a1 = new double[3];
	this->a2 = new double[3];
	this->a3 = new double[3];
	this->a1_star = new double[3];
	this->a2_star = new double[3];
	this->a3_star = new double[3];
	this->MT = new MathTools;
	this->IsCharge = false;
	// Read the crystal database
	char *database_env = getenv("CRYSTAL_DATABASE");
	string database;
	unsigned int crystal_index;
	if (database_env) {
		database = database_env;
	} else {
		#ifdef CRYSTAL_DATABASE
		database = CRYSTAL_DATABASE;
		#endif
	}
	if( database.empty() ) cout << "Warning database environment for crystal is empty" << endl;
	else{
		DIR *dir;
		struct dirent *diread;
		vector<string> AvailableCrystals;
		const char *env = database.c_str();
		string buffer_s;
		size_t pos;
		if( (dir = opendir(env) ) != nullptr ){
			while( (diread = readdir(dir)) != nullptr ){
				buffer_s = diread->d_name;
				pos = buffer_s.find(this->database_extension);
				if(pos!=string::npos) AvailableCrystals.push_back(buffer_s.erase(buffer_s.size()-this->database_extension.size()));
			}
			closedir(dir);
		}else{
			perror("opendir");
		}

		// search into the database if crystalName is present
		bool crystalok = false;
		for(unsigned int i=0;i<AvailableCrystals.size();i++){
			if( crystalName == AvailableCrystals[i] ){
				crystal_index = i;
				this->name = crystalName; 
				crystalok = true;
				break;
			}
		}
		if( !crystalok ){
			cout << "The crystal \"" << crystalName << "\" does not exist in the crystal database, please create \"" << crystalName << ".dat\" in the /data/Crystal/ directory of AtomHic or use an other Crystal constructor" << endl;
			cout << "List of available crystals :" << endl;
			for(unsigned int i=0;i<AvailableCrystals.size();i++) cout << AvailableCrystals[i] << endl;
		}else{
			this->path2database = database+"/"+AvailableCrystals[crystal_index]+this->database_extension;
		}
		read_database();
		computeReciproqual();
	}
}

// construct a z oriented (hkl) unit cell plane
void Crystal::ConstructOrientedPlane(const int& h_p, const int& k_p, const int& l_p){
	int arr[3] = {abs(h_p),abs(k_p),abs(l_p)};
	int maxhkl = MT->max(arr,3)*5;
	double tolScalarProd = 1e-3;
	vector<double> buffer_vector_d;
	vector<int> buffer_vector_i;
	double *normalDir = new double[3];
	
	// direction normal to the wanted plane
	normalDir[0] = h_p*this->a1_star[0] + k_p*this->a2_star[0] + l_p*this->a3_star[0];	
	normalDir[1] = h_p*this->a1_star[1] + k_p*this->a2_star[1] + l_p*this->a3_star[1];	
	normalDir[2] = h_p*this->a1_star[2] + k_p*this->a2_star[2] + l_p*this->a3_star[2];

	// search the smallest orthogonal direction with integer Miller indices to the wanted one => to be aligned with the cartesian x axis
	for(int i=-maxhkl;i<maxhkl;i++){
		for(int j=-maxhkl;j<maxhkl;j++){
			for(int k=-maxhkl;k<maxhkl;k++){
				if( fabs(i*h_p+j*k_p+k*k_p) < tolScalarProd && ( i != 0 || j != 0 || k!= 0 ) ){
					buffer_vector_d.push_back(pow(i,2.)+pow(j,2.)+pow(k,2.));
					buffer_vector_i.push_back(i);
					buffer_vector_i.push_back(j);
					buffer_vector_i.push_back(k);
				}
			}
		}
	}
	unsigned int minhkl = MT->min(buffer_vector_d);
	// compute the rotation angles needed to align the wanted direction with the cartesian z axis and have the x cartesian axis corresponding to the smallest crystallographic direction possible
	double theta_z = -acos(normalDir[0]/sqrt(pow(normalDir[0],2.)+pow(normalDir[1],2.)));	
	double theta_y = -acos(normalDir[0]/sqrt(pow(normalDir[0],2.)+pow(normalDir[1],2.)));	

	delete[] normalDir;
}

void Crystal::computeReciproqual(){
	// compute cell volume
	double *mixedProd = new double[3];
	MT->mixedProd(this->a2,this->a3,mixedProd);
	this->V = MT->dotProd(this->a1,mixedProd);
	// compute reciproqual vectors
	MT->mixedProd(this->a2,this->a3,this->a1_star);
	MT->mixedProd(this->a3,this->a1,this->a2_star);
	MT->mixedProd(this->a1,this->a2,this->a3_star);
	for(unsigned int i=0;i<3;i++){
		this->a1_star[i] /= this->V;
		this->a2_star[i] /= this->V;
		this->a3_star[i] /= this->V;
	}
	delete[] mixedProd;
}

void Crystal::read_database(){
	ifstream file(this->path2database, ios::in);
	size_t pos_at, pos_x, pos_y, pos_z, pos_attype, pos_Mass, pos_At, pos_Crystal, pos_tilt;
	unsigned int line_Mass(1000), line_At(1000), buffer_uint, buffer_uint_1, count(0);
	double buffer_1, buffer_2, buffer_3, buffer_4;
	string buffer_s, buffer_s_1, buffer_s_2, line;
	if(file){
		while(file){
			getline(file,line);
			// find crystallogaphy
			pos_Crystal=line.find("CRYSTAL");
			if(pos_Crystal!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1;
				this->crystallo = buffer_s_1;
			}
			// find number of atom
			pos_at=line.find("atoms");
			if(pos_at!=string::npos){
				istringstream text(line);
				text >> buffer_uint;
				this->nbAtom = buffer_uint;
				Motif = new Atom[this->nbAtom];
			}

			// find H1 vector
			pos_x=line.find("xlo xhi");
			if(pos_x!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				this->a1[0] = buffer_2-buffer_1;
			}

			// find H2 vector
			pos_y=line.find("ylo yhi");
			if(pos_y!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				this->a2[1] = buffer_2-buffer_1;
			}

			// find H3 vector
			pos_z=line.find("zlo zhi");
			if(pos_z!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				this->a3[2] = buffer_2-buffer_1;
			}
			// search tilt component for non cubic crystal
			if( this->crystallo == "Hexagonal"){
				pos_tilt=line.find("xy xz yz");
				if(pos_tilt!=string::npos){
					istringstream text(line);
					text >> buffer_1 >> buffer_2 >> buffer_3;
					this->a2[0] = buffer_1;
				}
			}

			// find nb atom type
			pos_attype=line.find("atom types");
			if(pos_attype!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_s >> buffer_s_1;
				this->nbAtomType = buffer_1;
				AtomType = new string[this->nbAtomType];
				AtomType_uint = new unsigned int[this->nbAtomType];
				NbAtomSite = new unsigned int[this->nbAtomType];
				AtomMass = new double[this->nbAtomType];
			}

			// get lines where are the keywords Masses and Atoms to get atom type masses and positions
			pos_Mass=line.find("Masses");
			if(pos_Mass!=string::npos) line_Mass = count;
			if( count > line_Mass+1 && count < line_Mass+2+this->nbAtomType ){
				istringstream text(line);
				text >> buffer_uint >> buffer_1 >> buffer_s_1 >> buffer_s >> buffer_uint_1;
				this->AtomMass[buffer_uint-1] = buffer_1;
				this->AtomType[buffer_uint-1] = buffer_s;
				this->AtomType_uint[buffer_uint-1] = buffer_uint;
				this->NbAtomSite[buffer_uint-1] = buffer_uint_1;
			}
			pos_At=line.find("Atoms");
			if(pos_At!=string::npos){
				istringstream text(line);
				text >> buffer_s_1 >> buffer_s_2 >> buffer_s;
				if( buffer_s == "charge" ) this->IsCharge = true;
			       	line_At = count;
			}
			if( count > line_At+1 ){
				istringstream text(line);
				if( this->IsCharge ){
					text >> buffer_uint >> buffer_uint_1 >> buffer_1 >> buffer_2 >> buffer_3 >> buffer_4;
					this->Motif[buffer_uint-1].pos.x = buffer_2;
					this->Motif[buffer_uint-1].pos.y = buffer_3;
					this->Motif[buffer_uint-1].pos.z = buffer_4;
					this->Motif[buffer_uint-1].type_uint = buffer_uint_1;
					this->Motif[buffer_uint-1].charge = buffer_1;
				}else{
					text >> buffer_uint >> buffer_uint_1 >> buffer_2 >> buffer_3 >> buffer_4;
					this->Motif[buffer_uint-1].pos.x = buffer_2;
					this->Motif[buffer_uint-1].pos.y = buffer_3;
					this->Motif[buffer_uint-1].pos.z = buffer_4;
					this->Motif[buffer_uint-1].type_uint = buffer_uint_1;
				}
			}
			count += 1;
		}
		file.close();
		this->a1[1] = 0;
		this->a1[2] = 0;
		this->a2[2] = 0;
		this->a3[0] = 0;
		this->a3[1] = 0;
		if( this->crystallo == "Cubic" || this->crystallo == "Orthorhombic" ){
			this->a2[0] = 0;
		}
	}else{
		cout << "The file " << this->path2database << " cannot be openned" << endl;
	}
}

Crystal::~Crystal(){
	delete[] a1;
	delete[] a2;
	delete[] a3;
	delete[] a1_star;
	delete[] a2_star;
	delete[] a3_star;
	delete MT;
	delete[] AtomType;
	delete[] AtomType_uint;
	delete[] AtomMass;
	delete[] NbAtomSite;
	delete[] this->Motif;
}
