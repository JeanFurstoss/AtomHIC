// AtomHic library files
#include "AtomicSystem.h"
#include "SteinhardtDescriptors.h"
#include "Crystal.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

void print2database_knownsites(string path, unsigned int *AtomSite, vector<vector<double>> UniqueBondOri, int l_sph, double rc){
	vector<string> file_content;
	ifstream file(path);
	if( file ){
		string line;
		while(getline(file,line)) file_content.push_back(line);
	}else{
		cerr << "Database file cannot be opened, aborting" << endl;
		exit(EXIT_FAILURE);
	}

	path += "_tmp";
	ofstream writefile(path);
	for(unsigned int f=0;f<file_content.size();f++) writefile << file_content[f] << endl;
	writefile << endl;
	writefile << "REFERENCE_BOND_ORIENTATIONAL_PARAMETERS" << endl;
	writefile << "STEINHARDT_MODE OneL" << endl;
	writefile << "NUMBER_OF_DIMENSION "+to_string(l_sph) << endl;
	writefile << "CUTOFF_RADIUS "+to_string(rc) << endl;
	for(unsigned int t=0;t<UniqueBondOri.size();t++){
		
		for(unsigned int s=0;s<UniqueBondOri[t].size();s++) writefile << t+1 << " " << s+1 << " " << UniqueBondOri[t][s] << endl;
	}
	writefile.close();
	cout << "We have printed the new database file in " << path << endl;
	cout << "Verify its content and it is ok you can replace the old file by the new one" << endl;
}

void print2database_unknownsites(string path, unsigned int *AtomSite, vector<vector<double>> UniqueBondOri, int l_sph, double rc){
	vector<string> file_content;
	ifstream file(path);
	if( file ){
		string line;
		unsigned int count(0), line_masses(1000), line_at(1000), nbAt(0), nbAtType(0);
		size_t pos_masses, pos_nbat, pos_at, pos_nbattype;
		string buffer_s;
		while(getline(file,line)){
			bool StoreLine = true;
			pos_nbat=line.find("atoms");
			if(pos_nbat!=string::npos){
				istringstream text(line);
				text >> nbAt;
			}
			pos_nbattype=line.find("atom types");
			if(pos_nbattype!=string::npos){
				istringstream text(line);
				text >> nbAtType;
			}
			pos_masses=line.find("Masses");
			if(pos_masses!=string::npos) line_masses = count;
			if( count > line_masses+1 && count < line_masses+nbAtType+2 ){
				StoreLine = false;
				istringstream text(line);
				unsigned int type;
				text >> type;
				string to_print = "\t"+to_string(type)+"\t";
				for(unsigned int i=0;i<3;i++){
					text >> buffer_s;
					to_print += buffer_s+"\t";
				}
				to_print += to_string(UniqueBondOri[type-1].size());
				file_content.push_back(to_print);
			}
			pos_at=line.find("Atoms #");
			if(pos_at!=string::npos) line_at = count;
			if( count > line_at+1 && count < line_at+nbAt+2 ){
				StoreLine = false;
				istringstream text(line);
				unsigned int index;
				text >> index;
				string to_print = "\t"+to_string(index)+"\t";
				for(unsigned int i=0;i<5;i++){
					text >> buffer_s;
					to_print += buffer_s+"\t";
				}
				to_print += to_string(1+AtomSite[index-1]);
				file_content.push_back(to_print);
			}
			if( StoreLine ) file_content.push_back(line);
			count++;
		}
		file.close();
	}else{
		cerr << "Database file cannot be opened, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	
	path += "_tmp";
	ofstream writefile(path);
	
	for(unsigned int f=0;f<file_content.size();f++) writefile << file_content[f] << endl;
	writefile << endl;
	writefile << "REFERENCE_BOND_ORIENTATIONAL_PARAMETERS" << endl;
	writefile << "STEINHARDT_MODE OneL" << endl;
	writefile << "NUMBER_OF_DIMENSION "+to_string(l_sph) << endl;
	writefile << "CUTOFF_RADIUS "+to_string(rc) << endl;
	for(unsigned int t=0;t<UniqueBondOri.size();t++){
		
		for(unsigned int s=0;s<UniqueBondOri[t].size();s++) writefile << t+1 << " " << s+1 << " " << UniqueBondOri[t][s] << endl;
	}
	writefile.close();
	cout << "We have printed the new database file in " << path << endl;
	cout << "Verify its content and it is ok you can replace the old file by the new one" << endl;
}


int main(int argc, char *argv[])
{
	if( argc < 2 ){
		cerr << "Usage: ./SaveNonCSCrystalBondOriParam CrystalName (lo_L hi_L lo_rc hi_rc min_diff min_val)" << endl;
		cerr << "This executable will store reference bond orientational parameters for the different site of a non centrosymmetric crystal in the Crystal database of AtomHIC" << endl;
		cerr << "The CrystalName argument should be the name of a .ath file present in /data/Crystal/ containing the informations of the crystal" << endl;
		cerr << "An example of such file can be found in /data/ExampleFiles/Crystal.ath" << endl;
		cerr << "If only CrystalName is provided, we will use the STEINHARDT_DESCRIPTORS_RCUT and STEINHARDT_DESCRIPTORS_L_SPH parameters in /data/FixedParameters/FixedParameters.dat (in this case you also have to put STEINHARDT_DESCRIPTORS_MODE to OneL" << endl;
		cerr << "In the other case, the program will search from lo_L to hi_L (step=1) spherical harmonics degree and lo_rc to hi_rc (step=0.5) cutoff radius, the lowest values allowing to separate the different crystal site by min_diff and with all values higher than min_val" << endl;
		cerr << "Typical arguments could be lo_L = (number of atom type)*2, lo_rc = (smallest_cell_vector), min_diff = 0.05, min_val = 0.5" << endl;
		cerr << "If the automatic founding of parameters does not success you can use the output of the executable to choose your prefered values, put it in /data/FixedParameters/FixedParameters.dat and then use this executable with CrystalName only to force these values" << endl;
		return EXIT_FAILURE;
	}
	
	cout << "Calculation running using " << omp_get_max_threads() << " threads" << endl;
	
	auto start = high_resolution_clock::now();
	
	MathTools MT;

	string CrystalType = argv[1];
	bool AutoFitOptiParam;
	unsigned int lo_L, hi_L;
	double lo_rc, hi_rc, min_diff, min_val;
	if( argc != 2 ){	
		istringstream iss_lo_L(argv[2]);
		iss_lo_L >> lo_L;
		istringstream iss_hi_L(argv[3]);
		iss_hi_L >> hi_L;
		istringstream iss_lo_rc(argv[4]);
		iss_lo_rc >> lo_rc;
		istringstream iss_hi_rc(argv[5]);
		iss_hi_rc >> hi_rc;
		istringstream iss_min_diff(argv[6]);
		iss_min_diff >> min_diff;
		istringstream iss_min_val(argv[7]);
		iss_min_val >> min_val;
		AutoFitOptiParam = true;
	}else AutoFitOptiParam = false;

	double zero_num = 1e-5;
	
	Crystal MyCrystal(CrystalType);
	if( MyCrystal.getIsReferenceBondOriParam() ){
		cerr << "The reference bond orientional parameters have already been computed for this crystal, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	string path = MyCrystal.getDatabasePath(CrystalType);
	bool AreSitesKnown = true;
	if( MyCrystal.getNbAtomSite(1) == 0 ) AreSitesKnown = false;
	
	AtomicSystem MySystem(MyCrystal.getMotif(),MyCrystal.getNbAtom(),&MyCrystal,MyCrystal.getA1(),MyCrystal.getA2(),MyCrystal.getA3());
		
	unsigned int *AtomSite = new unsigned int[MyCrystal.getNbAtom()];
	vector<string> Properties;

	if( !AutoFitOptiParam ){
		Properties.clear();
		Properties.push_back("STEINHARDT_MODE OneL");
		SteinhardtDescriptors MyDescriptor(&MySystem,Properties);
		vector<vector<double>> UniqueBondOri;
		for(unsigned int t=0;t<MyCrystal.getNbAtomType();t++) UniqueBondOri.push_back(vector<double>());
		for(unsigned int n=0;n<MyCrystal.getNbAtom();n++){
			bool already = false;
			unsigned int current_t = 0;
			for(unsigned int t=0;t<UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1].size();t++){
				current_t = t;
				if( fabs(MyDescriptor.getDescriptors()[n]-UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1][t]) < zero_num ){
					already = true;
					break;
				}
			}
			if( !already ){
				if( UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1].size() != 0 ) current_t++;
				UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1].push_back(MyDescriptor.getDescriptors()[n]);
			}
			AtomSite[n] = current_t;
		}
		if( AreSitesKnown ){
			bool ok = true;
			for(unsigned int t=0;t<MyCrystal.getNbAtomType();t++){
		       		if( MyCrystal.getNbAtomSite(t+1) != UniqueBondOri[t].size() ){
					ok = false;
					break;
				}
			}
			if( !ok ){
				cerr << "The number of crystallographic site does not correspond to the number given in the database file, aborting" << endl;
				cerr << "If you don't know the crystallographic sites of the crystal, please indicate it in the database file" << endl;
				exit(EXIT_FAILURE);
			}
			bool ok_site = true;
			for(unsigned int n=0;n<MyCrystal.getNbAtom();n++){
				if( AtomSite[n] != MyCrystal.getAtomSite(n) ){
					ok_site = false;
					break;
				}
			}
			if( !ok_site ){
				cerr << "The crystallographic sites do not correspond to the provided ones in the database file, aborting" << endl;
				cerr << "If you don't know the crystallographic sites of the crystal, please indicate it in the database file" << endl;
				exit(EXIT_FAILURE);
			}

			print2database_knownsites(path,AtomSite,UniqueBondOri,MyDescriptor.get_l_sph(),MyDescriptor.get_rc());
		}else print2database_unknownsites(path,AtomSite,UniqueBondOri,MyDescriptor.get_l_sph(),MyDescriptor.get_rc());

	}else{
		bool ok_all;
		unsigned int current_l;
		double current_rc;
		unsigned int nb_l_steps = hi_L-lo_L+1;
		unsigned int nb_r_steps = 2*(hi_rc-lo_rc)+1;
		for(unsigned int l_inc=0;l_inc<nb_l_steps;l_inc++){
			for(unsigned int rc_inc=0;rc_inc<nb_r_steps;rc_inc++){
				current_l = l_inc+lo_L;
				current_rc = ((double) rc_inc*.5)+lo_rc;
				cout << "L = " << current_l << ", RCUT = " << current_rc << endl;
				Properties.clear();
				Properties.push_back("STEINHARDT_MODE OneL");
				Properties.push_back("NUMBER_OF_DIMENSION "+to_string(current_l));
				Properties.push_back("CUTOFF_RADIUS "+to_string(current_rc));
				
				SteinhardtDescriptors MyDescriptor(&MySystem,Properties);
				
				vector<vector<double>> UniqueBondOri;
				for(unsigned int t=0;t<MyCrystal.getNbAtomType();t++) UniqueBondOri.push_back(vector<double>());
				for(unsigned int n=0;n<MyCrystal.getNbAtom();n++){
					bool already = false;
					unsigned int current_t = 0;
					for(unsigned int t=0;t<UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1].size();t++){
						current_t = t;
						if( fabs(MyDescriptor.getDescriptors()[n]-UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1][t]) < zero_num ){
							already = true;
							break;
						}
					}
					if( !already ){
						if( UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1].size() != 0 ) current_t++;
						UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1].push_back(MyDescriptor.getDescriptors()[n]);
					}
					AtomSite[n] = current_t;
				}
				vector<bool> opti_ok;
				for(unsigned int t=0;t<MyCrystal.getNbAtomType();t++){
				       	if( AreSitesKnown && UniqueBondOri[t].size() != MyCrystal.getNbAtomSite(t+1) ){
						opti_ok.push_back(false);
						break;
					}
					opti_ok.push_back(true);
					for(unsigned int i1=0;i1<UniqueBondOri[t].size();i1++){
					        if( UniqueBondOri[t][i1] < min_val ){
							opti_ok[t] = false;
							break;
						}	
						for(unsigned int i2=i1+1;i2<UniqueBondOri[t].size();i2++){
							if( fabs(UniqueBondOri[t][i1]-UniqueBondOri[t][i2]) < min_diff ){
								opti_ok[t] = false;
								break;
							}
						}
						if( !opti_ok[t] ) break;
					}
					if( !opti_ok[t] ) break;
				}
				ok_all = true;
				for(unsigned int t=0;t<opti_ok.size();t++) ok_all *= opti_ok[t];
				for(unsigned int i=0;i<MyCrystal.getNbAtomType();i++){
					cout << "Number of site for " << MyCrystal.getAtomType()[i] << " ions, " << UniqueBondOri[i].size() << ", with values : " << endl;
					for(unsigned int t=0;t<UniqueBondOri[i].size();t++) cout << UniqueBondOri[i][t] << " ";
					cout << endl;
				}
				if( ok_all ){
					if( AreSitesKnown ){
						bool ok_site = true;
						for(unsigned int n=0;n<MyCrystal.getNbAtom();n++){
							if( AtomSite[n] != MyCrystal.getAtomSite(n) ){
								ok_site = false;
								break;
							}
						}
						if( ok_site ) break;
						else{
							cout << "Found atom sites do not correspond to the one given in the database file" << endl;
							ok_all = false;
						}
					}else break;
				}
				//MySystem.setAux_vec(MyDescriptor.getDescriptors(),1,"BondOri_"+to_string(rc_inc+lo_rcut)+"_"+to_string(l_inc+lo_L));
				//MySystem.printSystem_aux("Out_"+to_string(rc_inc+lo_rcut)+"_"+to_string(l_inc+lo_L)+".cfg","BondOri_"+to_string(rc_inc+lo_rcut)+"_"+to_string(l_inc+lo_L));
			}
			if( ok_all ) break;
		}
		if( ok_all ){
			cout << "Optimal parameters found" << endl;
			cout << "l_sph = " << current_l << endl;
			cout << "rcut = " << current_rc << endl;
			Properties.clear();
			Properties.push_back("STEINHARDT_MODE OneL");
			Properties.push_back("NUMBER_OF_DIMENSION "+to_string(current_l));
			Properties.push_back("CUTOFF_RADIUS "+to_string(current_rc));
			
			SteinhardtDescriptors MyDescriptor(&MySystem, Properties);
			vector<vector<double>> UniqueBondOri;
			for(unsigned int t=0;t<MyCrystal.getNbAtomType();t++) UniqueBondOri.push_back(vector<double>());
			for(unsigned int n=0;n<MyCrystal.getNbAtom();n++){
				bool already = false;
				unsigned int current_t = 0;
				for(unsigned int t=0;t<UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1].size();t++){
					current_t = t;
					if( fabs(MyDescriptor.getDescriptors()[n]-UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1][t]) < zero_num ){
						already = true;
						break;
					}
				}
				if( !already ){
					if( UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1].size() != 0 ) current_t++;
					UniqueBondOri[MyCrystal.getMotif()[n].type_uint-1].push_back(MyDescriptor.getDescriptors()[n]);
				}
				AtomSite[n] = current_t;
			}
			if( AreSitesKnown ) print2database_knownsites(path,AtomSite,UniqueBondOri,MyDescriptor.get_l_sph(),MyDescriptor.get_rc());
			else print2database_unknownsites(path,AtomSite,UniqueBondOri,MyDescriptor.get_l_sph(),MyDescriptor.get_rc());
		}else{
			cout << "Optimal parameters were not found because either:" << endl;
			cout << "1. The difference between bond orientational parameters for the different sites are lower than min_diff" << endl;
			cout << "2. The computed values are higher than min_val" << endl;
			cout << "3. The parameters do not allow to retrieve the index of crystallographic sites provided in the crystal database" << endl;
			cout << "You can either, modify the values max_l_inc and max_rc_inc to increase the parameter space for the reasearch" << endl;
			cout << "Or search in the output which values are ok for you and force these values by providing only the CrystalName to this executable" << endl;
			cout << "For case 3. you can indicate that the crystallographic sites are unknown in the crystal database file" << endl;
		}
	}
	
	delete[] AtomSite;

	auto end = high_resolution_clock::now();	
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Resolution time = " << duration.count() << endl;
	
	return 0;
}
