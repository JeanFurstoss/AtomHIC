//**********************************************************************************
//*   GenerateYAML4ACE/main.cpp                                                    *
//**********************************************************************************
//* This file contains the implementation of the GenerateYAML4ACE executable.      *
//* It allows to generate a yaml file for computing ACE descriptors ...  *
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
//* What is still needed to do here:                                               *
//*	- modify the code when other descriptors are implemented         	   *
//**********************************************************************************

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <Displays.h>

using namespace std;

void exec_msg(){
	cerr << "Usage: ./GenerateYAML4ACE NumberOfBlock \"Species\" rcut K NRad1 NRad.. NRadK L2 L.. LK (Lint) OutputFilename" << endl;
	cerr << "This executable allows to generate a yaml file containing the informations for computing ACE descriptors." << endl;
	cerr << "It allows to specify for each combination of atomic species the parameters of the descriptors." << endl;
	cerr << "The created yaml file can then be used to compute ACE descriptors using the ComputeDescriptors executable" << endl;
	cerr << "The \"Species rcut K NRad1 NRad.. NRadK L1 L.. LK (Lint)\" parameters constitute a speciesblock for a given combination of specie(s). The number of provided speciesblock should be given as a first argument of the executable (NumberOfBlock). For each set of parameters, a new speciesblock section will be added to the output yaml file" << endl;
	cerr << " - \"Species\": the quote must be used for providing species, examples : \"Al\", \"Al O\", \"Al Mg O\", etc." << endl;
	cerr << " - K: interaction rank for this block (1 => 2-body, 2 => 2-body + 3-body, 3 => 2-body + 3-body + 4-body, etc. with maximum of 6-body interaction)" << endl;	
	cerr << " - rcut: cutoff radius " << endl; 
	cerr << " - NRad1 NRad.. NRadK: maximum order of radial functions used for N-body descriptors (K values must be provided)" << endl;
	cerr << " - L2 L.. LK : maximum order of spherical harmonics used for N-body descriptors (K-1 values must be provided as L is always 0 for 2-body interactions)" << endl;
	cerr << " - Lint : maximum degree of spherical harmonics for 5-body descriptors (it means that Lint should only be provided if K >= 4)" << endl;
	cerr << "Total number of descriptors can be decrease by removing lines in the nbody section" << endl;
}

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc == 0 || argc < 9 ){
		exec_msg();
		return EXIT_FAILURE;
	}

	unsigned int nbblock;

	string OutputFilename;

	unsigned int current_read_ind = 1;
	istringstream iss_nbblock(argv[current_read_ind]);
	iss_nbblock >> nbblock;
	current_read_ind++;
	
	vector<vector<string>> species(nbblock);
	vector<double> rcut(nbblock);
        vector<unsigned int> Rank(nbblock);
	vector<vector<unsigned int>> NRad(nbblock), L(nbblock), Lint(nbblock);
	unsigned int maxNRad1(0), maxNRadK(0), maxL(0);
	
	double buffer_d;
	unsigned int buffer_uint;
	string buffer_s;
	
	for(unsigned int b=0;b<nbblock;b++){
		buffer_s = argv[current_read_ind];
		istringstream spe(buffer_s);
		while( spe >> buffer_s ) species[b].push_back(buffer_s);
		current_read_ind++;

		istringstream iss_rc(argv[current_read_ind]);
		iss_rc >> buffer_d;
		rcut[b] = buffer_d;
		current_read_ind++;
		
		istringstream iss_rank(argv[current_read_ind]);
		iss_rank >> buffer_uint;
		Rank[b] = buffer_uint;
		current_read_ind++;

		istringstream iss_nrad1(argv[current_read_ind]);
		iss_nrad1 >> buffer_uint;
		NRad[b].push_back(buffer_uint);
		current_read_ind++;
		if( buffer_uint > maxNRad1 ) maxNRad1 = buffer_uint;
		for(unsigned int n=1;n<Rank[b];n++){
			istringstream iss_nrad(argv[current_read_ind]);
			iss_nrad >> buffer_uint;
			NRad[b].push_back(buffer_uint);
			current_read_ind++;
			if( buffer_uint > maxNRadK ) maxNRadK = buffer_uint;
		}
	
		L[b].push_back(0);	
		for(unsigned int n=1;n<Rank[b];n++){
			istringstream iss_l(argv[current_read_ind]);
			iss_l >> buffer_uint;
			L[b].push_back(buffer_uint);
			current_read_ind++;
			if( buffer_uint > maxL ) maxL = buffer_uint;
		}

		if( Rank[b] >= 4 ){
			istringstream iss_lint(argv[current_read_ind]);
			iss_lint >> buffer_uint;
			Lint[b].push_back(buffer_uint);
			current_read_ind++;
		}
		
		cout << "Block " << b+1 << ": ";
		for(unsigned int i=0;i<species[b].size();i++) cout << species[b][i] << " ";
		cout << ", rcut = " << rcut[b] << ", rank = " << Rank[b] << ", NRad = [" << NRad[b][0];
		for(unsigned int n=1;n<Rank[b];n++) cout << ", " << NRad[b][n];
		cout << "], l = [0";
		for(unsigned int n=1;n<Rank[b];n++) cout << ", " << L[b][n];
		cout << "]";
		if( Rank[b] >= 4 ) cout << ", lint = [" << Lint[b][0] << "]";
	        cout << endl;
	}
	
	OutputFilename = argv[current_read_ind];
	ofstream writefile(OutputFilename);

	writefile << "global:" << endl;
	writefile << "  DeltaSplineBins: 0.001" << endl;
	writefile << "species:" << endl;
	for(unsigned int b=0;b<nbblock;b++){
		unsigned int nbspec = species[b].size();
		writefile << "  - speciesblock: ";
		for(unsigned int s=0;s<nbspec;s++) writefile << species[b][s] << " ";
	        writefile << endl;
		writefile << "    ndensityi: 0" << endl;
		writefile << "    npoti: " << endl;
		writefile << "    parameters: []" << endl;
		writefile << "    nradmaxi: " << maxNRadK << endl;
		writefile << "    lmaxi: " << maxL << endl;
		writefile << "    rcutij: " << rcut[b] << endl;
		writefile << "    dcutij: 0.01 # does not impact values of descriptors" << endl;
		writefile << "    nradbaseij: " << maxNRad1 << endl;
		writefile << "    radbase: SBessel" << endl;
		writefile << "    radparameters: [1.] # does not impact values of descriptors" << endl;
		writefile << "    radcoefficients: [";
		for(unsigned int i=0;i<maxNRadK;i++){
			writefile << "[";
			for(unsigned int j=0;j<maxL+1;j++){
				writefile << "[1.";
				for(unsigned int k=1;k<maxNRad1;k++) writefile << ", 1.";
				writefile << "]";
				if( j != maxL ) writefile << ", ";
			}
			writefile << "]";
			if( i+1 != maxNRadK ) writefile << ", ";
		}
		writefile << "]" << endl;
		writefile << "    nbody:" << endl;
		// 2-body
		if( nbspec <= 2 && Rank[b] >= 1 ){
			for(unsigned int i=0;i<NRad[b][0];i++){
				writefile << "      - {type:";
				for(unsigned int s=0;s<nbspec;s++) writefile << " " << species[b][s];
				if( nbspec == 1 ) writefile << " " << species[b][0]; 
				writefile << ", nr: [" << i+1 << "], nl: [0]}" << endl;
			}
		}
		// 3-body
		if( Rank[b] >= 2 && nbspec <= 3 ){
			vector<string> spec_comb;
			if( nbspec == 1 ) spec_comb.push_back(species[b][0]+" "+species[b][0]+" "+species[b][0]);
			else if( nbspec == 3 ) spec_comb.push_back(species[b][0]+" "+species[b][1]+" "+species[b][2]);
			else{
				spec_comb.push_back(species[b][0]+" "+species[b][0]+" "+species[b][1]);
				spec_comb.push_back(species[b][0]+" "+species[b][1]+" "+species[b][1]);
			}
			for(unsigned int s=0;s<spec_comb.size();s++){
				for(unsigned int nr1=0;nr1<NRad[b][1];nr1++){
					for(unsigned int nr2=0;nr2<NRad[b][1];nr2++){
						for(unsigned int l=0;l<L[b][1]+1;l++){
							writefile << "      - {type: " << spec_comb[s];
							writefile << ", nr: [" << nr1+1 <<", " << nr2+1 << "], nl: [" << l << ", " << l << "]}" << endl;
						}
					}
				}
			}
		}
		// 4-body
		if( Rank[b] >= 3 && nbspec <= 4 ){
			vector<string> spec_comb;
			if( nbspec == 1 ) spec_comb.push_back(species[b][0]+" "+species[b][0]+" "+species[b][0]+" "+species[b][0]);
			else if( nbspec == 2 ){
				spec_comb.push_back(species[b][0]+" "+species[b][1]+" "+species[b][1]+" "+species[b][1]);
				spec_comb.push_back(species[b][0]+" "+species[b][0]+" "+species[b][1]+" "+species[b][1]);
				spec_comb.push_back(species[b][0]+" "+species[b][0]+" "+species[b][0]+" "+species[b][1]);
			}else if( nbspec == 3 ){
				spec_comb.push_back(species[b][0]+" "+species[b][1]+" "+species[b][2]+" "+species[b][2]);
				spec_comb.push_back(species[b][0]+" "+species[b][1]+" "+species[b][1]+" "+species[b][2]);
				spec_comb.push_back(species[b][0]+" "+species[b][0]+" "+species[b][1]+" "+species[b][2]);
			}else spec_comb.push_back(species[b][0]+" "+species[b][1]+" "+species[b][2]+" "+species[b][3]);
			for(unsigned int s=0;s<spec_comb.size();s++){
				for(unsigned int nr1=0;nr1<NRad[b][2];nr1++){
					for(unsigned int nr2=0;nr2<NRad[b][2];nr2++){
						for(unsigned int nr3=0;nr3<NRad[b][2];nr3++){
							for(unsigned int l1=0;l1<L[b][2]+1;l1++){
								for(unsigned int l2=0;l2<L[b][2]+1;l2++){
									for(unsigned int l3=0;l3<L[b][2]+1;l3++){
										if( (l1+l2+l3)%2 == 0 ){
											writefile << "      - {type: " << spec_comb[s];
											writefile << ", nr: [" << nr1+1 <<", " << nr2+1 << ", " << nr3+1 << "], nl: [" << l1 << ", " << l2 << ", " << l3 << "], lint: [" << l3 << "]}" << endl;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		// 5-body
		if( Rank[b] >= 4 && nbspec <= 5 ){
			vector<string> spec_comb;
			spec_comb.push_back(species[b][0]+" "+species[b][0]+" "+species[b][0]+" "+species[b][0]+" "+species[b][0]);
			for(unsigned int s=0;s<spec_comb.size();s++){
				for(unsigned int nr1=0;nr1<NRad[b][3];nr1++){
					for(unsigned int nr2=0;nr2<NRad[b][3];nr2++){
						for(unsigned int nr3=0;nr3<NRad[b][3];nr3++){
							for(unsigned int nr4=0;nr4<NRad[b][3];nr4++){
								for(unsigned int l1=0;l1<L[b][3]+1;l1++){
									for(unsigned int l2=0;l2<L[b][3]+1;l2++){
										for(unsigned int l3=0;l3<L[b][3]+1;l3++){
											for(unsigned int l4=0;l4<L[b][3]+1;l4++){
												if( (l1+l2+l3+l4)%2 == 0 ){
													for(unsigned int lint=0;lint<Lint[b][0]+1;lint++){
														writefile << "      - {type: " << spec_comb[s];
														writefile << ", nr: [" << nr1+1 <<", " << nr2+1 << ", " << nr3+1 << ", " << nr4+1 << "], nl: [" << l1 << ", " << l2 << ", " << l3 << ", " << l4 << "], lint: [" << lint << ", " << lint << "]}" << endl;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	Dis.ExecutionTime();
	// parameters having default values :
	// - r_in: 0	
	return 0;
}
