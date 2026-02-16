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

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc == 0 || (argc-2)%6 != 0 ){
		cerr << "Usage: ./GenerateYAML4ACE \"Species\" rcut NRadRank1 NRadRankN LRankN MaxNBody OutputFilename" << endl;
		cerr << "This executable allows to generate a yaml file containing the informations for computing ACE descriptors." << endl;
		cerr << "It allows to specify for each combination of atomic species the parameters of the descriptors." << endl;
		cerr << "The created yaml file can then be used to compute ACE descriptors using the ComputeDescriptors executable" << endl;
		cerr << "The first 6 parameters can be specified (in order) as many times as desired. For each set of 6 parameters, a new speciesblock section will be added to the output yaml file" << endl;
		cerr << " - \"Species\" : the quote must be used for providing species, examples : \"Al\", \"Al O\", \"Al Mg O\", etc." << endl;
		cerr << " - rcut : cutoff radius " << endl; 
		cerr << " - NRadRank1 : maximum order of the 2-body radial functions used for 2-body descriptors (if the number of species in the block is higher than 2, NradRank1 will not be used)" << endl;
		cerr << " - NRadRankN : maximum order of radial functions used for N-body descriptors" << endl; 
		cerr << " - LRankN : maximum order of spherical harmonics used for N-body descriptors" << endl;
	        cerr << " - MaxNBody : maximum order of interactions" << endl;	
		cerr << "For all combinations of the three later parameters a nbody section will be generated in the yaml file containing:" << endl;
		cerr << "\t - NRadRank1 2-body lines leading to NRadRank1 radial descriptors" << endl;
		cerr << "\t -  N-body lines leading to more than descriptors (depending on possible combinations of m-spherical harmonics degree leading to rotation invariance)" << endl;
		cerr << "Total number of descriptors can be decrease by removing lines in the nbody section" << endl;
		return EXIT_FAILURE;
	}
	unsigned int nbblock = (argc-2)/6;	
	vector<double> rcut(nbblock);
	vector<unsigned int> NRadRank1(nbblock), NRadRankN(nbblock), LRankN(nbblock), MaxNBody(nbblock);
	vector<vector<string>> species(nbblock);

	string OutputFilename;

	double buffer_d;
	unsigned int buffer_uint;
	string buffer_s;
	unsigned int current_read_ind = 1;

	for(unsigned int b=0;b<nbblock;b++){
		buffer_s = argv[current_read_ind];
		istringstream spe(buffer_s);
		while( spe >> buffer_s ) species[b].push_back(buffer_s);
		current_read_ind++;

		istringstream iss_rc(argv[current_read_ind]);
		iss_rc >> buffer_d;
		rcut[b] = buffer_d;
		current_read_ind++;
		
		istringstream iss_nrad1(argv[current_read_ind]);
		iss_nrad1 >> buffer_uint;
		NRadRank1[b] = buffer_uint;
		current_read_ind++;
		
		istringstream iss_nradN(argv[current_read_ind]);
		iss_nradN >> buffer_uint;
		NRadRankN[b] = buffer_uint;
		current_read_ind++;
		
		istringstream iss_l(argv[current_read_ind]);
		iss_l >> buffer_uint;
		LRankN[b] = buffer_uint;
		current_read_ind++;
		
		istringstream iss_maxnbody(argv[current_read_ind]);
		iss_maxnbody >> buffer_uint;
		MaxNBody[b] = buffer_uint;
		current_read_ind++;
		
		cout << "Block " << b+1 << ": ";
		for(unsigned int i=0;i<species[b].size();i++) cout << species[b][i] << " ";
		cout << ", rcut = " << rcut[b] << ", NRadRank1 = " << NRadRank1[b] << ", NRadRankN = " << NRadRankN[b] << ", LRankN = " << LRankN[b] << ", interaction order = " << MaxNBody[b] << endl;
	}
	
	OutputFilename = argv[current_read_ind];
	ofstream writefile(OutputFilename);

	for(unsigned int current_write_ind = 0;current_write_ind<nbblock;current_write_ind++){
		unsigned int nbspec = species[current_write_ind].size();
		writefile << "global:" << endl;
		writefile << "  DeltaSplineBins: 0.001" << endl;
		writefile << "species:" << endl;
		writefile << "  - speciesblock: ";
		for(unsigned int s=0;s<nbspec;s++) writefile << species[current_write_ind][s] << " ";
	        writefile << endl;
		writefile << "    ndensityi: 0" << endl;
		writefile << "    npoti: " << endl;
		writefile << "    parameters: []" << endl;
		writefile << "    nradmaxi: " << NRadRankN[current_write_ind] << endl;
		writefile << "    lmaxi: " << LRankN[current_write_ind] << endl;
		writefile << "    rcutij: " << rcut[current_write_ind] << endl;
		writefile << "    dcutij: 0.01 # does not impact values of descriptors" << endl;
		writefile << "    nradbaseij: " << NRadRank1[current_write_ind] << endl;
		writefile << "    radbase: SBessel" << endl;
		writefile << "    radparameters: [1.] # does not impact values of descriptors" << endl;
		writefile << "    radcoefficients: [";
		for(unsigned int i=0;i<NRadRankN[current_write_ind];i++){
			writefile << "[";
			for(unsigned int j=0;j<LRankN[current_write_ind]+1;j++){
				writefile << "[1.";
				for(unsigned int k=1;k<NRadRank1[current_write_ind];k++) writefile << ", 1.";
				writefile << "]";
				if( j != LRankN[current_write_ind] ) writefile << ", ";
			}
			writefile << "]";
			if( i+1 != NRadRankN[current_write_ind] ) writefile << ", ";
		}
		writefile << "]" << endl;
		writefile << "    nbody:" << endl;
		// 2-body
		if( nbspec <= 2 ){
			for(unsigned int i=0;i<NRadRank1[current_write_ind];i++){
				writefile << "      - {type:";
				for(unsigned int s=0;s<nbspec;s++) writefile << " " << species[current_write_ind][s];
				if( nbspec == 1 ) writefile << " " << species[current_write_ind][0]; 
				writefile << ", nr: [" << i+1 << "], nl: [0]}" << endl;
			}
		}
		// 3-body
		if( MaxNBody[current_write_ind] >= 3 && nbspec <= 3 ){
			//unsigned int nb_nr_comb = 0;
			//for(unsigned int i=0;i<NRadRankN[current_write_ind];i++) nb_nr_comb += i+1;
			//unsigned int nb_lines = nb_nr_comb*(LRankN+1);
			for(unsigned int nr1=0;nr1<NRadRankN[current_write_ind];nr1++){
				for(unsigned int nr2=0;nr2<NRadRankN[current_write_ind];nr2++){
					for(unsigned int l=0;l<LRankN[current_write_ind]+1;l++){
						writefile << "      - {type:";
						for(unsigned int s=0;s<3;s++) writefile << " " << species[current_write_ind][0];
						writefile << ", nr: [" << nr1+1 <<", " << nr2+1 << "], nl: [" << l << ", " << l << "]}" << endl;
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
