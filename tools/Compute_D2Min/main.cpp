//**********************************************************************************
//*   Compute_D2Min/main.cpp                                                       *
//**********************************************************************************
//* This file contains the implementation of the Compute_D2Min executable.         *
//* It allows to compute the D2Min of given dump files relative to a reference 	   *
//* dump file following Delbecq et al. 2023					   *
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
//*	- 									   *
//**********************************************************************************

#include <AtomicSystem.h>
#include <ComputeAuxiliary.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	string ReferenceFilename, ListOfFilename, buffer_string, line, StrainFilename;
	double rc;
	if( argc == 4 ){
		ReferenceFilename = argv[1];
		ListOfFilename = argv[2];
		istringstream iss_rc(argv[3]);
		iss_rc >> rc;
	}else{
		cerr << "Usage: ./Compute_D2Min ReferenceFilename Filename_of_list_of_analyzed_filename cutoff_radius" << endl;
		cerr << "The program will return the D2Min value according to the formulae in Delbecq et al. 2023" << endl;
		return EXIT_FAILURE;
	}
	Dis.Printer_NoFixedParams();
	
	AtomicSystem ReferenceSystem(ReferenceFilename);
	ComputeAuxiliary CA(&ReferenceSystem);

	unsigned int nbAt = ReferenceSystem.getNbAtom();

	vector<string> AnalyzedFilenames;	
	
	ifstream file(ListOfFilename, ios::in);
	if(file){
                while(file){
        		getline(file,line);
			if( !file ) break;
			istringstream text(line);
			text >> buffer_string;
			AnalyzedFilenames.push_back(buffer_string);
		}
	}else{
		cerr << "Could not open " << ListOfFilename << " file" << endl;
	}
	file.close();

	string outfilename;
	string suf="D2Min_";
	string ext = ".cfg";
	if( AnalyzedFilenames.size() != 0 ){
		AtomicSystem AnalyzedSystem(AnalyzedFilenames[0]);
		AnalyzedSystem.setAux(CA.Compute_D2Min(AnalyzedSystem,rc),"D2Min");
		outfilename = suf+AnalyzedFilenames[0].erase(AnalyzedFilenames[0].size()-4)+ext;
		AnalyzedSystem.printSystem_aux(outfilename,"D2Min");
		for(unsigned int i=1;i<AnalyzedFilenames.size();i++){
			AnalyzedSystem.FilenameConstructor(AnalyzedFilenames[i]);
			AnalyzedSystem.modifyAux_vec(CA.Compute_D2Min(AnalyzedSystem,rc),"D2Min");
			outfilename = suf+AnalyzedFilenames[i].erase(AnalyzedFilenames[i].size()-4)+ext;
			AnalyzedSystem.printSystem_aux(outfilename,"D2Min");
		}
	}
	Dis.ExecutionTime();	
	return 0;
}
