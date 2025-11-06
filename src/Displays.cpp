//**********************************************************************************
//*   Displays.cpp                                                                   *
//**********************************************************************************
//* This file contains the implementation of the Displays class                      *
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


#include "Displays.h"
#include "AtomHicConfig.h"
#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include <iomanip>

using namespace std; 
using namespace std::chrono;

Displays::Displays(){
	time_beg = high_resolution_clock::now();	
}

void Displays::Logo(){

    cout << R"(
       █████╗ ████████╗ ██████╗ ███╗   ███╗██╗  ██╗██╗ ██████╗
      ██╔══██╗╚══██╔══╝██╔═══██╗████╗ ████║██║  ██║██║██╔════╝
      ███████║   ██║   ██║   ██║██╔████╔██║███████║██║██║
      ██╔══██║   ██║   ██║   ██║██║╚██╔╝██║██╔══██║██║██║
      ██║  ██║   ██║   ╚██████╔╝██║ ╚═╝ ██║██║  ██║██║╚██████╗
      ╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝ ╚═════╝
    )";
    cout << "                  (C) 2025 Jean Furstoss " << endl;

    cout << "              =====================================" << endl;
    cout << endl;
    cout << "Calculation running using " << omp_get_max_threads() << " threads" << endl;
    cout << endl;
}

void Displays::ExecutionTime(){
	std::chrono::high_resolution_clock::time_point end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - time_beg);
	cout << "Execution time : " << duration.count()/1.e3 << " s" << endl;
}

string Displays::center(const string& text, int width) {
	int padding = width - text.length();
	int leftPadding = padding / 2;
	int rightPadding = padding - leftPadding;
	return string(leftPadding, ' ') + text + string(rightPadding, ' ');
}

void Displays::DisplayArray(const vector<vector<string>>& elements, const vector<vector<unsigned int>>& fusion) {
	unsigned int nbRow = fusion.size();
	unsigned int nbCols = 0;
	for(unsigned int i=0;i<fusion[0].size();i++) nbCols += fusion[0][i];
	for(unsigned int i=1;i<nbRow;i++){
	        unsigned int buf = 0;
		for(unsigned int ii=0;ii<fusion[i].size();ii++) buf += fusion[i][ii];
	    	if( buf != nbCols ){
	    		cerr << "The number of columns is not consistent for printing array, aborting" << endl;
	    		exit(EXIT_FAILURE);
	    	}
	}
	
	// compute column width
	unsigned int width = 0;
	for(unsigned int i=0;i<elements.size();i++){
		for(unsigned int j=0;j<elements[i].size();j++){
	    		if( elements[i][j].length() > width ) width = elements[i][j].length();
	    }
	}
	width += 4;
	
	string line = "+";
	for(unsigned int i=0;i<nbCols;i++) line += string(width, '-') + "+";
	
	cout << line << endl;
	for (unsigned int i = 0; i < nbRow; i++) {
		cout << "|";
		for(unsigned int j=0;j<fusion[i].size();j++){
			unsigned int colspan = fusion[i][j];
			unsigned int mergedWidth = width*colspan+(colspan-1);
			cout << center(elements[i][j], mergedWidth) << "|";
		}
		cout << endl << line << endl;
	}
}

void Displays::DisplayArray(const vector<vector<string>>& elements, const vector<vector<unsigned int>>& fusion, ofstream& filetoprint) {
	unsigned int nbRow = fusion.size();
	unsigned int nbCols = 0;
	for(unsigned int i=0;i<fusion[0].size();i++) nbCols += fusion[0][i];
	for(unsigned int i=1;i<nbRow;i++){
	        unsigned int buf = 0;
		for(unsigned int ii=0;ii<fusion[i].size();ii++) buf += fusion[i][ii];
	    	if( buf != nbCols ){
	    		cerr << "The number of columns is not consistent for printing array, aborting" << endl;
	    		exit(EXIT_FAILURE);
	    	}
	}
	
	// compute column width
	unsigned int width = 0;
	for(unsigned int i=0;i<elements.size();i++){
		for(unsigned int j=0;j<elements[i].size();j++){
	    		if( elements[i][j].length() > width ) width = elements[i][j].length();
	    }
	}
	width += 4;
	
	string line = "+";
	for(unsigned int i=0;i<nbCols;i++) line += string(width, '-') + "+";
	
	filetoprint << line << endl;
	for (unsigned int i = 0; i < nbRow; i++) {
		filetoprint << "|";
		for(unsigned int j=0;j<fusion[i].size();j++){
			unsigned int colspan = fusion[i][j];
			unsigned int mergedWidth = width*colspan+(colspan-1);
			filetoprint << center(elements[i][j], mergedWidth) << "|";
		}
		filetoprint << endl << line << endl;
	}
}

void Displays::DisplayOrthogonalCell(Crystal *Crystal){
	string *Dirs1 = new string[3];
	string *Planes1 = new string[3];
	unsigned int width = 0;
	string *i_1_d = new string[3];
	string *i_1_p = new string[3];
	for(unsigned int i=0;i<3;i++){
		i_1_d[i] = "";
		i_1_p[i] = "";
	}
	if( Crystal->getCrystallo() == "Hexagonal" ){
		for(unsigned int i=0;i<3;i++){
			i_1_d[i] = "_"+to_string(-Crystal->getOrthogonalDirs()[i*3]-Crystal->getOrthogonalDirs()[i*3+1]);
			i_1_p[i] = "_"+to_string(-Crystal->getOrthogonalPlanes()[i*3]-Crystal->getOrthogonalPlanes()[i*3+1]);
		}
	}

	for(unsigned int i=0;i<3;i++){
		Dirs1[i] = "["+to_string(Crystal->getOrthogonalDirs()[i*3]);
		Planes1[i] = "("+to_string(Crystal->getOrthogonalPlanes()[i*3]);
		
		Dirs1[i] += "_"+to_string(Crystal->getOrthogonalDirs()[i*3+1]);
		Planes1[i] += "_"+to_string(Crystal->getOrthogonalPlanes()[i*3+1]);
		
		Dirs1[i] += i_1_d[i];
		Planes1[i] += i_1_p[i];
		
		Dirs1[i] += "_"+to_string(Crystal->getOrthogonalDirs()[i*3+2]);
		Planes1[i] += "_"+to_string(Crystal->getOrthogonalPlanes()[i*3+2]);
		
		Dirs1[i] += "]";
		Planes1[i] += ")";
		if( Dirs1[i].length() > width ) width = Dirs1[i].length();
		if( Planes1[i].length() > width ) width = Planes1[i].length();
	}

	width += 12;
	string line = string(width+2, '-');
	string side = "|";
	side += string(width, ' ')+'|';
	cout << "The system is constructed following :" << endl << endl;
	// Display grain 1
	cout << "        " << "   /" << string(width, ' ') << "/" << endl;
	cout << "        " << "  /" << string(width, ' ') << "/" << endl;
	cout << "        " << " /" << center(Planes1[2], width) << "/ " << Planes1[0] << endl;
	cout << "        " << line << endl;
	cout << "        " << "|" << center(Planes1[1], width) << "|" << endl;
	cout << "        " << side << endl;
	cout << "        " << side << endl;
	cout << "        " << "|" << center("", width) << "|        ^ z   " << endl;
	cout << "        " << side << "        |     " << endl;
	cout << "        " << side << "        ---> x" << endl;
	cout << "        " << "|     ^ " << Dirs1[2] << string(width-7-Dirs1[2].length(), ' ') << "|     y /   " << endl;
	cout << "        " << "|     |" << string(width-6,' ') << "|      v      " << endl;
	cout << "        " << "|     ---> " << Dirs1[0] << string(width-10-Dirs1[0].length(), ' ') << "|" << endl;
	cout << "        " << "|    /" << string(width-5,' ') << "|" << endl;
	cout << "        " << "|   v " << Dirs1[1] << string(width-5-Dirs1[1].length(), ' ') << "|  /" << endl;
	cout << "        " << side << " /" << endl;
	cout << "        " << side << "/" << endl;
	cout << "        " << line << endl;


	delete[] Dirs1;
	delete[] Planes1;
	delete[] i_1_d;
	delete[] i_1_p;

}

void Displays::DisplayGB(Crystal *Crystal1, Crystal *Crystal2){
	string *Dirs1 = new string[3];
	string *Dirs2 = new string[3];
	string *Planes1 = new string[3];
	string *Planes2 = new string[3];
	unsigned int width = 0;
	string *i_1_d = new string[3];
	string *i_2_d = new string[3];
	string *i_1_p = new string[3];
	string *i_2_p = new string[3];
	for(unsigned int i=0;i<3;i++){
		i_1_d[i] = "";
		i_2_d[i] = "";
		i_1_p[i] = "";
		i_2_p[i] = "";
	}
	if( Crystal1->getCrystallo() == "Hexagonal" ){
		for(unsigned int i=0;i<3;i++){
			i_1_d[i] = "_"+to_string(-Crystal1->getOrthogonalDirs()[i*3]-Crystal1->getOrthogonalDirs()[i*3+1]);
			i_1_p[i] = "_"+to_string(-Crystal1->getOrthogonalPlanes()[i*3]-Crystal1->getOrthogonalPlanes()[i*3+1]);
		}
	}
	if( Crystal2->getCrystallo() == "Hexagonal" ){
		for(unsigned int i=0;i<3;i++){
			i_2_d[i] = "_"+to_string(-Crystal2->getOrthogonalDirs()[i*3]-Crystal2->getOrthogonalDirs()[i*3+1]);
			i_2_p[i] = "_"+to_string(-Crystal2->getOrthogonalPlanes()[i*3]-Crystal2->getOrthogonalPlanes()[i*3+1]);
		}
	}
	
	for(unsigned int i=0;i<3;i++){
		Dirs1[i] = "["+to_string(Crystal1->getOrthogonalDirs()[i*3]);
		Planes1[i] = "("+to_string(Crystal1->getOrthogonalPlanes()[i*3]);
		Dirs2[i] = "["+to_string(Crystal2->getOrthogonalDirs()[i*3]);
		Planes2[i] = "("+to_string(Crystal2->getOrthogonalPlanes()[i*3]);
		
		Dirs1[i] += "_"+to_string(Crystal1->getOrthogonalDirs()[i*3+1]);
		Planes1[i] += "_"+to_string(Crystal1->getOrthogonalPlanes()[i*3+1]);
		Dirs2[i] += "_"+to_string(Crystal2->getOrthogonalDirs()[i*3+1]);
		Planes2[i] += "_"+to_string(Crystal2->getOrthogonalPlanes()[i*3+1]);
		
		Dirs1[i] += i_1_d[i];
		Planes1[i] += i_1_p[i];
		Dirs2[i] += i_2_d[i];
		Planes2[i] += i_2_p[i];
		
		Dirs1[i] += "_"+to_string(Crystal1->getOrthogonalDirs()[i*3+2]);
		Planes1[i] += "_"+to_string(Crystal1->getOrthogonalPlanes()[i*3+2]);
		Dirs2[i] += "_"+to_string(Crystal2->getOrthogonalDirs()[i*3+2]);
		Planes2[i] += "_"+to_string(Crystal2->getOrthogonalPlanes()[i*3+2]);
		
		Dirs1[i] += "]";
		Planes1[i] += ")";
		Dirs2[i] += "]";
		Planes2[i] += ")";
		if( Dirs1[i].length() > width ) width = Dirs1[i].length();
		if( Dirs2[i].length() > width ) width = Dirs2[i].length();
		if( Planes1[i].length() > width ) width = Planes1[i].length();
		if( Planes2[i].length() > width ) width = Planes2[i].length();
	}

	width += 12;
	string line = string(width+2, '-');
	string side = "|";
	side += string(width, ' ')+'|';
	cout << "The GB is constructed following :" << endl << endl;
	// Display grain 1
	cout << "        " << "   /" << string(width, ' ') << "/" << endl;
	cout << "        " << "  /" << string(width, ' ') << "/" << endl;
	cout << "        " << " /" << center(Planes1[2], width) << "/ " << Planes1[0] << endl;
	cout << "        " << line << endl;
	cout << "        " << "|" << center(Planes1[1], width) << "|" << endl;
	cout << "        " << side << endl;
	cout << "        " << side << endl;
	cout << "        " << "|" << center("Grain 1", width) << "|        ^ z   " << endl;
	cout << "        " << side << "        |     " << endl;
	cout << "        " << side << "        ---> x" << endl;
	cout << "        " << "|     ^ " << Dirs1[2] << string(width-7-Dirs1[2].length(), ' ') << "|     y /   " << endl;
	cout << "        " << "|     |" << string(width-6,' ') << "|      v      " << endl;
	cout << "        " << "|     ---> " << Dirs1[0] << string(width-10-Dirs1[0].length(), ' ') << "|" << endl;
	cout << "        " << "|    /" << string(width-5,' ') << "|" << endl;
	cout << "        " << "|   v " << Dirs1[1] << string(width-5-Dirs1[1].length(), ' ') << "|  /" << endl;
	cout << "        " << side << " /" << endl;
	cout << "        " << side << "/" << endl;
	cout << "        " << line << endl;

	// Display grain 2
	cout << "        " << "   /" << string(width, ' ') << "/" << endl;
	cout << "        " << "  /" << string(width, ' ') << "/" << endl;
	cout << "        " << " /" << center(Planes2[2], width) << "/ " << Planes2[0] << endl;
	cout << "        " << line << endl;
	cout << "        " << "|" << center(Planes2[1], width) << "|" << endl;
	cout << "        " << side << endl;
	cout << "        " << side << endl;
	cout << "        " << "|" << center("Grain 2", width) << "|" << endl;
	cout << "        " << side << endl;
	cout << "        " << side << endl;
	cout << "        " << "|     ^ " << Dirs2[2] << string(width-7-Dirs2[2].length(), ' ') << "|" << endl;
	cout << "        " << "|     |" << string(width-6,' ') << "|" << endl;
	cout << "        " << "|     ---> " << Dirs2[0] << string(width-10-Dirs2[0].length(), ' ') << "|" << endl;
	cout << "        " << "|    /" << string(width-5,' ') << "|" << endl;
	cout << "        " << "|   v " << Dirs2[1] << string(width-5-Dirs2[1].length(), ' ') << "|  /" << endl;
	cout << "        " << side << " /" << endl;
	cout << "        " << side << "/" << endl;
	cout << "        " << line << endl;

	delete[] Dirs1;
	delete[] Dirs2;
	delete[] Planes1;
	delete[] Planes2;
	delete[] i_1_d;
	delete[] i_2_d;
	delete[] i_1_p;
	delete[] i_2_p;
}


void Displays::ProgressBar(unsigned int &Final, unsigned int &current){
	const int bar_length = 30;
	if( current+1 == Final ){
		cout << endl;
		return;
	}
	double prog = double(current+1)/double(Final);
	cout << "\r[" << string(floor(bar_length*prog),'X') << string(ceil(bar_length*(1-prog)),'-') << "] " << setprecision(3) << 100*prog << "%";
}
