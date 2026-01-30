//**********************************************************************************
//*   BasicManips/main.cpp                                                         *
//**********************************************************************************
//* This file contains the implementation of the BasicManips executable.      	   *
//* It allows to do the duplicate, shift and merge options of atomsk, the  	   *
//* difference with atomsk is that this program support bonds, angles and	   *
//* atom_style full								   *
//* This executable also allow to convert files (lmp, data, cfg)
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
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <Displays.h>

using namespace std;

void ExecMsg(){
	cerr << "Usage: ./BasicManips Option ParametersThatDependOnOption OutputFilename" << endl;
	cerr << "It allows to do the duplicate, shift and merge options of atomsk, the difference with atomsk is that this program support bonds, angles and atom_style full" << endl;
	cerr << "Possible options and related parameters :" << endl;
	cerr << "-duplicate InputFilename nx ny nz => duplicate the system nx, ny and nz times in x, y and z directions" << endl;
	cerr << "-shift InputFilename sx sy sz => shift the system of sx, sy and sz Angstrom in x, y and z directions" << endl;
	cerr << "-merge dir N Filename_1 Filename_2 .. Filename_N => merge N atomic systems in a given dir direction (which can be either x, y or z), Filename_i are the filename of the atomic systems to merge" << endl;
	cerr << "This executable also allows to convert files from .lmp (or .cfg) in .cfg (or .lmp) for this option, provide only the name of the inputfile and the name of the outpufile" << endl;
	exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	unsigned int nmerge;
	if( argc == 1 ) ExecMsg();
	string firstarg = argv[1];
	if( firstarg == "-merge" ){
		if( argc < 4 ) ExecMsg();
		istringstream iss_np(argv[3]);
		iss_np >> nmerge;
	}
	if( argc != 7 && argc != 3 && argc != nmerge+5 ) ExecMsg();
	Dis.Printer_OnlyAuxProp();

	if( argc == 7 && firstarg == "-duplicate" ){
		AtomicSystem AtSys(argv[2]);
		unsigned int nx,ny,nz;
		istringstream iss_nx(argv[3]);
		istringstream iss_ny(argv[4]);
		istringstream iss_nz(argv[5]);
		iss_nx >> nx;
		iss_ny >> ny;
		iss_nz >> nz;
		AtSys.duplicate(nx,ny,nz);
		AtSys.printSystem(argv[6]);
	}else if( argc == 7 && firstarg == "-shift" ){
		AtomicSystem AtSys(argv[2]);
		double sx,sy,sz;
		istringstream iss_sx(argv[3]);
		istringstream iss_sy(argv[4]);
		istringstream iss_sz(argv[5]);
		iss_sx >> sx;
		iss_sy >> sy;
		iss_sz >> sz;
		AtSys.ApplyShift(sx,sy,sz);
		AtSys.printSystem(argv[6]);
	}else if( firstarg == "-merge" ){
		AtomicSystem *AtSys = new AtomicSystem[nmerge];
		string mergedir = argv[2];
		for(unsigned int n=0;n<nmerge;n++) AtSys[n].FilenameConstructor(argv[n+4]);
		AtomicSystem newAtSys(AtSys,nmerge,mergedir);
		newAtSys.printSystem(argv[nmerge+4]);
		delete[] AtSys;
	}else if( argc == 3 ){
		AtomicSystem AtSys(argv[1]);
		AtSys.printSystem(argv[2]);
	}else ExecMsg();


	Dis.ExecutionTime();	
	return 0;
}
