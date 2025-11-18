//**********************************************************************************
//*   BondOriParam/main.cpp                                                        *
//**********************************************************************************
//* This file contains the implementation of the BondOriParam executable.          *
//* It allows to compute the bond orientational parameter (Furstoss et al 2024) of *
//* a given crystal								   *
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
#include <Crystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include <chrono>
#include <omp.h>
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	
	string InputFilename, CrystalType, outfilename;
	if( argc == 4 ){
		InputFilename = argv[1];
		CrystalType = argv[2];
		outfilename = argv[3];
		Crystal MyCrystal(CrystalType);
		if( MyCrystal.getIsReferenceBondOriParam() ) Dis.Printer_NoFixedParams();
		else Dis.Printer_BondOriParam();
	}else if( argc == 3 ){
		InputFilename = argv[1];
		outfilename = argv[2];
		Dis.Printer_BondOriParam();
	}else{
		cerr << "Usage: ./BondOriParam AtomicInputFilename CrystalType(The crystal type has to be defined in /data/Crystal/) outputFilename" << endl;
		cerr << "or: ./BondOriParam AtomicInputFilename outputFilename" << endl;
		cerr << "In the latter case, the crystal will be considered as a monosite crystal and the value of spherical harmonic degree and cutoff radius and SteinhardtStyle will be read in FixedParameters.ath file (if the file does not exist the default values will be used)" << endl; 
		cerr << "The first case could be used for multisite and non-centrosymmetric crystals for which the reference bond orientational parameters have been stored in the AtomHic database using the ./SaveNonCSCrystalBondOriParam executable" << endl;
		cerr << "In the first case, if the crystal is not multisite the values of cutoff radius will be computed from the cell parameters of crystal and the spherical harmonics degree and SteinhardtStyle will be read in FixedParameters.ath file if it exists" << endl;
		cerr << "More details on the computation of this order parameter can be found in Furstoss et al. (2024) Comp. Mat. Science." << endl;
		return EXIT_FAILURE;
	}
	AtomicSystem MySystem(InputFilename);
	if( argc == 4 ) MySystem.setCrystal(CrystalType);
	ComputeAuxiliary CA(&MySystem);
	MySystem.setAux(CA.BondOrientationalParameter(), "BondOriParam");
	MySystem.printSystem_aux(outfilename, "BondOriParam");
	
	Dis.ExecutionTime();	
	return 0;
}
