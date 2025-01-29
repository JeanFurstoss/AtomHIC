//**********************************************************************************
//*   AnalyzeCSL/main.cpp                                                          *
//**********************************************************************************
//* This file contains the implementation of the AnalyzeCSL executable.            *
//* It allows to compute the CSL of a given misorientation giving the crystal	   *
//* structure 									   *
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
//*	- output a dump file with the CSL and DSC lattices                         *
//*	- add option for choosing ouputs					   *
//**********************************************************************************

#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "MathTools.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	// check the number of argument
	if( argc < 7 ){
		cerr << "Usage: AnalyzeCSL RotAxis_X_Component RotAxis_Y_Component RotAxis_Z_Component Rot_Angle CrystalType Output_Filename" << endl;
		cerr << "This executable compute the CSL lattice for a given misorientation and crystal type (which has to be defined in /data/Crystal/)" << endl;
	        cerr << "It will print in Output_Filename : SigmaValue SallestCSLVector ShapeFactorCSLBasis_1 ShapeFactorCSLBasis_2 c1 c2 c3" << endl;
		cerr << "Where ShapeFactorCSLBasis_1 = LargestCSLVector / SmallestCSLVector" << endl;
		cerr << "Where ShapeFactorCSLBasis_2 = LargestCSLVector^2 / Product of the two other CSLVector" << endl;
		cerr << "c1, c2, c3 are the three CSL basis vectors" << endl;
		cerr << "The algorithm used for this calculation is based on a iterative procedure developped by Bonnet and Rolland 1975 and permitting to compute near-CSL for low symmetry crystals (for these type of crystals, the numerical parameters defined for CSL in /data/FixedParameters/FixedParameters.dat can be important)" << endl;
		return EXIT_FAILURE;
	}
	double ra_x, ra_y, ra_z, rot_angle;
	istringstream iss_ra_x(argv[1]);
	iss_ra_x >> ra_x;
	istringstream iss_ra_y(argv[2]);
	iss_ra_y >> ra_y;
	istringstream iss_ra_z(argv[3]);
	iss_ra_z >> ra_z;
	istringstream iss_rot_ang(argv[4]);
	iss_rot_ang >> rot_angle;
	string CrystalType = argv[5];
	string OutputFilename = argv[6];
	
	// compute the nearest [hkl] (with integer h, k and l)
	MathTools *MT = new MathTools;
	Crystal *MyCrystal = new Crystal(CrystalType);
	double *rot_Mil = new double[3];
	int *MilInd = new int[3];
	double *Rat_RotAx = new double[3];
	int *CSL_vec = new int[3];
	bool found, CSL_found;
	double Rat_RotAng;
	double tolIntVec = 1e-1;
	unsigned int hkl_max = 25000;
	rot_Mil[0] = ra_x;// / MyCrystal->getALength()[0];
	rot_Mil[1] = ra_y;// / MyCrystal->getALength()[1];
	rot_Mil[2] = ra_z;// / MyCrystal->getALength()[2];
	MT->find_integer_vector(rot_Mil,tolIntVec,hkl_max,MilInd,found);
	if( found ){
		// Create a bicrystal object as it has CSL tools
		Bicrystal *MyBicrystal = new Bicrystal(CrystalType, MilInd[0], MilInd[1], MilInd[2], rot_angle);
		// Rationalize the misorientation axis and angle and get a CSL vector
		Rat_RotAng = MyBicrystal->RationalizeOri(MilInd[0], MilInd[1], MilInd[2], rot_angle, Rat_RotAx, CSL_vec);
		// Compute the CSL basis
		CSL_found = MyBicrystal->searchCSL(Rat_RotAx, Rat_RotAng, CSL_vec, 1);
		if( CSL_found ){
			// compute interesting values : smallest CSL vector, shape factor of the basis
			double SmallestL, SF_1, SF_2;
			vector<double> list;
			for(unsigned int i=0;i<3;i++){
				list.push_back(0.);
				for(unsigned int j=0;j<3;j++) list[i] += pow(MyBicrystal->getCSL_Basis()[j*3+i],2.);
				list[i] = sqrt(list[i]);
			}
			MT->sort(list,0,1,list);
			SmallestL = list[0];
			SF_1 = list[2]/list[0];
			SF_2 = list[2]*list[2]/(list[0]*list[1]);
			// write the sigma value and the CSL basis
			ofstream writefile(OutputFilename);
			writefile << MyBicrystal->getSigma() << " " << SmallestL << " " << SF_1 << " " << SF_2;
	        	for(unsigned int i=0;i<9;i++) writefile << " " << MyBicrystal->getCSL_Basis()[i];
			writefile.close();
		}else{
			ofstream writefile(OutputFilename);
			writefile << "Failed";
			writefile.close();
			cout << "No CSL basis has been found" << endl;
		}
		
		delete MyBicrystal;
	}else{
		ofstream writefile(OutputFilename);
		writefile << "Failed";
		writefile.close();
		cout << "No integer vector has been found corresponding to this rotation angle" << endl;
	}
	
	delete[] rot_Mil;
	delete[] MilInd;
	delete[] Rat_RotAx;
	delete[] CSL_vec;
	delete MyCrystal;
	delete MT;
	Dis.ExecutionTime();	
	return 0;
}
