//**********************************************************************************
//*   FitAndSaveGMM/main.cpp                                                       *
//**********************************************************************************
//* This file contains the implementation of the FitAndSaveGMM executable.         *
//* It allows to fit and save in the AtomHIC database a gaussian mixture model	   *
//* fitted on given descriptors and using the SGMA method			   *
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
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include "AtomicSystem.h"
#include "MathTools.h"
#include "Crystal.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 2 ){
		cerr << "Usage: ./SymmetrizeSurfaces InputFilename Bottom/Upper OutputFilename" << endl;
		cerr << "InputFilename should be an atomic system file" << endl;
		cerr << "Bottom/Upper specify which surface should be duplicated" << endl;
		Dis.Printer_FitAndSaveGMM();	
		return EXIT_FAILURE;
	}

	Dis.Printer_FitAndSaveGMM();	
	string InputFilename = argv[1];
	string surf2dup = argv[2];
	if( surf2dup != "Bottom" && surf2dup != "Upper" ){
		cerr << "Second argument should be either \"Bottom\" or \"Upper\"" << endl;
		exit(EXIT_FAILURE);
	}
	string OutputFilename = argv[3];

	double slab_width = 10.;
	
	AtomicSystem AtSys(InputFilename);
	unsigned int nbAt = AtSys.getNbAtom();
	
	// search min and max z pos of the system
	double min_z = std::numeric_limits<double>::max();
	double max_z = std::numeric_limits<double>::min();
	for(unsigned int i=0;i<nbAt;i++){
		if( AtSys.getAtom(i).pos.z < min_z ) min_z = AtSys.getAtom(i).pos.z;
		if( AtSys.getAtom(i).pos.z > max_z ) max_z = AtSys.getAtom(i).pos.z;
	}
	if( max_z-min_z < 3*slab_width ) cout << "Warning the system is very thick, it may cause issues" << endl;

	// construct subsystems to make coincide
	vector<Atom> SubSys_ref; // reference one (the one to be duplicated)
	vector<Atom> SubSys_work;
	vector<unsigned int> work_ind;
	for(unsigned int i=0;i<nbAt;i++){
		if( surf2dup == "Bottom" ){
			if( AtSys.getAtom(i).pos.z < min_z+slab_width ) SubSys_ref.push_back(AtSys.getAtom(i));
			else if( AtSys.getAtom(i).pos.z > max_z-(2.*slab_width) ){
				SubSys_work.push_back(AtSys.getAtom(i));
				work_ind.push_back(i);
			}
		}else{
			if( AtSys.getAtom(i).pos.z > max_z-slab_width ) SubSys_ref.push_back(AtSys.getAtom(i));
			else if( AtSys.getAtom(i).pos.z < min_z+(2.*slab_width) ){
				SubSys_work.push_back(AtSys.getAtom(i));
				work_ind.push_back(i);
			}
		}
	}
	//cout << "Nb atom in work sys = " << SubSys_work.size() << ", in ref sys = " << SubSys_ref.size() << endl;

	// Apply rotation to the ref system
	for(unsigned int i=0;i<SubSys_ref.size();i++) SubSys_ref[i].pos.z *= -1.;
	MathTools *MT = new MathTools;
	double *vec = new double[3];
	vec[0] = 0.;
	vec[1] = 0.;
	vec[2] = 1.;
	double theta = M_PI;
	double *rotmat = new double[9];
	MT->Vec2rotMat(vec,theta,rotmat);

	double small_shift = 0.1;
	for(unsigned int i=0;i<SubSys_ref.size();i++){
		MT->MatDotAt(rotmat,SubSys_ref[i],SubSys_ref[i]);
		SubSys_ref[i].pos.x += small_shift;
		SubSys_ref[i].pos.y += small_shift;
	}
	for(unsigned int i=0;i<SubSys_work.size();i++){
		SubSys_work[i].pos.x += small_shift;
		SubSys_work[i].pos.y += small_shift;
	}

	// wrap the two systems
	double x,y,z;
	for(unsigned int i=0;i<SubSys_ref.size();i++){
		// compute reduced coordinates
		x = SubSys_ref[i].pos.x*AtSys.getG1()[0]+SubSys_ref[i].pos.y*AtSys.getG2()[0]+SubSys_ref[i].pos.z*AtSys.getG3()[0];
		y = SubSys_ref[i].pos.x*AtSys.getG1()[1]+SubSys_ref[i].pos.y*AtSys.getG2()[1]+SubSys_ref[i].pos.z*AtSys.getG3()[1];
		z = SubSys_ref[i].pos.x*AtSys.getG1()[2]+SubSys_ref[i].pos.y*AtSys.getG2()[2]+SubSys_ref[i].pos.z*AtSys.getG3()[2];
		if( x >= 1. || x < 0. ) x = x-floor(x);
		if( y >= 1. || y < 0. ) y = y-floor(y);
		// cartesian coordinates
		SubSys_ref[i].pos.x = x*AtSys.getH1()[0]+y*AtSys.getH2()[0]+z*AtSys.getH3()[0];
		SubSys_ref[i].pos.y = x*AtSys.getH1()[1]+y*AtSys.getH2()[1]+z*AtSys.getH3()[1];
	}
	for(unsigned int i=0;i<SubSys_work.size();i++){
		// compute reduced coordinates
		x = SubSys_work[i].pos.x*AtSys.getG1()[0]+SubSys_work[i].pos.y*AtSys.getG2()[0]+SubSys_work[i].pos.z*AtSys.getG3()[0];
		y = SubSys_work[i].pos.x*AtSys.getG1()[1]+SubSys_work[i].pos.y*AtSys.getG2()[1]+SubSys_work[i].pos.z*AtSys.getG3()[1];
		z = SubSys_work[i].pos.x*AtSys.getG1()[2]+SubSys_work[i].pos.y*AtSys.getG2()[2]+SubSys_work[i].pos.z*AtSys.getG3()[2];
		if( x >= 1. || x < 0. ) x = x-floor(x);
		if( y >= 1. || y < 0. ) y = y-floor(y);
		// cartesian coordinates
		SubSys_work[i].pos.x = x*AtSys.getH1()[0]+y*AtSys.getH2()[0]+z*AtSys.getH3()[0];
		SubSys_work[i].pos.y = x*AtSys.getH1()[1]+y*AtSys.getH2()[1]+z*AtSys.getH3()[1];
	}


	// search translation vectors making systems coinciding (this may allow also to find cell vectors)
	double tol = 1.e-1;
	vector<double> TransVecs;
	vector<double> TransVecsNorm;
	vector<unsigned int> count_TransVecs;
	
	// search for first atom of work sys to construct the vectors and then parallelized for all atoms
	//for(unsigned int j=0;j<SubSys_ref.size();j++){
	//	double dx = SubSys_ref[j].pos.x - SubSys_work[0].pos.x;
	//	double dy = SubSys_ref[j].pos.y - SubSys_work[0].pos.y;
	//	double dz = SubSys_ref[j].pos.z - SubSys_work[0].pos.z;
	//	double dist = dx*dx + dy*dy + dz*dz;
	//	TransVecs.push_back(dx);
	//	TransVecs.push_back(dy);
	//	TransVecs.push_back(dz);
	//	TransVecsNorm.push_back(dist);
	//	count_TransVecs.push_back(1);
	//}
	//
	//unsigned int nbAt_work = SubSys_work.size();
	//for(unsigned int i=1;i<SubSys_work.size();i++){
	//	Dis.ProgressBar(nbAt_work,i);
	//	for(unsigned int j=0;j<SubSys_ref.size();j++){
	//		double dx = SubSys_ref[j].pos.x - SubSys_work[i].pos.x;
	//		double dy = SubSys_ref[j].pos.y - SubSys_work[i].pos.y;
	//		double dz = SubSys_ref[j].pos.z - SubSys_work[i].pos.z;
	//		double dist = dx*dx + dy*dy + dz*dz;
	//		bool already = false;
	//		for(unsigned int k=0;k<TransVecsNorm.size();k++){
	//			if( fabs(TransVecsNorm[k]-dist) < tol ){
	//				already = true;
	//				count_TransVecs[k]++;
	//				break;
	//			}
	//		}
	//	}
	//}

	// better to do brut force to be sure to find the right place
	unsigned int nbAt_work = SubSys_work.size();
	for(unsigned int i=0;i<SubSys_work.size();i++){
		Dis.ProgressBar(nbAt_work,i);
		for(unsigned int j=0;j<SubSys_ref.size();j++){
			double dx = SubSys_ref[j].pos.x - SubSys_work[i].pos.x;
			double dy = SubSys_ref[j].pos.y - SubSys_work[i].pos.y;
			double dz = SubSys_ref[j].pos.z - SubSys_work[i].pos.z;
			double dist = dx*dx + dy*dy + dz*dz;
			bool already = false;
			for(unsigned int k=0;k<TransVecsNorm.size();k++){
				if( fabs(TransVecsNorm[k]-dist) < tol ){
					already = true;
					count_TransVecs[k]++;
					break;
				}
			}
			if( !already ){
				TransVecs.push_back(dx);
				TransVecs.push_back(dy);
				TransVecs.push_back(dz);
				TransVecsNorm.push_back(dist);
				count_TransVecs.push_back(1);
			}
		}
	}

	// shift the ref system and verify that the ref system fully coincide with the working sys when apply x and y bc
	bool coincide = false;
	unsigned int max_test = 10;
	unsigned int count = 0;
	double *at_sp = new double[SubSys_ref.size()];
	while( !coincide ){
		coincide = true;
		unsigned int opt_ind = MT->max(count_TransVecs);
		for(unsigned int j=0;j<SubSys_ref.size();j++){
			SubSys_ref[j].pos.x -= TransVecs[opt_ind*3];
			SubSys_ref[j].pos.y -= TransVecs[opt_ind*3+1];
			SubSys_ref[j].pos.z -= TransVecs[opt_ind*3+2];
		}

		double sigma = 0.25;
		double min_sp = 1.e-1;
		tol = 1.e-3;
		for(unsigned int i=0;i<SubSys_ref.size();i++){
			at_sp[i] = 0.;
			double xr = SubSys_ref[i].pos.x;
			double yr = SubSys_ref[i].pos.y;
			double zr = SubSys_ref[i].pos.z;
			for(unsigned int j=0;j<SubSys_work.size();j++){
				for(int bx=-1;bx<2;bx++){
					for(int by=-1;by<2;by++){
						double xw = SubSys_work[j].pos.x + AtSys.getH1()[0]*bx + AtSys.getH2()[0]*by;
						double yw = SubSys_work[j].pos.y + AtSys.getH1()[1]*bx + AtSys.getH2()[1]*by;
						double zw = SubSys_work[j].pos.z;
						at_sp[i] += MT->gaussian(xw, yw, zw, xr, yr, zr, sigma);
					}
				}
			}
			if( ( i != 0 && fabs(at_sp[i]-at_sp[i-1]) > tol ) || at_sp[i] < min_sp ){
				coincide = false;
				break;
			}
		}
		if( !coincide ){
			if( count >= max_test ) coincide = true;
			else{
				// remove the shift and erase the corresponding elements in vectors
				for(unsigned int j=0;j<SubSys_ref.size();j++){
					SubSys_ref[j].pos.x += TransVecs[opt_ind*3];
					SubSys_ref[j].pos.y += TransVecs[opt_ind*3+1];
					SubSys_ref[j].pos.z += TransVecs[opt_ind*3+2];
				}
				count_TransVecs.erase(count_TransVecs.begin()+opt_ind);
				TransVecs.erase(TransVecs.begin()+opt_ind);
				TransVecs.erase(TransVecs.begin()+opt_ind);
				TransVecs.erase(TransVecs.begin()+opt_ind);
				TransVecsNorm.erase(TransVecsNorm.begin()+opt_ind);
				count++;
			}
		}
	}
	if( count >= max_test ){
		cerr << "Cannot align systems, aborting.." << endl;
		exit(EXIT_FAILURE);
	}
	
	// Now remove ions of the working sys not coinciding with the ref sys
	// keep ions above/bellow mean z pos of ref system
	vector<unsigned int> index2rm;
	double mean_z_pos = 0.;
	for(unsigned int i=0;i<SubSys_ref.size();i++) mean_z_pos += SubSys_ref[i].pos.z;
	mean_z_pos /= SubSys_ref.size();

	for(unsigned int i=0;i<SubSys_work.size();i++){
		if( surf2dup == "Bottom" && SubSys_work[i].pos.z < mean_z_pos ) continue;
		else if( surf2dup == "Upper" && SubSys_work[i].pos.z > mean_z_pos ) continue;
		else{
			double sigma = 0.25;
			double min_sp = 1.e-1;
			double cur_at_sp = 0.;
			double xr = SubSys_work[i].pos.x;
			double yr = SubSys_work[i].pos.y;
			double zr = SubSys_work[i].pos.z;
			for(unsigned int j=0;j<SubSys_ref.size();j++){
				for(int bx=-1;bx<2;bx++){
					for(int by=-1;by<2;by++){
						double xw = SubSys_ref[j].pos.x + AtSys.getH1()[0]*bx + AtSys.getH2()[0]*by;
						double yw = SubSys_ref[j].pos.y + AtSys.getH1()[1]*bx + AtSys.getH2()[1]*by;
						double zw = SubSys_ref[j].pos.z;
						cur_at_sp += MT->gaussian(xw, yw, zw, xr, yr, zr, sigma);
					}
				}
			}
			if( cur_at_sp < min_sp ) index2rm.push_back(work_ind[i]);
		}
	}
	cout << "Removing " << index2rm.size() << " atoms to symmetrize the surfaces" << endl;
	AtSys.RemoveAtoms(index2rm);

	

	//Atom *AtList_ref = new Atom[SubSys_ref.size()];
	//Atom *AtList_work = new Atom[SubSys_work.size()];
	//Atom *AtList_tot = new Atom[SubSys_work.size()+SubSys_ref.size()];
	//for(unsigned int i=0;i<SubSys_ref.size();i++){
	//	AtList_ref[i] = SubSys_ref[i];
	//	AtList_tot[i] = SubSys_ref[i];
	//}
	//for(unsigned int i=0;i<SubSys_work.size();i++){
	//	AtList_work[i] = SubSys_work[i];
	//	AtList_tot[i+SubSys_ref.size()] = SubSys_work[i];
	//}
	//string crystalName = "Forsterite";
	//Crystal *MyC = new Crystal(crystalName);
	//AtomicSystem AtSys_ref(AtList_ref,SubSys_ref.size(),MyC,AtSys.getH1(),AtSys.getH2(),AtSys.getH3());
	//AtSys_ref.setAux(at_sp,"at_sp");
	//AtSys_ref.printSystem_aux("RefSys_2.lmp","at_sp");
	//AtomicSystem AtSys_work(AtList_work,SubSys_work.size(),MyC,AtSys.getH1(),AtSys.getH2(),AtSys.getH3());
	//AtSys_work.printSystem("WorkSys.lmp");
	//AtomicSystem AtSys_tot(AtList_tot,SubSys_work.size()+SubSys_ref.size(),MyC,AtSys.getH1(),AtSys.getH2(),AtSys.getH3());
	//AtSys_tot.printSystem("TotSys.lmp");





	//vector<unsigned int> opt_ind(4);
	//// find the 4 first more represented vectors
	//cout << "searching min " << count_TransVecs.size() << endl;
	//for(int k=0;k<4;k++){
	//	opt_ind[k] = 0;
	//	for(unsigned int i=1;i<TransVecsNorm.size();i++){
	//		if( count_TransVecs[i] > count_TransVecs[opt_ind[k]] ){
	//			bool already = false;
	//			for(int j=0;j<k;j++){
	//				cout << j << endl;
	//				if( i == opt_ind[j] ){
	//					already = true;
	//					break;
	//				}
	//			}
	//			if( !already ) opt_ind[k] = i;
	//		}
	//	}
	//}

	//for(unsigned int k=0;k<4;k++){
	//	cout << "Find vecs (" << opt_ind[k] << ", nb coincide = " << count_TransVecs[opt_ind[k]] << "): " << setprecision(6) << TransVecs[opt_ind[k]*3] << " " << TransVecs[opt_ind[k]*3+1] << " " << TransVecs[opt_ind[k]*3+2] << endl;
	//}	




	AtSys.printSystem(OutputFilename);
	Dis.ExecutionTime();	
	return 0;
}
