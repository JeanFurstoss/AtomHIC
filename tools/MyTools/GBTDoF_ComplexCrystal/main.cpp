//**********************************************************************************
//*   GBTDoF_ComplexCrystal/main.cpp                                               *
//**********************************************************************************
//* This file contains the implementation of the GBTDoF_ComplexCrystal executable. *
//* it allows to generate atomic systems sampling translational DoF along CSL      *
//* for complex crystals (multispecies and charged)				   *
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
//*	- implement control on h_p_x, k_p_x, l_p_x				   *
//**********************************************************************************
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include "AtomicSystem.h"
#include "MathTools.h"
#include "Bicrystal.h"
#include "Crystal.h"
#include "ComputeAuxiliary.h"
#include <Displays.h>

using namespace std;

void ExecMsg(){
	cerr << "Usage: ./GBTDoF_ComplexCrystal DiscretizationMode nx ny h_RotAxis k_RotAxis (i_RotAxis) l_RotAxis RotAngle(in degree) h_GBPlane k_GBPlane (i_GBPlane) l_GBPlane CrystalName(has to be defined in /data/Crystal/) Rationalize" << endl;
	cerr << "This executable creates atomic systems containing a GB with a given misorientation and GB plane and varying the translational (i.e. shift) degree of freedom" << endl << endl;
	cerr << "DiscretizationMode can be \"NbPts\" or \"Distance\":" << endl;
	cerr << "\t - NbPts => nx, ny specify the number of shifts (sampling points) along each of the two vectors defining the GB plane" << endl;
	cerr << "\t - Distance =>  nx, ny specify the distances between two shifts in each direction (then the number of shifts will be computed from these distances)" << endl << endl;
	cerr << "The i Miller indexes (for rotation axis and GB plane) should only be used if the crystal is hexagonal" << endl;
	cerr << "Rationalize can be either 0 or 1" << endl;
	cerr << "1 => rationalize the GB (i.e. search the closest CSL GB to the provided parameters, in this case the fixed parameters for CSL calculation can be important, they are read from FixedParameters.ath file if exist if not defaults values are used)" << endl;
	cerr << "0 => do not rationalize the GB" << endl;
	cerr << "Each configuration is saved as a separate dump file (e.g., GB_Shift_0_0.lmp or GB_DSC_Shift_0_0_0.lmp or GB_CSL_Shift_0_0_0.lmp, etc), also containing the values of the applied shift" << endl;
	cerr << "In addition the program will also return 3 dump files containing the CSL lattice and the two grains (the two latters can be used to change shift between crystals), if vacuum = 1 Grain1_vacuum.lmp and Grain2_vacuum.lmp will be also generated" << endl;
	
	exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc != 15 && argc != 13 ) ExecMsg();

        string DisMode = argv[1];
	if( DisMode != "NbPts" && DisMode != "Distance" ){
		cout << "\t *** Error : Unknown DiscretizationMode !!! ***" << endl << endl;
		ExecMsg();
	}

	Dis.Printer_SampleGB_GammaSurface();

	double z_step = 0.05; // step for searching different surfaces

	// get parameters
	int h_a, k_a ,l_a, h_p, k_p, l_p, i_a, i_p;
	double theta;
	unsigned int rat;
	unsigned int nx, ny;
	double dx, dy;
	bool rat_b;
	string crystalName;
	unsigned int current_readind = 2;
	if( DisMode == "NbPts" ){
		istringstream iss_nx(argv[current_readind]);
		iss_nx >> nx;
		current_readind++;
		istringstream iss_ny(argv[current_readind]); 
		iss_ny >> ny;
		current_readind++;
	}else{
		istringstream iss_dx(argv[current_readind]);
		iss_dx >> dx;
		current_readind++;
		istringstream iss_dy(argv[current_readind]); 
		iss_dy >> dy;
		current_readind++;
	}

	istringstream iss_ha(argv[current_readind]);
	iss_ha >> h_a;
	current_readind++;
	istringstream iss_ka(argv[current_readind]);
	iss_ka >> k_a;
	current_readind++;
	if( argc == 15 ){
		istringstream iss_ia(argv[current_readind]);
		iss_ia >> i_a;
		current_readind++;
		if( i_a != (-h_a-k_a) ) cout << "Warning, i index of rotation axis is different than -h-k, i will be considered as equal to " << -h_a-k_a << endl;
	}
	istringstream iss_la(argv[current_readind]);
	iss_la >> l_a;
	current_readind++;
	istringstream iss_theta(argv[current_readind]);
	iss_theta >> theta;
	current_readind++;
	theta *= M_PI/180.;
	istringstream iss_hp(argv[current_readind]);
	iss_hp >> h_p;
	current_readind++;
	istringstream iss_kp(argv[current_readind]);
	iss_kp >> k_p;
	current_readind++;
	if( argc == 15 ){
		istringstream iss_ip(argv[current_readind]);
		iss_ip >> i_p;
		current_readind++;
		if( i_p != (-h_p-k_p) ) cout << "Warning, i index of GB plane is different than -h-k, i will be considered as equal to " << -h_p-k_p << endl;
	}
	istringstream iss_lp(argv[current_readind]);
	iss_lp >> l_p;
	current_readind++;
	istringstream iss_cn(argv[current_readind]);
	iss_cn >> crystalName;
	current_readind++;
	istringstream iss_rat(argv[current_readind]);
	iss_rat >> rat;
	current_readind++;
	if( rat == 0 ) rat_b = false;
	else rat_b = true;
	
	
	// set minimum aside to very small value to be sure that we only have the UC of the GB
	vector<string> Properties;
	Properties.push_back("MIN_BOX_ASIDE 2.");
	Properties.push_back("FULL_GRAINS 0");
	
	// set the initial z shift of upper and lower grain to 0 
	Properties.push_back("MOTIF_SHIFT 0 0 0");
	Properties.push_back("LOWER_GRAIN_SHIFT 0 0 0");
	Properties.push_back("UPPER_GRAIN_SHIFT 0 0 0");
	Bicrystal *refGB = new Bicrystal(crystalName,h_a,k_a,l_a,theta,h_p,k_p,l_p,rat_b,Properties,0,0,0,true);

	if( DisMode == "Distance" ){
		nx = round(refGB->getH1()[0]/dx);
		ny = round(refGB->getH2()[1]/dy);
		if( nx < 1 ){
			cout << "Warning from the provided dx value nx is null => nx has been set equal to 1" << endl;
			nx = 1;
		}
		if( ny < 1 ){
			cout << "Warning from the provided dy value ny is null => ny has been set equal to 1" << endl;
			ny = 1;
		}
	}
	double true_dx = refGB->getH1()[0] / ((double) nx);
	double true_dy = refGB->getH2()[1] / ((double) ny);
	
	vector<int> GBPlane1, GBPlane2;
	vector<int> xPlane1, xPlane2;
	for(unsigned int d=0;d<3;d++){
		GBPlane1.push_back(refGB->getCrystal()->getOrthogonalPlanes()[6+d]);
		GBPlane2.push_back(refGB->getCrystal2()->getOrthogonalPlanes()[6+d]);
		xPlane1.push_back(refGB->getCrystal()->getOrthogonalPlanes()[d]);
		xPlane2.push_back(refGB->getCrystal2()->getOrthogonalPlanes()[d]);
	}
	bool SymmGB = true;
	for(unsigned int i=0;i<3;i++){
		if( fabs(GBPlane1[i]) != fabs(GBPlane2[i]) ){
			SymmGB = false;
			break;
		}
	}
	double dhkl1 = refGB->getCrystal()->ComputeD_hkl(GBPlane1[0],GBPlane1[1],GBPlane1[2]);
	double dhkl2 = refGB->getCrystal2()->ComputeD_hkl(GBPlane2[0],GBPlane2[1],GBPlane2[2]);
	//unsigned int nb_z_step1 = ceil( dhkl1 / z_step );
	unsigned int nb_z_step1 = 150;
	double true_z_step1 = dhkl1 / nb_z_step1;
	//unsigned int nb_z_step2 = ceil( dhkl2 / z_step );
	unsigned int nb_z_step2 = 150;
	double true_z_step2 = dhkl2 / nb_z_step2;
	double zero = 0.;
	ComputeAuxiliary CA;
	
	cout << "\t * * * Searching number of different surfaces of upper grain (z <=> (" << GBPlane1[0] << " " << GBPlane1[1] << " " << GBPlane1[2] << ") x <=> (" << xPlane1[0] << " " << xPlane1[1] << " " << xPlane1[2] << ")) * * *" << endl;
	vector<Crystal*> Cryst1(nb_z_step1);
	vector<unsigned int> CrystWithDiffSurf1;
	vector<double> z_shift_1;
	for(unsigned int i=0;i<nb_z_step1;i++){
		Dis.ProgressBar(nb_z_step1,i);
		Cryst1[i] = new Crystal(crystalName);
		Cryst1[i]->ReadProperties(Properties);
		Cryst1[i]->RotateCrystal(GBPlane1[0],GBPlane1[1],GBPlane1[2],xPlane1[0],xPlane1[1],xPlane1[2]);
		Cryst1[i]->ShiftMotif(zero,zero,i*true_z_step1);
		Cryst1[i]->ConstructOrthogonalCell();
		Cryst1[i]->getOrientedSystem()->MakeSurfaceNeutral(false);
		bool already = false;
		for(unsigned int j=0;j<CrystWithDiffSurf1.size();j++){
			if( CA.AreSurfacesSame(Cryst1[i]->getOrientedSystem(),"down",Cryst1[CrystWithDiffSurf1[j]]->getOrientedSystem(),"down") ){
				already = true;
				break;
			}
		}
		if( !already ){
		       CrystWithDiffSurf1.push_back(i);
		       //Cryst1[i]->getOrientedSystem()->print_lmp("G1_debug_"+to_string(z_shift_1.size())+".lmp");
	       	       z_shift_1.push_back(i*true_z_step1);
		}
	}
	cout << "\t\tDone ! " << z_shift_1.size() << " different neutral surfaces have been found for upper grain" << endl;
	// symmetrize the upper grain free surface and print the systems
	cout << "\t * * * Symmetrizing the free surfaces of the different upper grain systems * * *" << endl;
	vector<unsigned int> surf2rm;
	unsigned int count = 0;
	for(unsigned int i=0;i<z_shift_1.size();i++){
		if( CA.AreSurfacesSame(Cryst1[CrystWithDiffSurf1[i]]->getOrientedSystem(),"up",Cryst1[CrystWithDiffSurf1[i]]->getOrientedSystem(),"down") ){
			//Cryst1[CrystWithDiffSurf1[i]]->getOrientedSystem()->print_lmp("G1_"+to_string(count)+".lmp");
			count++;
		}else{
			if( Cryst1[CrystWithDiffSurf1[i]]->getOrientedSystem()->SymmetrizeSurfaces("up") ){
				//Cryst1[CrystWithDiffSurf1[i]]->getOrientedSystem()->print_lmp("G1_"+to_string(count)+".lmp");
				count++;
			}else{
				cout << "One free surface cannot be symmetrized, removing it from this system for the exploration of TDoF" << endl;
				surf2rm.push_back(i);
			}
		}
	}
	for(unsigned int i=0;i<surf2rm.size();i++){
		z_shift_1.erase(z_shift_1.begin()+surf2rm[surf2rm.size()-1-i]);
		CrystWithDiffSurf1.erase(CrystWithDiffSurf1.begin()+surf2rm[surf2rm.size()-1-i]);
	}
	cout << "\t\tDone ! " << endl;

	cout << "\t * * * Searching number of different surfaces of lower grain (z <=> (" << GBPlane2[0] << " " << GBPlane2[1] << " " << GBPlane2[2] << ") x <=> (" << xPlane2[0] << " " << xPlane2[1] << " " << xPlane2[2] << ")) * * *" << endl;
	vector<Crystal*> Cryst2(nb_z_step2);
	vector<unsigned int> CrystWithDiffSurf2;
	vector<double> z_shift_2;
	for(unsigned int i=0;i<nb_z_step2;i++){
		Dis.ProgressBar(nb_z_step2,i);
		Cryst2[i] = new Crystal(crystalName);
		Cryst2[i]->ReadProperties(Properties);
		//Cryst2[i]->RotateCrystal(GBPlane2[0],GBPlane2[1],GBPlane2[2],xPlane2[0],xPlane2[1],xPlane2[2]);
		Cryst2[i]->RotateCrystal(refGB->getRotMatG2());
		Cryst2[i]->ShiftMotif(zero,zero,i*true_z_step2);
		Cryst2[i]->ConstructOrthogonalCell();
		Cryst2[i]->getOrientedSystem()->MakeSurfaceNeutral(false);
		bool already = false;
		for(unsigned int j=0;j<CrystWithDiffSurf2.size();j++){
			if( CA.AreSurfacesSame(Cryst2[i]->getOrientedSystem(),"down",Cryst2[CrystWithDiffSurf2[j]]->getOrientedSystem(),"down") ){
				already = true;
				break;
			}
		}
		if( !already ){
		       CrystWithDiffSurf2.push_back(i);
		       //Cryst2[i]->getOrientedSystem()->print_lmp("G2_debug_"+to_string(z_shift_2.size())+".lmp");
	       	       z_shift_2.push_back(i*true_z_step2);
		}
	}
	cout << "\t\tDone ! " << z_shift_2.size() << " different neutral surfaces have been found for lower grain" << endl;
	// symmetrize the lower grain free surface and print the systems
	cout << "\t * * * Symmetrizing the free surfaces of the different lower grain systems * * *" << endl;
	surf2rm.clear();
	count = 0;
	for(unsigned int i=0;i<z_shift_2.size();i++){
		if( CA.AreSurfacesSame(Cryst2[CrystWithDiffSurf2[i]]->getOrientedSystem(),"up",Cryst2[CrystWithDiffSurf2[i]]->getOrientedSystem(),"down") ){
			//Cryst2[CrystWithDiffSurf2[i]]->getOrientedSystem()->print_lmp("G2_"+to_string(count)+".lmp");
			count++;
		}else{
			if( Cryst2[CrystWithDiffSurf2[i]]->getOrientedSystem()->SymmetrizeSurfaces("down") ){
				//Cryst2[CrystWithDiffSurf2[i]]->getOrientedSystem()->print_lmp("G2_"+to_string(count)+".lmp");
				count++;
			}else{
				cout << "One free surface cannot be symmetrized, removing it from this system for the exploration of TDoF" << endl;
				surf2rm.push_back(i);
			}
		}
	}
	for(unsigned int i=0;i<surf2rm.size();i++){
		z_shift_2.erase(z_shift_2.begin()+surf2rm[surf2rm.size()-1-i]);
		CrystWithDiffSurf2.erase(CrystWithDiffSurf2.begin()+surf2rm[surf2rm.size()-1-i]);
	}
	cout << "\t\tDone ! " << endl;

	vector<Crystal*> *CrystUp;
	vector<Crystal*> *CrystDown;
	vector<unsigned int> *indexesUp;
	vector<unsigned int> *indexesDown;
	bool inv = false;
	if( SymmGB ){
		bool inv = false;
		if( z_shift_1.size() != z_shift_2.size() ){
			cout << "Warning, the GB is symmetric but we found a different number of free surfaces for the 2 grains, it may lead to imcomplete TDoF exploration" << endl;
		}
		// rearrange order of the system having the highest number of surface to have the same surfaces of the other system
		if( z_shift_1.size() >= z_shift_2.size() ){
			CrystUp = &Cryst1;
			CrystDown = &Cryst2;
			indexesUp = &CrystWithDiffSurf1;
			indexesDown = &CrystWithDiffSurf2;
		}else{
			CrystUp = &Cryst2;
			CrystDown = &Cryst1;
			indexesUp = &CrystWithDiffSurf2;
			indexesDown = &CrystWithDiffSurf1;
			inv = true;
		}
		cout << "Searching correspondances between surfaces of G1 and G2 (as GB is symmetric)" << endl;
		vector<int> temp_ind((*indexesDown).size(),-1);
		vector<unsigned int> topushback;
		for(unsigned int i1=0;i1<(*indexesUp).size();i1++){
			bool found = false;
			for(unsigned int i2=0;i2<(*indexesDown).size();i2++){
				if( CA.AreSurfacesSame((*CrystUp)[(*indexesUp)[i1]]->getOrientedSystem(),"down",(*CrystDown)[(*indexesDown)[i2]]->getOrientedSystem(),"up",true) ){
					if( temp_ind[i2] != -1 ) cout << "WARNING, system found twice!" << endl;
					temp_ind[i2] = (*indexesUp)[i1];
					found = true;
					break;
				}
			}
			if( !found ) topushback.push_back((*indexesUp)[i1]);
		}
		// verify that all surfaces of sys 2 have been found
		(*indexesUp).clear();
		for(unsigned int i=0;i<(*indexesDown).size();i++)
			if( temp_ind[i] != -1 )	(*indexesUp).push_back(temp_ind[i]);
		if( (*indexesUp).size() != (*indexesDown).size() ) cout << "WARNING issues when comparing free surfaces of upper and lower grain" << endl;
		for(unsigned int i=0;i<topushback.size();i++) (*indexesUp).push_back(topushback[i]);
	}else{
		CrystUp = &Cryst1;
		CrystDown = &Cryst2;
		indexesUp = &CrystWithDiffSurf1;
		indexesDown = &CrystWithDiffSurf2;
	}

	// print the different systems (with free surfaces) and then apply the misfit, the duplication and reduce the z cell size
	double currentMx, currentMy, currentDupX, currentDupY;
	if( inv ){
		currentMx = refGB->getMx1();
		currentMy = refGB->getMy1();
		currentDupX = refGB->getDupX1();
		currentDupY = refGB->getDupY1();
	}else{
		currentMx = refGB->getMx2();
		currentMy = refGB->getMy2();
		currentDupX = refGB->getDupX2();
		currentDupY = refGB->getDupY2();
	}
	for(unsigned int il=0;il<(*indexesDown).size();il++){
		AtomicSystem *currentAtSys = (*CrystDown)[(*indexesDown)[il]]->getOrientedSystem();
		double min_z = std::numeric_limits<double>::max();
        	double max_z = -std::numeric_limits<double>::max();
		unsigned int current_nbat = currentAtSys->getNbAtom();
		for(unsigned int i=0;i<current_nbat;i++){
			if( currentAtSys->getAtomRef(i).pos.z > max_z ) max_z = currentAtSys->getAtomRef(i).pos.z;
			if( currentAtSys->getAtomRef(i).pos.z < min_z ) min_z = currentAtSys->getAtomRef(i).pos.z;
		}
		// set large H3 cell size and shift the system to be in the middle of the box
		currentAtSys->getH3()[2] *= 1.5;
		double shift_z = ((currentAtSys->getH3()[2]-(max_z-min_z))/2.)-min_z;
		currentAtSys->ApplyShift(zero,zero,shift_z);
		currentAtSys->print_lmp("LowerGrain_"+to_string(il)+".lmp");
		// set H3 cell size to exactly fit the size of the system and shift it to fit in the box
		shift_z = -((currentAtSys->getH3()[2]-(max_z-min_z))/2.);
		currentAtSys->getH3()[2] = max_z-min_z;
		currentAtSys->ApplyShift(zero,zero,shift_z);
		// Apply the misfit
		for(unsigned int i=0;i<current_nbat;i++){
			currentAtSys->getAtomRef(i).pos.x *= currentMx;
			currentAtSys->getAtomRef(i).pos.y *= currentMy;
		}
		currentAtSys->getH1()[0] *= currentMx;
		currentAtSys->getH1()[1] *= currentMy;
		currentAtSys->getH2()[0] *= currentMx;
		currentAtSys->getH2()[1] *= currentMy;
		// duplicate
		currentAtSys->duplicate(currentDupX,currentDupY,1);
	}
	if( inv ){
		currentMx = refGB->getMx2();
		currentMy = refGB->getMy2();
		currentDupX = refGB->getDupX2();
		currentDupY = refGB->getDupY2();
	}else{
		currentMx = refGB->getMx1();
		currentMy = refGB->getMy1();
		currentDupX = refGB->getDupX1();
		currentDupY = refGB->getDupY1();
	}
	for(unsigned int iu=0;iu<(*indexesUp).size();iu++){
		AtomicSystem *currentAtSys = (*CrystUp)[(*indexesUp)[iu]]->getOrientedSystem();
		double min_z = std::numeric_limits<double>::max();
        	double max_z = -std::numeric_limits<double>::max();
		unsigned int current_nbat = currentAtSys->getNbAtom();
		for(unsigned int i=0;i<current_nbat;i++){
			if( currentAtSys->getAtomRef(i).pos.z > max_z ) max_z = currentAtSys->getAtomRef(i).pos.z;
			if( currentAtSys->getAtomRef(i).pos.z < min_z ) min_z = currentAtSys->getAtomRef(i).pos.z;
		}
		// set large H3 cell size and shift the system to be in the middle of the box
		currentAtSys->getH3()[2] *= 1.5;
		double shift_z = ((currentAtSys->getH3()[2]-(max_z-min_z))/2.)-min_z;
		currentAtSys->ApplyShift(zero,zero,shift_z);
		currentAtSys->print_lmp("UpperGrain_"+to_string(iu)+".lmp");
		// set H3 cell size to exactly fit the size of the system and shift it to fit in the box
		shift_z = -((currentAtSys->getH3()[2]-(max_z-min_z))/2.);
		currentAtSys->getH3()[2] = max_z-min_z;
		currentAtSys->ApplyShift(zero,zero,shift_z);
		// Apply the misfit
		for(unsigned int i=0;i<current_nbat;i++){
			currentAtSys->getAtomRef(i).pos.x *= currentMx;
			currentAtSys->getAtomRef(i).pos.y *= currentMy;
		}
		currentAtSys->getH1()[0] *= currentMx;
		currentAtSys->getH1()[1] *= currentMy;
		currentAtSys->getH2()[0] *= currentMx;
		currentAtSys->getH2()[1] *= currentMy;
		// duplicate
		currentAtSys->duplicate(currentDupX,currentDupY,1);
	}

	// print bulk system
	if( inv ) refGB->print_Grains(false,"LowerGrainBulk.lmp","UpperGrainBulk.lmp");	
	else refGB->print_Grains(false,"UpperGrainBulk.lmp","LowerGrainBulk.lmp");	
	
	// Now create the GBs
	double rem_shift_x = -((double) nx)*true_dx;
	double rem_shift_y = -((double) ny)*true_dy;
	for(unsigned int il=0;il<(*indexesDown).size();il++){
		unsigned int iu_lo = 0;
		if( SymmGB ) iu_lo = il;
		for(unsigned int sx=0;sx<nx;sx++){
			(*CrystDown)[(*indexesDown)[il]]->getOrientedSystem()->ApplyShift(true_dx,zero,zero);
			for(unsigned int sy=0;sy<ny;sy++){
				(*CrystDown)[(*indexesDown)[il]]->getOrientedSystem()->ApplyShift(zero,true_dy,zero);
				for(unsigned int iu=iu_lo;iu<(*indexesUp).size();iu++){
					refGB->PasteGrains((*CrystUp)[(*indexesUp)[iu]]->getOrientedSystem(),(*CrystDown)[(*indexesDown)[il]]->getOrientedSystem());
					refGB->print_lmp("GB_PlaneLow_"+to_string(il)+"_PlaneUp_"+to_string(iu)+"_Shift_"+to_string(sx)+"_"+to_string(sy)+".lmp");
				}
			}
			(*CrystDown)[(*indexesDown)[il]]->getOrientedSystem()->ApplyShift(zero,rem_shift_y,zero);
		}
		(*CrystDown)[(*indexesDown)[il]]->getOrientedSystem()->ApplyShift(rem_shift_x,zero,zero);
	}



	Dis.ExecutionTime();
	delete refGB;
	for(unsigned int i=0;i<nb_z_step1;i++) delete Cryst1[i];
	for(unsigned int i=0;i<nb_z_step2;i++) delete Cryst2[i];

	return 0;
}
