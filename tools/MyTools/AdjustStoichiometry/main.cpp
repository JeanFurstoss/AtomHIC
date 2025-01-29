//**********************************************************************************
//*   AtomicStrain/main.cpp                                                        *
//**********************************************************************************
//* This file contains the implementation of the AtomicStrain executable.          *
//* It allows to compute the atomic strain of given dump files relative to a 	   *
//* reference dump file								   *
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
//*	- make a generic function in AtomicSystem to adjust stoichiometry and put  *
//* this executable in generic tools						   *
//**********************************************************************************


#include <AtomicSystem.h>
#include <Bicrystal.h>
#include <Crystal.h>
#include <ComputeAuxiliary.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include "MyStructs.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	Crystal *MyCrystal = new Crystal("Forsterite");
	
	string InputFilename, OutputFilename;
	double SiO_d;
	double SecFac = 2.;
	bool AuxPrint = false;
	if( argc >= 3 ){
		InputFilename = argv[1];
		OutputFilename = argv[2];
		if( argc > 3 ) AuxPrint = true;
	}else{
		cerr << "Usage: ./AdjustStoichiometry InputFilename OutputFilename AuxName2Print_1 AuxName2Print_2.." << endl;
		cerr << "Adjsut the stoichiometry of Mg2SiO4 system by removing Si, Mg and O ions, if we have to remove Si ions we assume that SiO4 tetrahedra have not been separated" << endl;
		cerr << "TODO AuxName2Print implementation" << endl;
		return EXIT_FAILURE;
	}
	AtomicSystem MySystem(InputFilename);
	MySystem.printSystem_aux("TEST.cfg","vx");
	const unsigned int nbAt = MySystem.getNbAtom();

	vector<string> AtTypes;
	for(unsigned int i=0;i<3;i++) AtTypes.push_back("");

	unsigned int type_uint_O, type_uint_Si, type_uint_Mg;
	for(unsigned int i=0;i<3;i++){
		if( MySystem.getAtomType(i) == "Mg" ){
			type_uint_Mg = i+1;
			AtTypes[i] = "Mg";
		} else if( MySystem.getAtomType(i) == "Si" ){
			type_uint_Si = i+1;
			AtTypes[i] = "Si";
		} else if( MySystem.getAtomType(i) == "O" ){
			type_uint_O = i+1;
			AtTypes[i] = "O";
		}
	}
	unsigned int nbO(0), nbSi(0), nbMg(0);
	vector<unsigned int> Index_O, Index_Si, Index_Mg;
	for(unsigned int i=0;i<nbAt;i++) {
		if( MySystem.getAtom(i).type_uint == type_uint_Mg ) Index_Mg.push_back(i);
		else if( MySystem.getAtom(i).type_uint == type_uint_Si ) Index_Si.push_back(i);
		else if( MySystem.getAtom(i).type_uint == type_uint_O ) Index_O.push_back(i);
	}
	nbO = Index_O.size();
	nbMg = Index_Mg.size();
	nbSi = Index_Si.size();

	int NbExtraMg, NbExtraO, NbExtraSi;
	NbExtraMg = nbMg - nbSi*2;
	NbExtraO = nbO - nbSi*4;
	if( NbExtraMg < 0 ){
		if( nbMg % 2 != 0 ) NbExtraMg = 1;
		else NbExtraMg = 0;
	}
	if( NbExtraO < 0 ) NbExtraO = 0;
	NbExtraSi = nbSi - (nbMg-NbExtraMg)/2.;
	if( NbExtraSi < 0 ) NbExtraSi = 0;

	MathTools MT;
	MySystem.searchNeighbours(4.5);
	const unsigned int nbNMax = MySystem.getNbMaxN();
	double dmin, xpos, ypos, zpos, xp, yp, zp, dtemp;
	unsigned int id;
	vector<double> MgList, OList, SiList;
	vector<unsigned int> RmO(NbExtraO), RmSi(NbExtraSi), RmMg(NbExtraMg);
	if( NbExtraMg > 0 ){
		for(unsigned int i=0;i<nbMg;i++) {
			xpos = MySystem.getWrappedPos(Index_Mg[i]).x;
			ypos = MySystem.getWrappedPos(Index_Mg[i]).y;
			zpos = MySystem.getWrappedPos(Index_Mg[i]).z;
			dmin = 3.;
			for(unsigned int j=0;j<MySystem.getNeighbours(Index_Mg[i]*(nbNMax+1));j++){
				id = MySystem.getNeighbours(Index_Mg[i]*(nbNMax+1)+j+1);
				xp = MySystem.getWrappedPos(id).x+MySystem.getCLNeighbours(Index_Mg[i]*nbNMax*3+j*3)*MySystem.getH1()[0]+MySystem.getCLNeighbours(Index_Mg[i]*nbNMax*3+j*3+1)*MySystem.getH2()[0]+MySystem.getCLNeighbours(Index_Mg[i]*nbNMax*3+j*3+2)*MySystem.getH3()[0]-xpos;
				yp = MySystem.getWrappedPos(id).y+MySystem.getCLNeighbours(Index_Mg[i]*nbNMax*3+j*3)*MySystem.getH1()[1]+MySystem.getCLNeighbours(Index_Mg[i]*nbNMax*3+j*3+1)*MySystem.getH2()[1]+MySystem.getCLNeighbours(Index_Mg[i]*nbNMax*3+j*3+2)*MySystem.getH3()[1]-ypos;
				zp = MySystem.getWrappedPos(id).z+MySystem.getCLNeighbours(Index_Mg[i]*nbNMax*3+j*3)*MySystem.getH1()[2]+MySystem.getCLNeighbours(Index_Mg[i]*nbNMax*3+j*3+1)*MySystem.getH2()[2]+MySystem.getCLNeighbours(Index_Mg[i]*nbNMax*3+j*3+2)*MySystem.getH3()[2]-zpos;
				dtemp = sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)); 
				if( dtemp < dmin ) dmin = dtemp;
			}
			if( dmin < 2.5 ){
				MgList.push_back(Index_Mg[i]);
				MgList.push_back(dmin);
			}
		}
		MT.sort(MgList,1,2,MgList);
		for(unsigned int i=0;i<NbExtraMg;i++) RmMg[i] = (unsigned int) MgList[i*2];
	}

	if( NbExtraO > 0 ){
		for(unsigned int i=0;i<nbO;i++) {
			xpos = MySystem.getWrappedPos(Index_O[i]).x;
			ypos = MySystem.getWrappedPos(Index_O[i]).y;
			zpos = MySystem.getWrappedPos(Index_O[i]).z;
			dmin = 3.;
			for(unsigned int j=0;j<MySystem.getNeighbours(Index_O[i]*(nbNMax+1));j++){
				id = MySystem.getNeighbours(Index_O[i]*(nbNMax+1)+j+1);
				xp = MySystem.getWrappedPos(id).x+MySystem.getCLNeighbours(Index_O[i]*nbNMax*3+j*3)*MySystem.getH1()[0]+MySystem.getCLNeighbours(Index_O[i]*nbNMax*3+j*3+1)*MySystem.getH2()[0]+MySystem.getCLNeighbours(Index_O[i]*nbNMax*3+j*3+2)*MySystem.getH3()[0]-xpos;
				yp = MySystem.getWrappedPos(id).y+MySystem.getCLNeighbours(Index_O[i]*nbNMax*3+j*3)*MySystem.getH1()[1]+MySystem.getCLNeighbours(Index_O[i]*nbNMax*3+j*3+1)*MySystem.getH2()[1]+MySystem.getCLNeighbours(Index_O[i]*nbNMax*3+j*3+2)*MySystem.getH3()[1]-ypos;
				zp = MySystem.getWrappedPos(id).z+MySystem.getCLNeighbours(Index_O[i]*nbNMax*3+j*3)*MySystem.getH1()[2]+MySystem.getCLNeighbours(Index_O[i]*nbNMax*3+j*3+1)*MySystem.getH2()[2]+MySystem.getCLNeighbours(Index_O[i]*nbNMax*3+j*3+2)*MySystem.getH3()[2]-zpos;
				dtemp = sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)); 
				if( dtemp < dmin ) dmin = dtemp;
			}
			if( dmin < 2.5 ){
				OList.push_back(Index_O[i]);
				OList.push_back(dmin);
			}
		}
		MT.sort(OList,1,2,OList);
		for(unsigned int i=0;i<NbExtraO;i++) RmO[i] = (unsigned int) OList[i*2];
	}

	if( NbExtraSi > 0 ){
		for(unsigned int i=0;i<nbSi;i++) {
			xpos = MySystem.getWrappedPos(Index_Si[i]).x;
			ypos = MySystem.getWrappedPos(Index_Si[i]).y;
			zpos = MySystem.getWrappedPos(Index_Si[i]).z;
			dmin = 3.;
			for(unsigned int j=0;j<MySystem.getNeighbours(Index_Si[i]*(nbNMax+1));j++){
				id = MySystem.getNeighbours(Index_Si[i]*(nbNMax+1)+j+1);
				xp = MySystem.getWrappedPos(id).x+MySystem.getCLNeighbours(Index_Si[i]*nbNMax*3+j*3)*MySystem.getH1()[0]+MySystem.getCLNeighbours(Index_Si[i]*nbNMax*3+j*3+1)*MySystem.getH2()[0]+MySystem.getCLNeighbours(Index_Si[i]*nbNMax*3+j*3+2)*MySystem.getH3()[0]-xpos;
				yp = MySystem.getWrappedPos(id).y+MySystem.getCLNeighbours(Index_Si[i]*nbNMax*3+j*3)*MySystem.getH1()[1]+MySystem.getCLNeighbours(Index_Si[i]*nbNMax*3+j*3+1)*MySystem.getH2()[1]+MySystem.getCLNeighbours(Index_Si[i]*nbNMax*3+j*3+2)*MySystem.getH3()[1]-ypos;
				zp = MySystem.getWrappedPos(id).z+MySystem.getCLNeighbours(Index_Si[i]*nbNMax*3+j*3)*MySystem.getH1()[2]+MySystem.getCLNeighbours(Index_Si[i]*nbNMax*3+j*3+1)*MySystem.getH2()[2]+MySystem.getCLNeighbours(Index_Si[i]*nbNMax*3+j*3+2)*MySystem.getH3()[2]-zpos;
				dtemp = sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.));
				if( dtemp < dmin ) dmin = dtemp;
			}
			if( dmin < 2.5 ){
				SiList.push_back(Index_Si[i]);
				SiList.push_back(dmin);
			}
		}
		MT.sort(SiList,1,2,SiList);
		vector<double> O2Rm;
		for(unsigned int i=0;i<NbExtraSi;i++){
			RmSi[i] = (unsigned int) SiList[i*2];
			xpos = MySystem.getWrappedPos(RmSi[i]).x;
			ypos = MySystem.getWrappedPos(RmSi[i]).y;
			zpos = MySystem.getWrappedPos(RmSi[i]).z;
			vector<double>().swap(O2Rm);
			for(unsigned int j=0;j<MySystem.getNeighbours(RmSi[i]*(nbNMax+1));j++){
				id = MySystem.getNeighbours(RmSi[i]*(nbNMax+1)+j+1);
				if( MySystem.getAtom(id).type_uint == type_uint_O ){        	
					xp = MySystem.getWrappedPos(id).x+MySystem.getCLNeighbours(RmSi[i]*nbNMax*3+j*3)*MySystem.getH1()[0]+MySystem.getCLNeighbours(RmSi[i]*nbNMax*3+j*3+1)*MySystem.getH2()[0]+MySystem.getCLNeighbours(RmSi[i]*nbNMax*3+j*3+2)*MySystem.getH3()[0]-xpos;
					yp = MySystem.getWrappedPos(id).y+MySystem.getCLNeighbours(RmSi[i]*nbNMax*3+j*3)*MySystem.getH1()[1]+MySystem.getCLNeighbours(RmSi[i]*nbNMax*3+j*3+1)*MySystem.getH2()[1]+MySystem.getCLNeighbours(RmSi[i]*nbNMax*3+j*3+2)*MySystem.getH3()[1]-ypos;
					zp = MySystem.getWrappedPos(id).z+MySystem.getCLNeighbours(RmSi[i]*nbNMax*3+j*3)*MySystem.getH1()[2]+MySystem.getCLNeighbours(RmSi[i]*nbNMax*3+j*3+1)*MySystem.getH2()[2]+MySystem.getCLNeighbours(RmSi[i]*nbNMax*3+j*3+2)*MySystem.getH3()[2]-zpos;
					O2Rm.push_back(id);
					O2Rm.push_back(sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)));
				}
			}
			MT.sort(O2Rm,1,2,O2Rm);
			unsigned int NbOStored = 0;
			unsigned int currentO_index = 0;
			while( NbOStored != 4 ){
				unsigned int sizeO = RmO.size();
				bool already = false;
				for(unsigned int k=0;k<sizeO;k++){
					if( (unsigned int) O2Rm[currentO_index*2] == RmO[k] ){
						already = true;
						break;
					}
				}
				if( !already ){
					RmO.push_back((unsigned int) O2Rm[currentO_index*2]);
					NbOStored += 1;
				}
				currentO_index += 1;
			}
			NbExtraO += 4;
		}
	}
	cout << "To adjust the stoichiometry of the system we will remove " << NbExtraMg << " Mg, " << NbExtraSi << " Si and " << NbExtraO << " O ions" << endl;
	Atom *NewAtoms = new Atom[nbAt-NbExtraO-NbExtraMg-NbExtraSi];
	double *Aux = new double[nbAt-NbExtraO-NbExtraMg-NbExtraSi];
	unsigned int buffer;
	unsigned int grainId_ind = MySystem.getAuxIdAndSize("grainID",buffer);
	bool Store;
	unsigned int count = 0;
	for(unsigned int i=0;i<nbAt;i++){
		if( MySystem.getAtom(i).type_uint == type_uint_Mg ){
			Store = true;
			for(unsigned int j=0;j<RmMg.size();j++){
				if( i == RmMg[j] ){
					Store = false;
					RmMg.erase(RmMg.begin()+j);
					break;
				}
			}
			if( Store ){
				NewAtoms[count] = MySystem.getAtom(i);
				if( AuxPrint ) Aux[count] = MySystem.getAux(grainId_ind)[i];
				count++;
			}
		}else if( MySystem.getAtom(i).type_uint == type_uint_Si ){
			Store = true;
			for(unsigned int j=0;j<RmSi.size();j++){
				if( i == RmSi[j] ){
					Store = false;
					RmSi.erase(RmSi.begin()+j);
					break;
				}
			}
			if( Store ){
				NewAtoms[count] = MySystem.getAtom(i);
				if( AuxPrint ) Aux[count] = MySystem.getAux(grainId_ind)[i];
				count++;
			}
		}else if( MySystem.getAtom(i).type_uint == type_uint_O ){
			Store = true;
			for(unsigned int j=0;j<RmO.size();j++){
				if( i == RmO[j] ){
					Store = false;
					RmO.erase(RmO.begin()+j);
					break;
				}
			}
			if( Store ){
				NewAtoms[count] = MySystem.getAtom(i);
				if( RmO.size() == 129 ){
					cout << "stop" << endl;
				}
				if( AuxPrint ) Aux[count] = MySystem.getAux(grainId_ind)[i];
				count++;
			}
		}
	}
	cout << "Number of atom in the new system : " << count << endl;

	AtomicSystem NewSys(NewAtoms,nbAt-NbExtraO-NbExtraMg-NbExtraSi,MyCrystal,MySystem.getH1(),MySystem.getH2(),MySystem.getH3());
	if( AuxPrint ){
		NewSys.setAux(Aux,"grainID");	
		NewSys.printSystem_aux(OutputFilename,"grainID");
	}else{
		NewSys.printSystem(OutputFilename);
	}
	cout << "New system printed" << endl;
	Dis.ExecutionTime();	
	return 0;
}
