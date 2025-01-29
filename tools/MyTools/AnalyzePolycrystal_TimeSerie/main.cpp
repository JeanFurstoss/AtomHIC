//**********************************************************************************
//*   AnalyzePolycrystal_TimeSerie/main.cpp		                           *
//**********************************************************************************
//* This file contains the implementation of the AnalyzePolycrystal_TimeSerie      *
//* executable.	   								   *
//* It is used to analyze a time serie of polycrystaline deformation simulation    *
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
//*	- use DBScan to identify grains						   *
//**********************************************************************************

#include <filesystem>
#include <AtomicSystem.h>
#include <Bicrystal.h>
#include <Crystal.h>
#include <ComputeAuxiliary.h>
#include <MathTools.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include "MyStructs.h"
#include "Descriptors.h"
#include "GaussianMixtureModel.h"
#include <filesystem>
#include <dirent.h>
#include <Displays.h>

using namespace std;

// For this I think it better to decrease tolerance for GMM fit to 1e-4 and increase the number of KMean initialization to 400 to have more reproductible results
void ReadFile_NbClustAndGB_Init(string Filename_NbClustAndGB_Init, vector<vector<unsigned int>> &NbClustAndGB_Init, vector<double> &GBId_ForNbClust_Init){
	ifstream file(Filename_NbClustAndGB_Init); // The file shouldn't have header and must be organized as: GBId NbClust_1 NbClust_2 NbGB
	if( file ){
		string line;
		unsigned int nbClust_1, nbClust_2, nbGBInit;
		double GBIdInit;
		unsigned int current_GB(0);
		NbClustAndGB_Init.push_back(vector<unsigned int>());
		getline(file,line);
		istringstream text_init(line);
		text_init >> GBIdInit >> nbClust_1 >> nbClust_2 >> nbGBInit;
		GBId_ForNbClust_Init.push_back(GBIdInit);
		NbClustAndGB_Init[current_GB].push_back(nbClust_1);
		NbClustAndGB_Init[current_GB].push_back(nbClust_2);
		NbClustAndGB_Init[current_GB].push_back(nbGBInit);
		while(getline(file,line)){
			istringstream text(line);
			text >> GBIdInit >> nbClust_1 >> nbClust_2 >> nbGBInit;
			if( GBIdInit != GBId_ForNbClust_Init[current_GB] ){
				current_GB++;
				NbClustAndGB_Init.push_back(vector<unsigned int>());
				GBId_ForNbClust_Init.push_back(GBIdInit);
			}
			NbClustAndGB_Init[current_GB].push_back(nbClust_1);
			NbClustAndGB_Init[current_GB].push_back(nbClust_2);
			NbClustAndGB_Init[current_GB].push_back(nbGBInit);
		}
	}else{
		cerr << "We cannot open " << Filename_NbClustAndGB_Init << endl;
		exit(EXIT_FAILURE);
	}
}

bool read_GTFile(string GTFilename, vector<unsigned int> &IdForGrains){
	ifstream file_i(GTFilename, ios::in);
	string line;
	if( file_i ){
		getline(file_i,line);
		while(file_i){
			getline(file_i,line);
			IdForGrains.push_back(0);
		}
	}else{
		cout << "cannot open grain tag file" << endl;
		return false;
	}
	ifstream file(GTFilename, ios::in);
	getline(file,line);
	unsigned int IdIon, vecRef;
	for(unsigned int i=0;i<IdForGrains.size();i++){
		getline(file,line);
		istringstream text(line);
		text >> IdIon >> vecRef;
		unsigned int idGrain = floor(vecRef/10.);
		unsigned int cellParam = vecRef-idGrain*10;
		IdForGrains[(idGrain-1)*4+cellParam-1] = IdIon;
	}
	return true;
			
}

void ComputeAndPrintCellVectors(AtomicSystem &MySystem, string GTFilename, string timestep){
	
	// TODO make global variables
	double tolspangle = 5.;
	double tol_sp = cos(tolspangle*M_PI/180.);
	
	vector<unsigned int> IdIonsGrain;
	bool ok = read_GTFile(GTFilename, IdIonsGrain);
	if( ok ) cout << "Grain tag file read successfully" << endl;
	else cout << "Issue when reading grain tag file" << endl;
	unsigned int nbGrains = IdIonsGrain.size()/4;
	double *CellVectors = new double[9*nbGrains];
	for(unsigned int i=0;i<nbGrains;i++){
		for(unsigned int dim=0;dim<3;dim++){
			CellVectors[i*9+dim*3] = MySystem.getAtom(IdIonsGrain[i*4+dim+1]).pos.x - MySystem.getAtom(IdIonsGrain[i*4]).pos.x;
			CellVectors[i*9+dim*3+1] = MySystem.getAtom(IdIonsGrain[i*4+dim+1]).pos.y - MySystem.getAtom(IdIonsGrain[i*4]).pos.y;
			CellVectors[i*9+dim*3+2] = MySystem.getAtom(IdIonsGrain[i*4+dim+1]).pos.z - MySystem.getAtom(IdIonsGrain[i*4]).pos.z;
		}
		// verify that the three vectors are almost normal
		double na(0.),nb(0.),nc(0.),sp(0.);
		for(unsigned int dim=0;dim<3;dim++){
			na += pow(CellVectors[i*9+dim],2.);
			nb += pow(CellVectors[i*9+3+dim],2.);
			nc += pow(CellVectors[i*9+6+dim],2.);
		}
		na = sqrt(na);
		nb = sqrt(nb);
		nc = sqrt(nc);
		// a,b
		for(unsigned int k=0;k<3;k++) sp += CellVectors[i*9+k]*CellVectors[i*9+3+k];
		if( sp/(na*nb) > tol_sp ) cout << "Warning is seems that a and b are not normals in grain " << i+1 << endl;
		// a,c
		sp = 0.;
		for(unsigned int k=0;k<3;k++) sp += CellVectors[i*9+k]*CellVectors[i*9+6+k];
		if( sp/(na*nc) > tol_sp ) cout << "Warning is seems that a and c are not normals in grain " << i+1 << endl;
		// b,c
		sp = 0.;
		for(unsigned int k=0;k<3;k++) sp += CellVectors[i*9+3+k]*CellVectors[i*9+6+k];
		if( sp/(nb*nc) > tol_sp ) cout << "Warning is seems that b and c are not normals in grain " << i+1 << endl;
	}
		
	string namef1 = "GrainCellVec_";
	string nameext = ".dat";
	std::ofstream f_out_cv(namef1+timestep+nameext);
	for(unsigned int i=0;i<nbGrains;i++){
		for(unsigned int dim1=0;dim1<3;dim1++){
			f_out_cv << CellVectors[i*9+dim1*3];
			for(unsigned int dim2=1;dim2<3;dim2++) f_out_cv << " " << CellVectors[i*9+dim1*3+dim2];
			f_out_cv << endl;
		}
	}
	f_out_cv.close();
	delete[] CellVectors;
}

void getAuxIdAndSize(AtomicSystem &MySystem, unsigned int &size_Struct, unsigned int &Struct_ind, unsigned int &size_GBId, unsigned int &GBId_ind, unsigned int &size_AtVol, unsigned int &AtVol_ind, unsigned int &size_Stress, unsigned int &Stress_ind, unsigned int &size_AtStrain, unsigned int &AtStrain_ind){
	Struct_ind = MySystem.getAuxIdAndSize("Struct",size_Struct);
	GBId_ind = MySystem.getAuxIdAndSize("GBId",size_GBId);
	AtVol_ind = MySystem.getAuxIdAndSize("AtVol",size_AtVol);
	Stress_ind = MySystem.getAuxIdAndSize("c_s",size_Stress);
	AtStrain_ind = MySystem.getAuxIdAndSize("AtomicStrain",size_AtStrain);
}

void ConstructGBIdAndGBIonsArrays(AtomicSystem &MySystem, vector<double> &GBId_arr, vector<vector<unsigned int>> &GBIons, unsigned int &size_Struct, unsigned int &Struct_ind, unsigned int &size_GBId, unsigned int &GBId_ind, unsigned int &Struct_GB, unsigned int &Struct_Amorph){
	cout << "Searching GB ions" << endl;
	// Clear arrays
	GBId_arr.clear();
	vector<double>().swap(GBId_arr);
	for(unsigned int i=0;i<GBIons.size();i++){
		GBIons[i].clear();
		vector<unsigned int>().swap(GBIons[i]);
	}
	GBIons.clear();
	vector<vector<unsigned int>>().swap(GBIons);
	
	unsigned int nbAt = MySystem.getNbAtom();
	for(unsigned int i=0;i<nbAt;i++) {
		bool IsGB = false;
		bool IsAmorph = false;
		unsigned int GBOrAmorph_1;
		unsigned int GBOrAmorph_2;
		if( (unsigned int) MySystem.getAux(Struct_ind)[i*size_Struct] == Struct_GB ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
			IsAmorph = false,
			IsGB = true;
		}else if( (unsigned int) MySystem.getAux(Struct_ind)[i*size_Struct] == Struct_Amorph ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
			IsAmorph = true,
			IsGB = true;
		}
		if( IsGB ){
			bool IsAlreadyStored = false;
			// We use to store the GBId as G1.G2 with G1<G2 (if it is already the case we are in interface 1, in the other case we are in interface 2)
			// TODO warning this method should not work if a grain index is higher than 10
			double IntPart_GBId = floor(MySystem.getAux(GBId_ind)[i*size_GBId]);
			double FloatPart_GBId = MySystem.getAux(GBId_ind)[i*size_GBId]-IntPart_GBId;
			double current_GBId;
			unsigned int Int1Or2OrAmorph;
			if( IntPart_GBId < 10.*FloatPart_GBId ){
				current_GBId = IntPart_GBId+FloatPart_GBId;
				if( IsAmorph ) Int1Or2OrAmorph = 2;
				else Int1Or2OrAmorph = 0;
			}else{
				current_GBId = (FloatPart_GBId*10.)+IntPart_GBId/10.;
				if( IsAmorph ) Int1Or2OrAmorph = 2;
				else Int1Or2OrAmorph = 1;
			}
			//cout << current_GBId << endl;
			for(unsigned int g=0;g<GBId_arr.size();g++){
				if( fabs(current_GBId-GBId_arr[g]) < 1e-5 ){
					GBIons[g*3+Int1Or2OrAmorph].push_back(i);
					IsAlreadyStored = true;
					break;
				}
			}
			if( !IsAlreadyStored ){
				GBId_arr.push_back(current_GBId);
				GBIons.push_back(vector<unsigned int>());
				GBIons.push_back(vector<unsigned int>());
				GBIons.push_back(vector<unsigned int>());
				GBIons[(GBId_arr.size()-1)*3+Int1Or2OrAmorph].push_back(i);
			}
		}
	}
	cout << "Done" << endl;
}

void ShiftAndWrapAmorph(AtomicSystem &MySystem, vector<vector<unsigned int>> &GBIons, vector<vector<double>> &coords, unsigned int &g, vector<double> &current_TransVec){
	cout << "TRANSVEC" << endl;
	for(unsigned int dim=0;dim<3;dim++) cout << current_TransVec[dim] << " ";
	cout << endl;
	for(unsigned int i=0;i<GBIons[g*3+2].size();i++){
		coords.push_back(vector<double>());
		coords[i].push_back(MySystem.getWrappedPos(GBIons[g*3+2][i]).x+current_TransVec[0]);
		coords[i].push_back(MySystem.getWrappedPos(GBIons[g*3+2][i]).y+current_TransVec[1]);
		coords[i].push_back(MySystem.getWrappedPos(GBIons[g*3+2][i]).z+current_TransVec[2]);
	}
	double *TempVec = new double[3];
	for(unsigned int i=0;i<GBIons[g*3+2].size();i++){
		// Compute reduced coordinates
		double xpos = coords[i][0]*MySystem.getG1()[0]+coords[i][1]*MySystem.getG2()[0]+coords[i][2]*MySystem.getG3()[0];
		double ypos = coords[i][0]*MySystem.getG1()[1]+coords[i][1]*MySystem.getG2()[1]+coords[i][2]*MySystem.getG3()[1];
		double zpos = coords[i][0]*MySystem.getG1()[2]+coords[i][1]*MySystem.getG2()[2]+coords[i][2]*MySystem.getG3()[2];
		// wrap
		if( xpos >= 1. || xpos < 0. ) coords[i][0] = xpos-floor(xpos);
		else coords[i][0] = xpos;
		if( ypos >= 1. || ypos < 0. ) coords[i][1] = ypos-floor(ypos);
		else coords[i][1] = ypos;
		if( zpos >= 1. || zpos < 0. ) coords[i][2] = zpos-floor(zpos);
		else coords[i][2] = zpos;
		// normal coordinates
		for(unsigned int k=0;k<3;k++){
			TempVec[k] = coords[i][0]*MySystem.getH1()[k] + coords[i][1]*MySystem.getH2()[k] + coords[i][2]*MySystem.getH3()[k];
			coords[i][k] = TempVec[k];
		}
	}
	delete[] TempVec;

}

vector<double> CenterGBSubsystems(AtomicSystem &MySystem, vector<vector<unsigned int>> &GBIons, vector<double*> &coords, unsigned int &g){
	
	cout << "Centering the system" << endl;

	double tolBound = 2.; // in A, for an ions to be consider in contact with box boundaries
	unsigned int max_count = 30; // number of shift of the entire system for centering the ions TODO warning here verify because in the old version with max_count=10 a lot of GB cannot be centered
	double tolRedX = 2./MySystem.getH1()[0]; 
	double tolRedY = 2./MySystem.getH2()[1]; 
	double tolRedZ = 2./MySystem.getH3()[2]; 
	bool XCentered = true;
	bool YCentered = true;
	bool ZCentered = true;
	vector<double> TransVec;

	for(unsigned int gbid=0;gbid<3;gbid++){
		for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
			coords[gbid][i*3] = MySystem.getWrappedPos(GBIons[g*3+gbid][i]).x;
			coords[gbid][i*3+1] = MySystem.getWrappedPos(GBIons[g*3+gbid][i]).y;
			coords[gbid][i*3+2] = MySystem.getWrappedPos(GBIons[g*3+gbid][i]).z;
		}
	}
	for(unsigned int gbid=0;gbid<3;gbid++){
		for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
			// Compute reduced coordinates
			double xpos = coords[gbid][i*3]*MySystem.getG1()[0]+coords[gbid][i*3+1]*MySystem.getG2()[0]+coords[gbid][i*3+2]*MySystem.getG3()[0];
			double ypos = coords[gbid][i*3]*MySystem.getG1()[1]+coords[gbid][i*3+1]*MySystem.getG2()[1]+coords[gbid][i*3+2]*MySystem.getG3()[1];
			double zpos = coords[gbid][i*3]*MySystem.getG1()[2]+coords[gbid][i*3+1]*MySystem.getG2()[2]+coords[gbid][i*3+2]*MySystem.getG3()[2];
			coords[gbid][i*3] = xpos;
			coords[gbid][i*3+1] = ypos;
			coords[gbid][i*3+2] = zpos;
			if( coords[gbid][i*3] < tolRedX || coords[gbid][i*3] > 1.-tolRedX ) XCentered = false;
			if( coords[gbid][i*3+1] < tolRedY || coords[gbid][i*3+1] > 1.-tolRedY ) YCentered = false;
			if( coords[gbid][i*3+2] < tolRedZ || coords[gbid][i*3+2] > 1.-tolRedZ ) ZCentered = false;
		}
	}
	double transvec_x = 0.;
	double transvec_y = 0.;
	double transvec_z = 0.;
	unsigned int count = 0;
	double *TempVec = new double[3];
	while( ( !XCentered || !YCentered || !ZCentered ) && ( count < max_count ) ){
		if( !XCentered ){
			XCentered = true;
			YCentered = true;
			ZCentered = true;
			transvec_x += 1./( (double) max_count );
			for(unsigned int gbid=0;gbid<3;gbid++){
				for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
					// Change reduced coordinates
					coords[gbid][i*3] += 1./( (double) max_count );
					// Wrapped reduced coordinates
					if( coords[gbid][i*3] >= 1. || coords[gbid][i*3] < 0. ) coords[gbid][i*3] = coords[gbid][i*3]-floor(coords[gbid][i*3]);
					// Wrapped normal coordinates
					for(unsigned int k=0;k<3;k++) TempVec[k] = coords[gbid][i*3]*MySystem.getH1()[k] + coords[gbid][i*3+1]*MySystem.getH2()[k] + coords[gbid][i*3+2]*MySystem.getH3()[k];
					// New reduced coordinates
					for(unsigned int k=0;k<3;k++) coords[gbid][i*3+k] = TempVec[0]*MySystem.getG1()[k] + TempVec[1]*MySystem.getG2()[k] + TempVec[2]*MySystem.getG3()[k];
					if( coords[gbid][i*3] < tolRedX || coords[gbid][i*3] > 1.-tolRedX ) XCentered = false;
					if( coords[gbid][i*3+1] < tolRedY || coords[gbid][i*3+1] > 1.-tolRedY ) YCentered = false;
					if( coords[gbid][i*3+2] < tolRedZ || coords[gbid][i*3+2] > 1.-tolRedZ ) ZCentered = false;
				}
			}
		}
		if( !YCentered ){
			XCentered = true;
			YCentered = true;
			ZCentered = true;
			transvec_y += 1./( (double) max_count );
			for(unsigned int gbid=0;gbid<3;gbid++){
				for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
					// Change reduced coordinates
					coords[gbid][i*3+1] += 1./( (double) max_count );
					// Wrapped reduced coordinates
					if( coords[gbid][i*3+1] >= 1. || coords[gbid][i*3+1] < 0. ) coords[gbid][i*3+1] = coords[gbid][i*3+1]-floor(coords[gbid][i*3+1]);
					// Wrapped normal coordinates
					for(unsigned int k=0;k<3;k++) TempVec[k] = coords[gbid][i*3]*MySystem.getH1()[k] + coords[gbid][i*3+1]*MySystem.getH2()[k] + coords[gbid][i*3+2]*MySystem.getH3()[k];
					// New reduced coordinates
					for(unsigned int k=0;k<3;k++) coords[gbid][i*3+k] = TempVec[0]*MySystem.getG1()[k] + TempVec[1]*MySystem.getG2()[k] + TempVec[2]*MySystem.getG3()[k];
					if( coords[gbid][i*3+0] < tolRedX || coords[gbid][i*3+0] > 1.-tolRedX ) XCentered = false;
					if( coords[gbid][i*3+1] < tolRedY || coords[gbid][i*3+1] > 1.-tolRedY ) YCentered = false;
					if( coords[gbid][i*3+2] < tolRedZ || coords[gbid][i*3+2] > 1.-tolRedZ ) ZCentered = false;
				}
			}
		}
		if( !ZCentered ){
			XCentered = true;
			YCentered = true;
			ZCentered = true;
			transvec_z +=1./( (double) max_count );
			for(unsigned int gbid=0;gbid<3;gbid++){
				for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
					// Change reduced coordinates
					coords[gbid][i*3+2] += 1./( (double) max_count );
					// Wrapped reduced coordinates
					if( coords[gbid][i*3+2] >= 1. || coords[gbid][i*3+2] < 0. ) coords[gbid][i*3+2] = coords[gbid][i*3+2]-floor(coords[gbid][i*3+2]);
					// Wrapped normal coordinates
					for(unsigned int k=0;k<3;k++) TempVec[k] = coords[gbid][i*3]*MySystem.getH1()[k] + coords[gbid][i*3+1]*MySystem.getH2()[k] + coords[gbid][i*3+2]*MySystem.getH3()[k];
					// New reduced coordinates
					for(unsigned int k=0;k<3;k++) coords[gbid][i*3+k] = TempVec[0]*MySystem.getG1()[k] + TempVec[1]*MySystem.getG2()[k] + TempVec[2]*MySystem.getG3()[k];
					if( coords[gbid][i*3+0] < tolRedX || coords[gbid][i*3+0] > 1.-tolRedX ) XCentered = false;
					if( coords[gbid][i*3+1] < tolRedY || coords[gbid][i*3+1] > 1.-tolRedY ) YCentered = false;
					if( coords[gbid][i*3+2] < tolRedZ || coords[gbid][i*3+2] > 1.-tolRedZ ) ZCentered = false;
				}
			}
		}
	count += 1;	
	} // end while
	if( count == max_count ){
		cout << "We didn't find translation" << endl;
		for(unsigned int i=0;i<3;i++) TransVec.push_back(0.);
	}else{
		// Wrapped normal coordinates
		for(unsigned int gbid=0;gbid<3;gbid++){
			for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
				for(unsigned int k=0;k<3;k++){
					TempVec[k] = coords[gbid][i*3]*MySystem.getH1()[k] + coords[gbid][i*3+1]*MySystem.getH2()[k] + coords[gbid][i*3+2]*MySystem.getH3()[k];
					coords[gbid][i*3+k] = TempVec[k];
				}
			}
		}
		for(unsigned int k=0;k<3;k++) TransVec.push_back(transvec_x*MySystem.getH1()[k] + transvec_y*MySystem.getH2()[k] + transvec_z*MySystem.getH3()[k]);
		//cout << "Translation find with vector : " << transvec_x << " " << transvec_y << " " << transvec_z << endl; 
	}
	delete[] TempVec;
	return TransVec;
}

unsigned int GetIndexNbClust(double GBId, vector<double> &GBId_ForNbClust_Init){
	unsigned int index=0;
	double tol=1e-8; 
	for(unsigned int i=0;i<GBId_ForNbClust_Init.size();i++){
		if( fabs(GBId-GBId_ForNbClust_Init[i]) < tol ){
			index=i;
			break;
		}
	}
	return index;
}

vector<unsigned int> GetIndexesPrevious(double &GBId, vector<double> &GBId_ForNbClust_Init, vector<vector<vector<unsigned int>>> &Previous_GBIons){
	vector<unsigned int> indexes;
	double tol=1e-8;
	for(unsigned int i=0;i<GBId_ForNbClust_Init.size();i++){
		if( fabs(GBId-GBId_ForNbClust_Init[i]) < tol && Previous_GBIons[i][0].size() != 0 && Previous_GBIons[i][1].size() != 0 ) indexes.push_back(i);
	}
	cout << "GBId " << GBId << endl;
	for(unsigned int i=0;i<indexes.size();i++) cout << indexes[i] << endl;
	return indexes;
}


void AverageStressAndStrain(AtomicSystem &MySystem, double *MeanStrain_1, double *MeanStrain_2, double *MeanStrain_am, long double *MeanStress_1, long double *MeanStress_2, long double *MeanStress_am, vector<vector<double>> &MeanAndCovarClust_Final, vector<vector<unsigned int>> &GBIons_Final, unsigned int &AtStrain_ind, unsigned int &size_AtStrain, unsigned int &Stress_ind, unsigned int &size_Stress, unsigned int &AtVol_ind, unsigned int &size_AtVol, unsigned int &i, double &am_thick){
	//cout << "GB " << i << " " << GBId_arr_Final[i] << endl;
	//cout << "N1 = " << MeanAndCovarClust_Final[i][0] << " " << MeanAndCovarClust_Final[i][1] << " " << MeanAndCovarClust_Final[i][2] << endl; 
	//cout << "N2 = " << MeanAndCovarClust_Final[i][6] << " " << MeanAndCovarClust_Final[i][7] << " " << MeanAndCovarClust_Final[i][8] << endl; 
	//cout << "Average normal = "  << MeanAndCovarClust_Final[i][12] << " " << MeanAndCovarClust_Final[i][13] << " " << MeanAndCovarClust_Final[i][14] << endl; 
	//cout << "M1 = " << MeanAndCovarClust_Final[i][3] << " " << MeanAndCovarClust_Final[i][4] << " " << MeanAndCovarClust_Final[i][5] << endl; 
	//cout << "M2 = " << MeanAndCovarClust_Final[i][9] << " " << MeanAndCovarClust_Final[i][10] << " " << MeanAndCovarClust_Final[i][11] << endl; 
	double d1 = 0.;
	double d2 = 0.;
	for(unsigned int dim=0;dim<3;dim++){
		//MeanAndCovarClust_Final[i][3+dim] -= TransVec[i*6+dim];
		//MeanAndCovarClust_Final[i][9+dim] -= TransVec[i*6+3+dim];
		d1 += MeanAndCovarClust_Final[i][3+dim]*MeanAndCovarClust_Final[i][12+dim];
		d2 += MeanAndCovarClust_Final[i][9+dim]*MeanAndCovarClust_Final[i][12+dim];
	}
	if( d1 > d2 ) am_thick = d1-d2;
	else am_thick = d2-d1;
	//cout << "Amorphous thickness : " << am_thick << endl;
	for(unsigned int ii=0;ii<size_AtStrain;ii++){
		MeanStrain_1[ii] = 0.;
		MeanStrain_2[ii] = 0.;
		MeanStrain_am[ii] = 0.;
	}
	for(unsigned int ii=0;ii<size_Stress;ii++){
		MeanStress_1[ii] = 0.;
		MeanStress_2[ii] = 0.;
		MeanStress_am[ii] = 0.;
	}
	// compute average strain and stress for interface 1
	double temp_vol = 0.;
	for(unsigned int n=0;n<GBIons_Final[i*3].size();n++){
		for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_1[ii] += MySystem.getAux(AtStrain_ind)[GBIons_Final[i*3][n]*size_AtStrain+ii];
		for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_1[ii] += (long double) MySystem.getAux(Stress_ind)[GBIons_Final[i*3][n]*size_Stress+ii];
	       	temp_vol += MySystem.getAux(AtVol_ind)[GBIons_Final[i*3][n]*size_AtVol];
	}
	for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_1[ii] /= temp_vol;
	temp_vol = 0.;
	for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_1[ii] /= GBIons_Final[i*3].size();
	// compute average strain and stress for interface 2
	for(unsigned int n=0;n<GBIons_Final[i*3+1].size();n++){
		for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_2[ii] += MySystem.getAux(AtStrain_ind)[GBIons_Final[i*3+1][n]*size_AtStrain+ii];
		for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_2[ii] += (long double) MySystem.getAux(Stress_ind)[GBIons_Final[i*3+1][n]*size_Stress+ii];
	        temp_vol += MySystem.getAux(AtVol_ind)[GBIons_Final[i*3+1][n]*size_AtVol];
	}
	for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_2[ii] /= temp_vol;
	temp_vol = 0.;
	for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_2[ii] /= GBIons_Final[i*3+1].size();
	// compute average strain and stress for amorph 
	for(unsigned int n=0;n<GBIons_Final[i*3+2].size();n++){
		for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_am[ii] += MySystem.getAux(AtStrain_ind)[GBIons_Final[i*3+2][n]*size_AtStrain+ii];
		for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_am[ii] += (long double) MySystem.getAux(Stress_ind)[GBIons_Final[i*3+2][n]*size_Stress+ii];
	        temp_vol += MySystem.getAux(AtVol_ind)[GBIons_Final[i*3+2][n]*size_AtVol];
	}
	for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_am[ii] /= temp_vol;
	for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_am[ii] /= GBIons_Final[i*3+2].size();
}

void printAtomicFile(AtomicSystem &MySystem, string &fullname, string &timestep, vector<vector<unsigned int>> &GBIons_Final, unsigned int &i){
	MathTools MT;
	std::ofstream f_out(fullname);
	f_out << "ITEM: TIMESTEP" << endl;
	f_out << timestep << endl;
	f_out << "ITEM: NUMBER OF ATOMS" << endl;
	f_out << GBIons_Final[i*3].size()+GBIons_Final[i*3+1].size()+GBIons_Final[i*3+2].size() << endl;
	double arr[4] = {0.,MySystem.getH2()[0],MySystem.getH3()[0],MySystem.getH2()[0]+MySystem.getH3()[0]};
        double arr_2[2] = {0.,MySystem.getH3()[1]};
	f_out << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n" << MT.min(arr,4) << "\t" << MySystem.getH1()[0]+MT.max(arr,4) << "\t" << MySystem.getH2()[0] << "\n" << MT.min(arr_2,2) << "\t" << MySystem.getH2()[1]+MT.max(arr_2,2) << "\t" << MySystem.getH3()[0] << "\n0\t" << MySystem.getH3()[2] << "\t" << MySystem.getH3()[1] << "\n";
	f_out << "ITEM: ATOMS id x y z clustId" << endl;
	for(unsigned int n=0;n<GBIons_Final[i*3].size();n++){
		f_out << n+1 << " " << MySystem.getWrappedPos(GBIons_Final[i*3][n]).x << " " << MySystem.getWrappedPos(GBIons_Final[i*3][n]).y << " " << MySystem.getWrappedPos(GBIons_Final[i*3][n]).z << " " << 0 << endl;
	}
	for(unsigned int n=0;n<GBIons_Final[i*3+1].size();n++){
		f_out << n+1+GBIons_Final[i*3].size() << " " << MySystem.getWrappedPos(GBIons_Final[i*3+1][n]).x << " " << MySystem.getWrappedPos(GBIons_Final[i*3+1][n]).y << " " << MySystem.getWrappedPos(GBIons_Final[i*3+1][n]).z << " " << 1 << endl;
	}
	for(unsigned int n=0;n<GBIons_Final[i*3+2].size();n++){
		f_out << n+1+GBIons_Final[i*3].size()+GBIons_Final[i*3+1].size() << " " << MySystem.getWrappedPos(GBIons_Final[i*3+2][n]).x << " " << MySystem.getWrappedPos(GBIons_Final[i*3+2][n]).y << " " << MySystem.getWrappedPos(GBIons_Final[i*3+2][n]).z << " " << 2 << endl;
	}
	f_out.close();
}

void IncreaseSmallestEigval(vector<vector<double>> &Matrix, double facIncrease){
	MathTools MT;
	double *eigval = new double[3];
	double *eigvec = new double[9];
	double *buffer_mat_1 = new double[9];
	double *buffer_mat_2 = new double[9];
	MT.EigenDecomposition(Matrix,eigval,eigvec);
	eigval[2] *= facIncrease;
	for(unsigned int i=0;i<9;i++) buffer_mat_1[i] = 0.;
	for(unsigned int i=0;i<3;i++){
		buffer_mat_1[i*3+i] = eigval[i];
		for(unsigned int j=0;j<3;j++) buffer_mat_2[i*3+j] = eigvec[j*3+i];
	}
	MT.MatDotMat(buffer_mat_2,buffer_mat_1,buffer_mat_1);
	MT.invert3x3(buffer_mat_2,buffer_mat_2);
	MT.MatDotMat(buffer_mat_1,buffer_mat_2,buffer_mat_1);
	for(unsigned int i=0;i<3;i++){
		for(unsigned int j=0;j<3;j++) Matrix[i][j] = buffer_mat_1[i*3+j];
	}
	delete[] eigval;
	delete[] eigvec;
	delete[] buffer_mat_1;
	delete[] buffer_mat_2;
}

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();

	// main steps of the algorithm
	// 1. for the first timestep, we provide the number of cluster and the number of GB to have a first precise evaluation of the GBs
	// 2. we store in Previous_variables for each GBId, the GBs TransVec (for wrapping at center of cell), GMM properties (means covars weights) and the atoms ids belonging to the GB
	// 3. For the other timesteps, we begin by searching the clusters with the previous number of cluster and search if we retrieve the GB stored in previous timestep
	// 	- if the number of ions in the found GB is lower than expected, we decrease the number of cluster
	// 	- if the number of ions is higher than expected, we increase the number of cluster

	// Here we use a dump file containing (1) the structure of a given ion (either interface, amorph or crystal), (2) the GB Id, (3) the atomic volume, (4) strain and (5) stress tensor
	// the objective is to return a file for each GB Id with:
	// nx ny nz AmorphousThickness GBSurface Interface1Strain(AndStress)Tensor Interface2Strain(AndStress)Tensor AmorphousStrain(AndStress)Tensor Interface1Strain(AndStress)Invariant Interface2Strain(AndStress)Invariants AmorphousStrain(AndStress)Invariants 
	string InputFilenamePrefix, OutputFilename, GTFilename, Filename_NbClustAndGB_Init, timesteps_list;
	vector<string> timesteps;

	if( argc == 6 ){
		InputFilenamePrefix = argv[1]; // Full_
		GTFilename = argv[2]; // GrainTags
		OutputFilename = argv[3];
		timesteps_list= argv[4]; //
		Filename_NbClustAndGB_Init = argv[5];
		//istringstream iss_rc(argv[4]);
		//iss_rc >> timestep;
		//istringstream iss_c(argv[4]);
		//iss_c >> tolspangle;
	}else{
		cerr << "Usage: ./AnalyzePolycrystal InputFilename GrainTagsFilename OutputFilename Timestep" << endl;
		cerr << "TODO description" << endl;
		return EXIT_FAILURE;
	}

	ifstream file(timesteps_list);
	if( file ){
		string line;
		while(getline(file,line)) timesteps.push_back(line);
	}else{
		cerr << "Cannot open file containing the timestep list" << endl;
		exit(EXIT_FAILURE);
	}
			
	// TODO make global variables
	double tolspangle = 10.;
	double tol_sp = cos(tolspangle*M_PI/180.);
	double facClust = 0.15;
	double facNbIons = .4; // increase (1.+fac) or decrease (1.-fac) factor for the number ions allowed for following a GB from one timestep to another (used to adapt number of cluster for GMM)
	long double minProbAmorph = 1e-6;
	unsigned int nbMinIonsInGB = 500;
	double facIncreaseSmallestEigval = 2.; // increasing factor of the smallest eigenvalue of average gaussian cluster for identification of amorphous ions 
	unsigned int Struct_GB = 5; //TODO add in argument of exe
	unsigned int Struct_Amorph = 0;
	unsigned int nclust_min_opt=1;
	unsigned int nclust_max_opt=4;
	
	unsigned int zero=0;
	MathTools MT;

	vector<double> GBIdOverTimestep; // GBIdOverTimestep[i] = GBId of the ith stored GB
	vector<vector<vector<double>>> Previous_MeanAndCovarClust; // Previous_MeanAndCovarClust[i][0][dim] = mean of interface 1 of the ith stored GB, [3+dim] = normal of interface 1 of the ith stored GB, [1][dim] = mean of interface 2 of the ith stored GB ..
	vector<vector<vector<unsigned int>>> Previous_GBIons; // Previous_GBIons[i][0][n] = id of the nth ion belonging to interface 1 of the ith GB stored, [i][1][n] = same but in interface 2
	vector<vector<double>> Previous_TransVec; // Previous_TransVec[i][dim] = TransVec of ith stored GB
	vector<vector<unsigned int>> Previous_NbClust; // Previous_NbClust[i][0] = number of cluster used for GMM for the ith GB in interface 1, [i][1] => interface 2

	// ### BEGIN FIRST TREATING FOR INITIAL TIMESTEP ### //
	bool init = true;
	vector<vector<unsigned int>> NbClustAndGB_Init;
	vector<double> GBId_ForNbClust_Init;
	if( init ){
		// Read initial nb cluster and GB
		ReadFile_NbClustAndGB_Init(Filename_NbClustAndGB_Init, NbClustAndGB_Init, GBId_ForNbClust_Init);
		
		AtomicSystem MySystem(InputFilenamePrefix+timesteps[0]+".cfg");
		
		// Search the cell vector of each grains
		ComputeAndPrintCellVectors(MySystem,GTFilename,timesteps[0]);
		
		// Get the different aux properties
		unsigned int size_Struct, Struct_ind, size_GBId, GBId_ind, size_AtVol, AtVol_ind, size_Stress, Stress_ind, size_AtStrain, AtStrain_ind;
		getAuxIdAndSize(MySystem, size_Struct, Struct_ind, size_GBId, GBId_ind, size_AtVol, AtVol_ind, size_Stress, Stress_ind, size_AtStrain, AtStrain_ind);

		// Search which ions belongs to GB and differentiate the one in amorph and the one in interface
		vector<double> GBId_arr;
		vector<vector<unsigned int>> GBIons;
		ConstructGBIdAndGBIonsArrays(MySystem, GBId_arr, GBIons, size_Struct, Struct_ind, size_GBId, GBId_ind, Struct_GB, Struct_Amorph);

		vector<vector<double>> MeanAndCovarClust_temp;
		vector<vector<unsigned int>> GBIons_temp;
		vector<vector<double>> MeanAndCovarClust_Final;
		vector<double> GBId_arr_Final;
		vector<vector<unsigned int>> GBIons_Final;
		vector<vector<double>> TransVec;
		unsigned int nbGBAnalyzed = 0;
		unsigned int current_nbGBAnalyzed;
		for(unsigned int g=0;g<GBId_arr.size();g++){
			cout << "Treating " << GBId_arr[g] << " GB" << endl;
			if( GBIons[g*3].size() < nbMinIonsInGB || GBIons[g*3+1].size() < nbMinIonsInGB ){
				cout << "skipping, not enough ions" << endl;
				continue;
			}

			unsigned int index_NbClust = GetIndexNbClust(GBId_arr[g],GBId_ForNbClust_Init);

			MeanAndCovarClust_temp.clear();
			MeanAndCovarClust_temp.push_back(vector<double>());
			MeanAndCovarClust_temp.push_back(vector<double>());

			// Variables for identifying amorph ions 
			vector<vector<vector<double>>> AverageCov; // compute before the avergae of covariance and then inverse for identification of amorph
			vector<vector<vector<double>>> AverageInvCov;
			vector<vector<double>> AverageMu;
			vector<vector<vector<double>>> AverageCov_temp;
			vector<vector<double>> AverageMu_temp;
			vector<double> AverageDetC;

			for(unsigned int i=0;i<GBIons_temp.size();i++){
				GBIons_temp[i].clear();
				vector<unsigned int>().swap(GBIons_temp[i]);
			}
			GBIons_temp.clear();
			vector<vector<unsigned int>>().swap(GBIons_temp);
			current_nbGBAnalyzed = 0;

			// Center the system and store the vector needed for it
			vector<double*> current_coords;
			for(unsigned int gbid=0;gbid<3;gbid++) current_coords.push_back(new double[GBIons[g*3+gbid].size()*3]);
			TransVec.push_back(CenterGBSubsystems(MySystem,GBIons,current_coords,g));
								
			// Shift and center the amorph ions with the same vector
			vector<vector<double>> current_coords_am;
			for(unsigned int id=0;id<GBIons[g*3+2].size();id++){
				current_coords_am.push_back(vector<double>());
				for(unsigned int kdim=0;kdim<3;kdim++) current_coords_am[id].push_back(current_coords[2][id*3+kdim]);
			}
			//ShiftAndWrapAmorph(MySystem,GBIons,current_coords_am,g,TransVec[g]);

			for(unsigned int gbid=0;gbid<2;gbid++){
				// find the GMM
				Descriptors GMM_des(current_coords[gbid],GBIons[g*3+gbid].size(),3);
				GaussianMixtureModel GMM;
				GMM.setDescriptors(&GMM_des);
				unsigned int zero=0;
				GMM.TrainModel(NbClustAndGB_Init[index_NbClust][gbid],zero);
				unsigned int nbClustFinal = GMM.getNbClust(zero);
				GMM.Classify();
				vector<unsigned int> *Ids = new vector<unsigned int>[nbClustFinal];
				for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
					if( GMM.getClassificator()[i*2+1] != 0 ) Ids[(unsigned int) GMM.getClassificator()[i*2]].push_back(GBIons[g*3+gbid][i]);
				}
				
				for(unsigned int i=0;i<nbClustFinal;i++){
					if( Ids[i].size() >  nbMinIonsInGB ){
						double *covar = new double[9];
						double *eigval = new double[3];
						double *eigvec = new double[9];
						for(unsigned int d1=0;d1<3;d1++){
							for(unsigned int d2=0;d2<3;d2++) covar[d1*3+d2] = GMM.getCov()[i*9+d1*3+d2];
						}
						MT.EigenDecomposition(covar,3,eigval,eigvec);
						bool tostore = false;
						if( ( eigval[0] < facClust*eigval[1] ) && ( eigval[0] < facClust*eigval[2] ) ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp[gbid].push_back(eigvec[dim]);	
							tostore = true;
						}else if( ( eigval[1] < facClust*eigval[0] ) && ( eigval[1] < facClust*eigval[2] ) ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp[gbid].push_back(eigvec[3+dim]);	
							tostore = true;
						}else if( ( eigval[2] < facClust*eigval[1] ) && ( eigval[2] < facClust*eigval[0] ) ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp[gbid].push_back(eigvec[6+dim]);	
							tostore = true;
						} 
						if( tostore && gbid == 0 ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp[gbid].push_back(GMM.getMu()[i*3+dim]);	
							GBIons_temp.push_back(vector<unsigned int>());
							for(unsigned int n=0;n<Ids[i].size();n++) GBIons_temp[GBIons_temp.size()-1].push_back(Ids[i][n]);
							// Store cov and mean for averaging and identifying amorph ions 
							AverageCov_temp.push_back(vector<vector<double>>(3,vector<double>(3)));
							AverageMu_temp.push_back(vector<double>(3));
							for(unsigned int am_g_1=0;am_g_1<3;am_g_1++){
								AverageMu_temp[AverageMu_temp.size()-1][am_g_1] = GMM.getMu()[i*3+am_g_1];
								for(unsigned int am_g_2=0;am_g_2<3;am_g_2++) AverageCov_temp[AverageCov_temp.size()-1][am_g_1][am_g_2] = GMM.getCov()[i*9+am_g_1*3+am_g_2];
							}
						}else if( tostore && gbid == 1 ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp[gbid].push_back(GMM.getMu()[i*3+dim]);	
							// search if there is corresponding GB 1 with almost the same normal
							bool isAssociated = false;
							double n0, n1, sp, fullsp;
							vector<unsigned int> indexes;
							vector<double> scalarprod;
							unsigned int index_clust_ass;
							for(unsigned int n=0;n<GBIons_temp.size();n++){
								sp = 0.;
								n0 = 0.;
								n1 = 0.;
								for(unsigned int dim=0;dim<3;dim++){
									sp += MeanAndCovarClust_temp[0][6*n+dim]*MeanAndCovarClust_temp[1][((MeanAndCovarClust_temp[1].size()/6)-1)*6+dim];
									n0 += pow(MeanAndCovarClust_temp[0][6*n+dim],2.);
									n1 += pow(MeanAndCovarClust_temp[1][((MeanAndCovarClust_temp[1].size()/6)-1)*6+dim],2.);
								}
								n0 = sqrt(n0);
								n1 = sqrt(n1);
								fullsp = fabs(sp/(n1*n0));
								if( fullsp > tol_sp ){
									isAssociated = true;
									scalarprod.push_back(fullsp);
									indexes.push_back(n);
								}
							}
							if( isAssociated ){

								index_clust_ass = indexes[MT.max(scalarprod)];
								nbGBAnalyzed += 1;
								current_nbGBAnalyzed += 1;
								GBIons_Final.push_back(vector<unsigned int>());
								GBIons_Final.push_back(vector<unsigned int>());
								GBIons_Final.push_back(vector<unsigned int>());

								// Average covariance and mean of the two interfaces and search amorph ions belonging to this GB 
								AverageCov.push_back(vector<vector<double>>(3,vector<double>(3)));
								AverageMu.push_back(vector<double>(3));
								for(unsigned int am_g_1=0;am_g_1<3;am_g_1++){
									AverageMu[current_nbGBAnalyzed-1][am_g_1] = GMM.getMu()[i*3+am_g_1];
									AverageMu[current_nbGBAnalyzed-1][am_g_1] += AverageMu_temp[index_clust_ass][am_g_1];
									AverageMu[current_nbGBAnalyzed-1][am_g_1] /= 2.;
									for(unsigned int am_g_2=0;am_g_2<3;am_g_2++){
										AverageCov[current_nbGBAnalyzed-1][am_g_1][am_g_2] = GMM.getCov()[i*9+am_g_1*3+am_g_2];
										AverageCov[current_nbGBAnalyzed-1][am_g_1][am_g_2] += AverageCov_temp[index_clust_ass][am_g_1][am_g_2];
										AverageCov[current_nbGBAnalyzed-1][am_g_1][am_g_2] /= 2.;
									}
								}
								AverageInvCov.push_back(vector<vector<double>>(3,vector<double>(3)));
								IncreaseSmallestEigval(AverageCov[current_nbGBAnalyzed-1], facIncreaseSmallestEigval);
								AverageDetC.push_back(MT.invert3x3(AverageCov[current_nbGBAnalyzed-1],AverageInvCov[current_nbGBAnalyzed-1]));
								// an amorph ion is considered to belong to this GB if its prob is higher than minProbAmorph
								for(unsigned int am=0;am<GBIons[g*3+2].size();am++){
									if( MT.Prob_MultidimGaussian(AverageInvCov[current_nbGBAnalyzed-1],AverageMu[current_nbGBAnalyzed-1],AverageDetC[current_nbGBAnalyzed-1],current_coords_am[am]) > minProbAmorph ) GBIons_Final[(nbGBAnalyzed-1)*3+2].push_back(GBIons[g*3+2][am]);
								}

								// Store normals and mean positions for this GB and average of normals
								MeanAndCovarClust_Final.push_back(vector<double>());
								for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back(MeanAndCovarClust_temp[0][6*index_clust_ass+dim]);
								for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back(MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]);
								double x1 = MeanAndCovarClust_temp[0][6*index_clust_ass+0];
								double y1 = MeanAndCovarClust_temp[0][6*index_clust_ass+1];
								double z1 = MeanAndCovarClust_temp[0][6*index_clust_ass+2];
								double x2 = MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+0];
								double y2 = MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+1];
								double z2 = MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+2];
								if( ( fabs(x1) > 1e-1 && fabs(x2) > 1e-1 && x1*x2 < 0. ) ||  ( fabs(y1) > 1e-1 && fabs(y2) > 1e-1 && y1*y2 < 0. ) || ( fabs(z1) > 1e-1 && fabs(z2) > 1e-1 && z1*z2 < 0. ) ){
									double norm = 0.;
									for(unsigned int dim=0;dim<3;dim++) norm += pow((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]-MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/2.,2.);
									norm = sqrt(norm);
									for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]-MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/(2.*norm));
								}else{
									double norm = 0.;
									for(unsigned int dim=0;dim<3;dim++) norm += pow((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]+MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/2.,2.);
									for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]+MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/(2.*norm));
								}
								// Store ions ID for this GB
								for(unsigned int nb=0;nb<GBIons_temp[index_clust_ass].size();nb++) GBIons_Final[(nbGBAnalyzed-1)*3].push_back(GBIons_temp[index_clust_ass][nb]);
								for(unsigned int nb=0;nb<Ids[i].size();nb++) GBIons_Final[(nbGBAnalyzed-1)*3+1].push_back(Ids[i][nb]);

								// Stored ids of grain
								GBId_arr_Final.push_back(GBId_arr[g]);

								// delete this cluster
								for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_temp[0].erase(MeanAndCovarClust_temp[0].begin()+6*index_clust_ass);
								GBIons_temp.erase(GBIons_temp.begin()+index_clust_ass);
							}
						}
						delete[] covar;
						delete[] eigval;
						delete[] eigvec;
					}
				} // end nclust loop
				for(unsigned int id=0;id<nbClustFinal;id++){
					Ids[id].clear();
					vector<unsigned int>().swap(Ids[id]);
				}
				delete[] Ids;
			} // end gbid loop
			// free memory of coordinates
			for(unsigned int gbid=0;gbid<2;gbid++) delete[] current_coords[gbid];
			vector<double*>().swap(current_coords);
			cout << current_nbGBAnalyzed << " interfaces analyzed for this GB" << endl;
		}

		for(unsigned int i=0;i<nbGBAnalyzed;i++){
			GBIdOverTimestep.push_back(GBId_arr_Final[i]);
			Previous_GBIons.push_back(vector<vector<unsigned int>>());
			Previous_NbClust.push_back(vector<unsigned int>());
			unsigned int index_NbClust = GetIndexNbClust(GBId_arr_Final[i],GBId_ForNbClust_Init);
			for(unsigned int gbid=0;gbid<2;gbid++){
				Previous_NbClust[i].push_back(NbClustAndGB_Init[index_NbClust][gbid]);
				Previous_GBIons[i].push_back(vector<unsigned int>());
				for(unsigned int n=0;n<GBIons_Final[i*3+gbid].size();n++) Previous_GBIons[i][gbid].push_back(GBIons_Final[i*3+gbid][n]);
			}
		}
	
		double *MeanStrain_1 = new double[8];
		double *MeanStrain_2 = new double[8];
		double *MeanStrain_am = new double[8];
		long double *MeanStress_1 = new long double[6];
		long double *MeanStress_2 = new long double[6];
		long double *MeanStress_am = new long double[6];
		double am_thick;

		for(unsigned int i=0;i<nbGBAnalyzed;i++){

			AverageStressAndStrain(MySystem, MeanStrain_1, MeanStrain_2, MeanStrain_am, MeanStress_1, MeanStress_2, MeanStress_am, MeanAndCovarClust_Final, GBIons_Final, AtStrain_ind, size_AtStrain, Stress_ind, size_Stress, AtVol_ind, size_AtVol, i, am_thick);

			string path_1 = "./GB_";
			string path_2 = "/";
			string fullpath = path_1+to_string(i)+path_2;
			std::filesystem::create_directory(fullpath);

			std::ofstream f_out_res(fullpath+OutputFilename);
			f_out_res << "timestep(0) nx(1) ny(2) nz(3) G1(4) G2(5) AmThick(6) Eps_xx_1(7) Eps_yy_1(8) Eps_zz_1(9) Eps_xy_1(10) Eps_xz_1(11) Eps_yz_1(12) ShearInv_1(13) SphInv_1(14) Eps_xx_2(15) Eps_yy_2(16) Eps_zz_2(17) Eps_xy_2(18) Eps_xz_2(19) Eps_yz_2(20) ShearInv_2(21) SphInv_2(22) Eps_xx_am(23) Eps_yy_am(24) Eps_zz_am(25) Eps_xy_am(26) Eps_xz_am(27) Eps_yz_am(28) ShearInv_am(29) SphInv_am(30) sigma_xx_1(31) sigma_yy_1(32) sigma_zz_1(33) sigma_xy_1(34) sigma_xz_1(35) sigma_yz_1(36) sigma_xx_2(37) sigma_yy_2(38) sigma_zz_2(39) sigma_xy_2(40) sigma_xz_2(41) sigma_yz_2(42) sigma_xx_am(43) sigma_yy_am(44) sigma_zz_am(45) sigma_xy_am(46) sigma_xz_am(47) sigma_yz_am(48) nbAt(49)" << endl;

			f_out_res << timesteps[0] << " " << MeanAndCovarClust_Final[i][12] << " " << MeanAndCovarClust_Final[i][13] << " " << MeanAndCovarClust_Final[i][14] << " " << GBId_arr_Final[i] - (GBId_arr_Final[i]-floor(GBId_arr_Final[i])) << " " << (int) 10.*(GBId_arr_Final[i]-floor(GBId_arr_Final[i])) << " " << am_thick ;
			for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_1[ii]; 
			for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_2[ii]; 
			for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_am[ii]; 
			for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_1[ii]; 
			for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_2[ii]; 
			for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_am[ii];
		        f_out_res << " " << GBIons_Final[3*i].size() + GBIons_Final[3*i+1].size() + GBIons_Final[3*i+2].size();
			f_out_res << endl;
			f_out_res.close();

			string name = "GB_";
			string ext = ".cfg";
			string fullname = fullpath+name+timesteps[0]+ext;
			printAtomicFile(MySystem, fullname, timesteps[0], GBIons_Final, i);
		}
		
		delete[] MeanStrain_1;
		delete[] MeanStrain_2;
		delete[] MeanStrain_am;
		delete[] MeanStress_1;
		delete[] MeanStress_2;
		delete[] MeanStress_am;
	}	
	// ### END FIRST TREATING FOR INITIAL TIMESTEP ### //
	//
	// ### BEGIN TREAT ALL TIMESTEPS
	for(unsigned int t=1;t<timesteps.size();t++){
		// read system
		AtomicSystem MySystem(InputFilenamePrefix+timesteps[t]+".cfg");
		
		// Search the cell vector of each grains
		ComputeAndPrintCellVectors(MySystem,GTFilename,timesteps[t]);
		
		// Get the different aux properties
		unsigned int size_Struct, Struct_ind, size_GBId, GBId_ind, size_AtVol, AtVol_ind, size_Stress, Stress_ind, size_AtStrain, AtStrain_ind;
		getAuxIdAndSize(MySystem, size_Struct, Struct_ind, size_GBId, GBId_ind, size_AtVol, AtVol_ind, size_Stress, Stress_ind, size_AtStrain, AtStrain_ind);

		// Search which ions belongs to GB and differentiate the one in amorph and the one in interface
		vector<double> GBId_arr;
		vector<vector<unsigned int>> GBIons;
		ConstructGBIdAndGBIonsArrays(MySystem, GBId_arr, GBIons, size_Struct, Struct_ind, size_GBId, GBId_ind, Struct_GB, Struct_Amorph);

		vector<vector<double>> MeanAndCovarClust_temp;
		vector<vector<unsigned int>> GBIons_temp;
		vector<vector<double>> MeanAndCovarClust_Final;
		vector<vector<double>> MeanAndCovarClust_temp_prev;
		vector<vector<unsigned int>> GBIons_temp_prev;
		vector<vector<double>> MeanAndCovarClust_Final_prev;
		vector<double> GBId_arr_Final;
		vector<vector<unsigned int>> GBIons_Final;
		vector<vector<unsigned int>> GBIons_Final_prev;
		vector<vector<double>> TransVec;
		vector<unsigned int> IndexesInPreviousVariables;
		vector<vector<unsigned int>> NbClustToSaveInPrevious;
		unsigned int nbGBAnalyzed = 0;
		unsigned int nbGBAnalyzed_prev = 0;
		unsigned int current_nbGBAnalyzed;
		unsigned int current_nbGBAnalyzed_prev;
		for(unsigned int g=0;g<GBId_arr.size();g++){
			cout << "Treating " << GBId_arr[g] << " GB" << endl;
			if( GBIons[g*3].size() < nbMinIonsInGB || GBIons[g*3+1].size() < nbMinIonsInGB ){
				cout << "skipping, not enough ions" << endl;
				continue;
			}

			cout << 1 << endl;
		cout << "Size of prevs : " << GBIdOverTimestep.size() << " " << Previous_GBIons.size() << endl;
			vector<unsigned int> indexes_Previous = GetIndexesPrevious(GBId_arr[g],GBIdOverTimestep, Previous_GBIons);
			vector<vector<double>> *Covars_Identified;
			vector<double> *Means_Identified;
			bool ClusterIdentifiedToDelete = false;
			unsigned int ClusterIdentifiedSize = indexes_Previous.size();
			bool *fitting_ok_all;
			if( ClusterIdentifiedSize != 0 ){
				Covars_Identified = new vector<vector<double>>[indexes_Previous.size()];
				Means_Identified = new vector<double>[indexes_Previous.size()];
				fitting_ok_all = new bool[2*indexes_Previous.size()];
				ClusterIdentifiedToDelete = true;
			}
		                                               	
			vector<unsigned int> NbClustToSaveInPrevious_temp;
			
			MeanAndCovarClust_temp.clear();
			MeanAndCovarClust_temp.push_back(vector<double>());
			MeanAndCovarClust_temp.push_back(vector<double>());
			for(unsigned int i=0;i<GBIons_temp.size();i++){
				GBIons_temp[i].clear();
				vector<unsigned int>().swap(GBIons_temp[i]);
			}
			GBIons_temp.clear();
			vector<vector<unsigned int>>().swap(GBIons_temp);

			MeanAndCovarClust_temp_prev.clear();
			MeanAndCovarClust_temp_prev.push_back(vector<double>());
			MeanAndCovarClust_temp_prev.push_back(vector<double>());
			for(unsigned int i=0;i<GBIons_temp_prev.size();i++){
				GBIons_temp_prev[i].clear();
				vector<unsigned int>().swap(GBIons_temp_prev[i]);
			}
			GBIons_temp_prev.clear();
			vector<vector<unsigned int>>().swap(GBIons_temp_prev);

			// Variables for identifying amorph ions 
			vector<vector<vector<double>>> AverageCov; // compute before the avergae of covariance and then inverse for identification of amorph
			vector<vector<vector<double>>> AverageInvCov;
			vector<vector<double>> AverageMu;
			vector<vector<vector<double>>> AverageCov_temp;
			vector<vector<double>> AverageMu_temp;
			vector<double> AverageDetC;

			vector<vector<vector<double>>> AverageCov_prev; // compute before the avergae of covariance and then inverse for identification of amorph
			vector<vector<vector<double>>> AverageInvCov_prev;
			vector<vector<double>> AverageMu_prev;
			vector<vector<vector<double>>> AverageCov_temp_prev;
			vector<vector<double>> AverageMu_temp_prev;
			vector<double> AverageDetC_prev;

			current_nbGBAnalyzed = 0;
			current_nbGBAnalyzed_prev = 0;

			cout << 2 << endl;
			// Center the system and store the vector needed for it
			vector<double*> current_coords;
			for(unsigned int gbid=0;gbid<3;gbid++) current_coords.push_back(new double[GBIons[g*3+gbid].size()*3]);
			TransVec.push_back(CenterGBSubsystems(MySystem,GBIons,current_coords,g));
			cout << 3 << endl;
			// Compute the difference with the previous TransVec
			// For the moment we will test to identify the GBs bsaed onyl on ion indexes
			//vector<double> TransVec_diff(3);
			//for(unsigned int dim=0;dim<3;dim++) TransVec_diff[dim] = TransVec[dim] - Previous_TransVec[indexes_Previous[0]][dim]; // use only indexes_Previous[0] because its all the same for a given GBId
			// Shift and center the amorph ions with the same vector
			vector<vector<double>> current_coords_am;
			for(unsigned int id=0;id<GBIons[g*3+2].size();id++){
				current_coords_am.push_back(vector<double>());
				for(unsigned int kdim=0;kdim<3;kdim++) current_coords_am[id].push_back(current_coords[2][id*3+kdim]);
			}
			//ShiftAndWrapAmorph(MySystem,GBIons,current_coords_am,g,TransVec[g]);

			cout << 4 << endl;


			for(unsigned int gbid=0;gbid<2;gbid++){
				//vector<unsigned int> indexes_Previous_temp = GetIndexesPrevious(GBId_arr[g],GBIdOverTimestep, Previous_GBIons); // get one more time this array because it may have change depending in gbid value
				//indexes_Previous.clear();
				//vector<unsigned int>().swap(indexes_Previous);
				//for(unsigned int ind=0;ind<indexes_Previous_temp.size();ind++) indexes_Previous.push_back(indexes_Previous_temp[ind]);
				GaussianMixtureModel GMM; 
				GaussianMixtureModel Adjusted_GMM;
				Descriptors GMM_des(current_coords[gbid],GBIons[g*3+gbid].size(),3);
				bool fitting_ok = false;
				unsigned int max_count = 10;
				unsigned int count = 0;
				unsigned int current_NbClust;
				unsigned int *ClusterIndexAssociated;
				vector<vector<unsigned int>> Final_Ids;
				vector<unsigned int> *Prev_Final_Ids;
				if( indexes_Previous.size() != 0 ){
					current_NbClust = Previous_NbClust[indexes_Previous[0]][gbid];
					ClusterIndexAssociated = new unsigned int[indexes_Previous.size()];
					Prev_Final_Ids = new vector<unsigned int>[indexes_Previous.size()];
				}
				if( indexes_Previous.size() == 0 ){
					Adjusted_GMM.setDescriptors(&GMM_des);
					cout << "Trying to find new GBs using optimal fitting of GMM" << endl;
					Adjusted_GMM.fitOptimalGMM(nclust_min_opt,nclust_max_opt);
					current_NbClust = Adjusted_GMM.getNbClust(zero);
					Adjusted_GMM.Classify();
					for(unsigned int c=0;c<current_NbClust;c++) Final_Ids.push_back(vector<unsigned int>());
					for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
						if( Adjusted_GMM.getClassificator()[i*2+1] != 0 ) Final_Ids[(unsigned int) Adjusted_GMM.getClassificator()[i*2]].push_back(GBIons[g*3+gbid][i]);
					}
				}else{
					GMM.setDescriptors(&GMM_des);
					double *covar = new double[9];
					double *eigval = new double[3];
					double *eigvec = new double[9];
					cout << "Trying to retrieve previous GBs" << endl;
					for(unsigned int Np=0;Np<indexes_Previous.size();Np++) fitting_ok_all[Np*2+gbid] = false;
					while( !fitting_ok && count < max_count ){
						// find the GMM
						unsigned int zero=0;
						cout << "NbClust= " << current_NbClust << endl;
						GMM.TrainModel(current_NbClust,zero);
						GMM.Classify();
						vector<unsigned int> *Ids = new vector<unsigned int>[current_NbClust];
						unsigned int init_Ids_size = current_NbClust;
						for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
							if( GMM.getClassificator()[i*2+1] != 0 ) Ids[(unsigned int) GMM.getClassificator()[i*2]].push_back(GBIons[g*3+gbid][i]);
						}
						// loop over the GB already find and still active (with non-null number of ions) for this id
						// search wich cluster have most comon ions with the previous GB
						bool breaked = false;
						for(unsigned int Np=0;Np<indexes_Previous.size();Np++){
							if( !fitting_ok_all[Np*2+gbid] ){
								vector<double> RatioComonIon;
								for(unsigned int Nc=0;Nc<current_NbClust;Nc++){
									RatioComonIon.push_back(0.);
									for(unsigned int indp=0;indp<Previous_GBIons[indexes_Previous[Np]][gbid].size();indp++){
										for(unsigned int indc=0;indc<Ids[Nc].size();indc++){
											if( Previous_GBIons[indexes_Previous[Np]][gbid][indp] == Ids[Nc][indc] ){
												RatioComonIon[Nc] += 1;
												break;
											}
										}
									}
									RatioComonIon[Nc] /= Previous_GBIons[indexes_Previous[Np]][gbid].size();
								}
								ClusterIndexAssociated[Np] = MT.max(RatioComonIon);
								
								unsigned int old_nbions = Previous_GBIons[indexes_Previous[Np]][gbid].size();
								unsigned int new_nbions = Ids[ClusterIndexAssociated[Np]].size();
								// if the number of ion is much lower, decrease current_nbClust
								if( new_nbions < (1.-facNbIons)*old_nbions ){
									cout << "Decreasing nb clust, old ions " << old_nbions << ", new " << new_nbions << endl;
									if( current_NbClust > indexes_Previous.size() ){
										current_NbClust--;
										breaked = true;
										break;
									}else count = max_count;
								// if the number of ion is much higher, increase current_nbClust
								}else if( new_nbions > (1.+facNbIons)*old_nbions ){
									cout << "Increasing nb clust, old ions " << old_nbions << ", new " << new_nbions << endl;
									if( current_NbClust < GMM.getNbMaxCluster() ){
										current_NbClust++;
										breaked = true;
										break;
									}else count = max_count;
								}else{
									// search for facClust
									for(unsigned int d1=0;d1<3;d1++){
										for(unsigned int d2=0;d2<3;d2++) covar[d1*3+d2] = GMM.getCov()[ClusterIndexAssociated[Np]*9+d1*3+d2];
									}
									MT.EigenDecomposition(covar,3,eigval,eigvec);
									//if( ( ( eigval[0] < facClust*eigval[1] ) && ( eigval[0] < facClust*eigval[2] ) ) || ( ( eigval[1] < facClust*eigval[0] ) && ( eigval[1] < facClust*eigval[2] ) ) || ( ( eigval[2] < facClust*eigval[1] ) && ( eigval[2] < facClust*eigval[0] ) ) ){
									if( ( ( eigval[0] < facClust*eigval[1] ) && ( eigval[0] < eigval[2] ) ) || ( ( eigval[1] < facClust*eigval[0] ) && ( eigval[1] < eigval[2] ) ) || ( ( eigval[2] < facClust*eigval[1] ) && ( eigval[2] < eigval[0] ) ) || ( ( eigval[0] < facClust*eigval[2] ) && ( eigval[0] < eigval[1] ) ) || ( ( eigval[1] < facClust*eigval[2] ) && ( eigval[1] < eigval[0] ) ) || ( ( eigval[2] < facClust*eigval[0] ) && ( eigval[2] < eigval[1] ) ) ){
										fitting_ok_all[Np*2+gbid] = true;
										// save cluster prop
										for(unsigned int d1=0;d1<3;d1++){
											Means_Identified[Np].push_back(GMM.getMu()[ClusterIndexAssociated[Np]*3+d1]);
											Covars_Identified[Np].push_back(vector<double>(3));
											for(unsigned int d2=0;d2<3;d2++) Covars_Identified[Np][d1][d2] = covar[d1*3+d2];
										}
										// save ids of ions
										for(unsigned int id=0;id<Ids[ClusterIndexAssociated[Np]].size();id++) Prev_Final_Ids[Np].push_back(Ids[ClusterIndexAssociated[Np]][id]);

									}else{
										cout << "Increasing nb clust to have more planar GB" << endl;
										if( current_NbClust < GMM.getNbMaxCluster() ){
											current_NbClust++;
											breaked = true;
											break;
										}else count = max_count;
									}
									if( breaked ) break;
								}
								if( breaked ) break;
							}
						}
						fitting_ok = true;
						for(unsigned int Np=0;Np<indexes_Previous.size();Np++) fitting_ok *= fitting_ok_all[Np*2+gbid];
						if( fitting_ok ){
							cout << "All GBs were retrieved" << endl;
							//for(unsigned int c=0;c<current_NbClust;c++){
							//	Final_Ids.push_back(vector<unsigned int>());
							//	for(unsigned int id=0;id<Ids[c].size();id++) Final_Ids[c].push_back(Ids[c][id]);
							//}
						}
						cout << "clearing vars" << endl;
						for(unsigned int id=0;id<init_Ids_size;id++){
							Ids[id].clear();
							vector<unsigned int>().swap(Ids[id]);
						}
						delete[] Ids;
						count++;
					} // end while
					if( !fitting_ok ){
						cout << "Issue when identifying GB from previous timesteps" << endl;
						// Update Previous_GBIons to zero meaning that this GB are no longer active
						vector<unsigned int> index_identified;
						for(unsigned int ind=0;ind<indexes_Previous.size();ind++){
							if( !fitting_ok_all[ind*2+gbid] ){
								Previous_GBIons[indexes_Previous[ind]][gbid].clear(); // set to zero ions
								vector<unsigned int>().swap(Previous_GBIons[indexes_Previous[ind]][gbid]);
							}else index_identified.push_back(ind);
						}
						// then fit an optimal GMM
						cout << "Trying to find new GBs using optimal fitting of GMM" << endl;
						// remove descriptors that have been already found
						vector<double> adjusted_coords_temp;
						for(unsigned int at1=0;at1<GBIons[g*3+gbid].size();at1++){
							bool to_store = true;
							for(unsigned int ind=0;ind<index_identified.size();ind++){
								for(unsigned int ind_inv=0;ind_inv<Prev_Final_Ids[index_identified[ind]].size();ind_inv++){
									if( GBIons[g*3+gbid][at1] == Prev_Final_Ids[index_identified[ind]][ind_inv] ){
										to_store = false;
										break;
									}
								}
								if( !to_store ) break;
							}
							if( to_store ){
								for(unsigned int kdim=0;kdim<3;kdim++) adjusted_coords_temp.push_back(current_coords[gbid][at1*3+kdim]);
							}
						}
						unsigned int new_size = adjusted_coords_temp.size()/3;
						double *adjusted_coords = new double[new_size*3];
						for(unsigned int at1=0;at1<new_size;at1++){
							for(unsigned int kdim=0;kdim<3;kdim++) adjusted_coords[at1*3+kdim] = adjusted_coords_temp[at1*3+kdim];
						}
						Descriptors Adjusted_Des(adjusted_coords,new_size,3);
						Adjusted_GMM.setDescriptors(&Adjusted_Des);
						Adjusted_GMM.fitOptimalGMM(nclust_min_opt,nclust_max_opt);
						current_NbClust = Adjusted_GMM.getNbClust(zero);
						Adjusted_GMM.Classify();
						for(unsigned int c=0;c<current_NbClust;c++) Final_Ids.push_back(vector<unsigned int>());
						for(unsigned int i=0;i<new_size;i++){
							if( Adjusted_GMM.getClassificator()[i*2+1] != 0 ) Final_Ids[(unsigned int) Adjusted_GMM.getClassificator()[i*2]].push_back(GBIons[g*3+gbid][i]);
						}
						delete[] adjusted_coords;
					}
					delete[] covar;
					delete[] eigval;
					delete[] eigvec;
				}
				// treat the existing GBs 
				cout << "treating prevs" << endl; 
				for(unsigned int Np=0;Np<indexes_Previous.size();Np++){
						for(unsigned int dimd=0;dimd<6;dimd++) MeanAndCovarClust_temp_prev[gbid].push_back(0.);;
					if( gbid == 0 ){
						GBIons_temp_prev.push_back(vector<unsigned int>());
						AverageCov_temp_prev.push_back(vector<vector<double>>(3,vector<double>(3)));
						AverageMu_temp_prev.push_back(vector<double>(3));
					}
					if( fitting_ok_all[Np*2+gbid] ){
						double *covar = new double[9];
						double *eigval = new double[3];
						double *eigvec = new double[9];
						for(unsigned int d1=0;d1<3;d1++){
							for(unsigned int d2=0;d2<3;d2++) covar[d1*3+d2] = Covars_Identified[Np][d1][d2];
						}
						MT.EigenDecomposition(covar,3,eigval,eigvec);
						bool tostore = false;
						if( ( ( eigval[0] < facClust*eigval[1] ) && ( eigval[0] < eigval[2] ) ) || ( ( eigval[0] < facClust*eigval[2] ) && ( eigval[0] < eigval[1] ) ) ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp_prev[gbid][Np*6+dim] = eigvec[dim];	
							tostore = true;
						}else if( ( eigval[1] < facClust*eigval[0] ) && ( eigval[1] < facClust*eigval[2] ) ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp_prev[gbid][Np*6+dim] = eigvec[3+dim];	
							tostore = true;
						}else if( ( eigval[2] < facClust*eigval[1] ) && ( eigval[2] < facClust*eigval[0] ) ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp_prev[gbid][Np*6+dim] = eigvec[6+dim];	
							tostore = true;
						}
						//else{
						//	fitting_ok_all[Np] = false;	
						//	for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp_prev[gbid].push_back(eigvec[6+dim]);	
						//	for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp_prev[gbid].push_back(Means_Identified[Np][dim]);	
						//	GBIons_temp_prev.push_back(vector<unsigned int>());
						//}
						//if( tostore && gbid == 0 ){
						if( gbid == 0 ){
							Previous_NbClust[indexes_Previous[Np]][0] = current_NbClust;
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp_prev[gbid][Np*6+3+dim] = Means_Identified[Np][dim];	
							for(unsigned int n=0;n<Prev_Final_Ids[Np].size();n++) GBIons_temp_prev[GBIons_temp_prev.size()-1].push_back(Prev_Final_Ids[Np][n]);
							// Store cov and mean for averaging and identifying amorph ions 
							for(unsigned int am_g_1=0;am_g_1<3;am_g_1++){
								AverageMu_temp_prev[AverageMu_temp_prev.size()-1][am_g_1] = Means_Identified[Np][am_g_1];
								for(unsigned int am_g_2=0;am_g_2<3;am_g_2++) AverageCov_temp_prev[AverageCov_temp_prev.size()-1][am_g_1][am_g_2] = Covars_Identified[Np][am_g_1][am_g_2];
							}
						//}else if( tostore && gbid == 1 && fitting_ok_all[Np] ){
						}else if( gbid == 1 && fitting_ok_all[Np*2] ){
							// we already know that the two interfaces are associated no need to compute scalar product
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp_prev[gbid][Np*6+3+dim] = Means_Identified[Np][dim];	
							nbGBAnalyzed_prev += 1;
							current_nbGBAnalyzed_prev += 1;
							GBIons_Final_prev.push_back(vector<unsigned int>());
							GBIons_Final_prev.push_back(vector<unsigned int>());
							GBIons_Final_prev.push_back(vector<unsigned int>());

							// Average covariance and mean of the two interfaces and search amorph ions belonging to this GB 
							AverageCov_prev.push_back(vector<vector<double>>(3,vector<double>(3)));
							AverageMu_prev.push_back(vector<double>(3));
							for(unsigned int am_g_1=0;am_g_1<3;am_g_1++){
								AverageMu_prev[current_nbGBAnalyzed_prev-1][am_g_1] = Means_Identified[Np][am_g_1];
								AverageMu_prev[current_nbGBAnalyzed_prev-1][am_g_1] += AverageMu_temp_prev[Np][am_g_1];
								AverageMu_prev[current_nbGBAnalyzed_prev-1][am_g_1] /= 2.;
								for(unsigned int am_g_2=0;am_g_2<3;am_g_2++){
									AverageCov_prev[current_nbGBAnalyzed_prev-1][am_g_1][am_g_2] = Covars_Identified[Np][am_g_1][am_g_2];
									AverageCov_prev[current_nbGBAnalyzed_prev-1][am_g_1][am_g_2] += AverageCov_temp_prev[Np][am_g_1][am_g_2];
									AverageCov_prev[current_nbGBAnalyzed_prev-1][am_g_1][am_g_2] /= 2.;
								}
							}
							AverageInvCov_prev.push_back(vector<vector<double>>(3,vector<double>(3)));
							IncreaseSmallestEigval(AverageCov_prev[current_nbGBAnalyzed_prev-1], facIncreaseSmallestEigval);
							AverageDetC_prev.push_back(MT.invert3x3(AverageCov_prev[current_nbGBAnalyzed_prev-1],AverageInvCov_prev[current_nbGBAnalyzed_prev-1]));
							// an amorph ion is considered to belong to this GB if its prob is higher than minProbAmorph
							for(unsigned int am=0;am<GBIons[g*3+2].size();am++){
								if( MT.Prob_MultidimGaussian(AverageInvCov_prev[current_nbGBAnalyzed_prev-1],AverageMu_prev[current_nbGBAnalyzed_prev-1],AverageDetC_prev[current_nbGBAnalyzed_prev-1],current_coords_am[am]) > minProbAmorph ) GBIons_Final_prev[(nbGBAnalyzed_prev-1)*3+2].push_back(GBIons[g*3+2][am]);
							}
							if( GBIons_Final_prev[(nbGBAnalyzed_prev-1)*3+2].size() == 0 ) cout << "No amorphous ions have been found" << endl;
							// Store normals and mean positions for this GB and average of normals
							MeanAndCovarClust_Final_prev.push_back(vector<double>());
							for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_Final_prev[nbGBAnalyzed_prev-1].push_back(MeanAndCovarClust_temp_prev[0][6*Np+dim]);
							for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_Final_prev[nbGBAnalyzed_prev-1].push_back(MeanAndCovarClust_temp_prev[1][6*((MeanAndCovarClust_temp_prev[1].size()/6)-1)+dim]);
							double x1 = MeanAndCovarClust_temp_prev[0][6*Np+0];
							double y1 = MeanAndCovarClust_temp_prev[0][6*Np+1];
							double z1 = MeanAndCovarClust_temp_prev[0][6*Np+2];
							double x2 = MeanAndCovarClust_temp_prev[1][6*((MeanAndCovarClust_temp_prev[1].size()/6)-1)+0];
							double y2 = MeanAndCovarClust_temp_prev[1][6*((MeanAndCovarClust_temp_prev[1].size()/6)-1)+1];
							double z2 = MeanAndCovarClust_temp_prev[1][6*((MeanAndCovarClust_temp_prev[1].size()/6)-1)+2];
							if( ( fabs(x1) > 1e-1 && fabs(x2) > 1e-1 && x1*x2 < 0. ) ||  ( fabs(y1) > 1e-1 && fabs(y2) > 1e-1 && y1*y2 < 0. ) || ( fabs(z1) > 1e-1 && fabs(z2) > 1e-1 && z1*z2 < 0. ) ){
								double norm = 0.;
								for(unsigned int dim=0;dim<3;dim++) norm += pow((MeanAndCovarClust_temp_prev[1][6*((MeanAndCovarClust_temp_prev[1].size()/6)-1)+dim]-MeanAndCovarClust_temp_prev[0][6*Np+dim])/2.,2.);
								norm = sqrt(norm);
								for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_Final_prev[nbGBAnalyzed_prev-1].push_back((MeanAndCovarClust_temp_prev[1][6*((MeanAndCovarClust_temp_prev[1].size()/6)-1)+dim]-MeanAndCovarClust_temp_prev[0][6*Np+dim])/(2.*norm));
							}else{
								double norm = 0.;
								for(unsigned int dim=0;dim<3;dim++) norm += pow((MeanAndCovarClust_temp_prev[1][6*((MeanAndCovarClust_temp_prev[1].size()/6)-1)+dim]+MeanAndCovarClust_temp_prev[0][6*Np+dim])/2.,2.);
								for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_Final_prev[nbGBAnalyzed_prev-1].push_back((MeanAndCovarClust_temp_prev[1][6*((MeanAndCovarClust_temp_prev[1].size()/6)-1)+dim]+MeanAndCovarClust_temp_prev[0][6*Np+dim])/(2.*norm));
							}
							// Store ions ID for this GB
							for(unsigned int nb=0;nb<GBIons_temp_prev[Np].size();nb++) GBIons_Final_prev[(nbGBAnalyzed_prev-1)*3].push_back(GBIons_temp_prev[Np][nb]);
							for(unsigned int n=0;n<Prev_Final_Ids[Np].size();n++) GBIons_Final_prev[(nbGBAnalyzed_prev-1)*3+1].push_back(Prev_Final_Ids[Np][n]);

							// Store ids of grain
							IndexesInPreviousVariables.push_back(indexes_Previous[Np]);

							// Update Previous_Variables
							Previous_GBIons[indexes_Previous[Np]][0].clear();
							vector<unsigned int>().swap(Previous_GBIons[indexes_Previous[Np]][0]);
							for(unsigned int n=0;n<GBIons_Final_prev[(nbGBAnalyzed_prev-1)*3].size();n++) Previous_GBIons[indexes_Previous[Np]][0].push_back(GBIons_Final_prev[(nbGBAnalyzed_prev-1)*3][n]);
							Previous_GBIons[indexes_Previous[Np]][1].clear();
							vector<unsigned int>().swap(Previous_GBIons[indexes_Previous[Np]][1]);
							for(unsigned int n=0;n<GBIons_Final_prev[(nbGBAnalyzed_prev-1)*3+1].size();n++) Previous_GBIons[indexes_Previous[Np]][1].push_back(GBIons_Final_prev[(nbGBAnalyzed_prev-1)*3+1][n]);
							Previous_NbClust[indexes_Previous[Np]][1] = current_NbClust;

						}
						//else if( !fitting_ok_all[Np] ){
						//	cout << "Loosing GB because not good shape" << endl;
						//	Previous_GBIons[indexes_Previous[Np]][0].clear(); // set to zero ions
						//	vector<unsigned int>().swap(Previous_GBIons[indexes_Previous[Np]][0]);
						//	Previous_GBIons[indexes_Previous[Np]][1].clear(); // set to zero ions
						//	vector<unsigned int>().swap(Previous_GBIons[indexes_Previous[Np]][1]);
						//	// delete this cluster
						//	//for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_temp_prev[0].erase(MeanAndCovarClust_temp_prev[0].begin());
						//	//GBIons_temp_prev.erase(GBIons_temp_prev.begin());
						//}
						delete[] covar;
						delete[] eigval;
						delete[] eigvec;
					}
					//}else{
					//	cout << "Loosing GB because not enough ions" << endl;
					//	Previous_GBIons[indexes_Previous[Np]][0].clear(); // set to zero ions
					//	vector<unsigned int>().swap(Previous_GBIons[indexes_Previous[Np]][0]);
					//	Previous_GBIons[indexes_Previous[Np]][1].clear(); // set to zero ions
					//	vector<unsigned int>().swap(Previous_GBIons[indexes_Previous[Np]][1]);
					//}
				} // end nclust loop, now we have found the GB and updated Previous_variables
				
				// Now search if new GB
				cout << "Searching if there is new GBs" << endl;	
				for(unsigned int i=0;i<Final_Ids.size();i++){
					bool to_treat = true;
					if( Final_Ids[i].size() >  nbMinIonsInGB && to_treat ){
						double *covar = new double[9];
						double *eigval = new double[3];
						double *eigvec = new double[9];
						for(unsigned int d1=0;d1<3;d1++){
							for(unsigned int d2=0;d2<3;d2++) covar[d1*3+d2] = Adjusted_GMM.getCov()[i*9+d1*3+d2];
						}
						MT.EigenDecomposition(covar,3,eigval,eigvec);
						bool tostore = false;
						if( ( eigval[0] < facClust*eigval[1] ) && ( eigval[0] < facClust*eigval[2] ) ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp[gbid].push_back(eigvec[dim]);	
							tostore = true;
						}else if( ( eigval[1] < facClust*eigval[0] ) && ( eigval[1] < facClust*eigval[2] ) ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp[gbid].push_back(eigvec[3+dim]);	
							tostore = true;
						}else if( ( eigval[2] < facClust*eigval[1] ) && ( eigval[2] < facClust*eigval[0] ) ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp[gbid].push_back(eigvec[6+dim]);	
							tostore = true;
						} 
						if( tostore && gbid == 0 ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp[gbid].push_back(Adjusted_GMM.getMu()[i*3+dim]);	
							GBIons_temp.push_back(vector<unsigned int>());
							for(unsigned int n=0;n<Final_Ids[i].size();n++) GBIons_temp[GBIons_temp.size()-1].push_back(Final_Ids[i][n]);
							NbClustToSaveInPrevious_temp.push_back(Final_Ids.size());
							// Store cov and mean for averaging and identifying amorph ions 
							AverageCov_temp.push_back(vector<vector<double>>(3,vector<double>(3)));
							AverageMu_temp.push_back(vector<double>(3));
							for(unsigned int am_g_1=0;am_g_1<3;am_g_1++){
								AverageMu_temp[AverageMu_temp.size()-1][am_g_1] = Adjusted_GMM.getMu()[i*3+am_g_1];
								for(unsigned int am_g_2=0;am_g_2<3;am_g_2++) AverageCov_temp[AverageCov_temp.size()-1][am_g_1][am_g_2] = Adjusted_GMM.getCov()[i*9+am_g_1*3+am_g_2];
							}
						}else if( tostore && gbid == 1 ){
						        for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_temp[gbid].push_back(Adjusted_GMM.getMu()[i*3+dim]);	
							// search if there is corresponding GB 1 with almost the same normal
							bool isAssociated = false;
							double n0, n1, sp, fullsp;
							vector<unsigned int> indexes;
							vector<double> scalarprod;
							unsigned int index_clust_ass;
							for(unsigned int n=0;n<GBIons_temp.size();n++){
								sp = 0.;
								n0 = 0.;
								n1 = 0.;
								for(unsigned int dim=0;dim<3;dim++){
									sp += MeanAndCovarClust_temp[0][6*n+dim]*MeanAndCovarClust_temp[1][((MeanAndCovarClust_temp[1].size()/6)-1)*6+dim];
									n0 += pow(MeanAndCovarClust_temp[0][6*n+dim],2.);
									n1 += pow(MeanAndCovarClust_temp[1][((MeanAndCovarClust_temp[1].size()/6)-1)*6+dim],2.);
								}
								n0 = sqrt(n0);
								n1 = sqrt(n1);
								fullsp = fabs(sp/(n1*n0));
								if( fullsp > tol_sp ){
									isAssociated = true;
									scalarprod.push_back(fullsp);
									indexes.push_back(n);
								}
							}
							if( isAssociated ){

								index_clust_ass = indexes[MT.max(scalarprod)];
								nbGBAnalyzed += 1;
								current_nbGBAnalyzed += 1;
								GBIons_Final.push_back(vector<unsigned int>());
								GBIons_Final.push_back(vector<unsigned int>());
								GBIons_Final.push_back(vector<unsigned int>());
								NbClustToSaveInPrevious.push_back(vector<unsigned int>());
								NbClustToSaveInPrevious[nbGBAnalyzed-1].push_back(NbClustToSaveInPrevious_temp[index_clust_ass]);
								NbClustToSaveInPrevious[nbGBAnalyzed-1].push_back(Final_Ids.size());

								// Average covariance and mean of the two interfaces and search amorph ions belonging to this GB 
								AverageCov.push_back(vector<vector<double>>(3,vector<double>(3)));
								AverageMu.push_back(vector<double>(3));
								for(unsigned int am_g_1=0;am_g_1<3;am_g_1++){
									AverageMu[current_nbGBAnalyzed-1][am_g_1] = Adjusted_GMM.getMu()[i*3+am_g_1];
									AverageMu[current_nbGBAnalyzed-1][am_g_1] += AverageMu_temp[index_clust_ass][am_g_1];
									AverageMu[current_nbGBAnalyzed-1][am_g_1] /= 2.;
									for(unsigned int am_g_2=0;am_g_2<3;am_g_2++){
										AverageCov[current_nbGBAnalyzed-1][am_g_1][am_g_2] = Adjusted_GMM.getCov()[i*9+am_g_1*3+am_g_2];
										AverageCov[current_nbGBAnalyzed-1][am_g_1][am_g_2] += AverageCov_temp[index_clust_ass][am_g_1][am_g_2];
										AverageCov[current_nbGBAnalyzed-1][am_g_1][am_g_2] /= 2.;
									}
								}
								AverageInvCov.push_back(vector<vector<double>>(3,vector<double>(3)));
								IncreaseSmallestEigval(AverageCov[current_nbGBAnalyzed-1], facIncreaseSmallestEigval);
								AverageDetC.push_back(MT.invert3x3(AverageCov[current_nbGBAnalyzed-1],AverageInvCov[current_nbGBAnalyzed-1]));
								// an amorph ion is considered to belong to this GB if its prob is higher than minProbAmorph
								for(unsigned int am=0;am<GBIons[g*3+2].size();am++){
									if( MT.Prob_MultidimGaussian(AverageInvCov[current_nbGBAnalyzed-1],AverageMu[current_nbGBAnalyzed-1],AverageDetC[current_nbGBAnalyzed-1],current_coords_am[am]) > minProbAmorph ) GBIons_Final[(nbGBAnalyzed-1)*3+2].push_back(GBIons[g*3+2][am]);
								}
								if( GBIons_Final[(nbGBAnalyzed-1)*3+2].size() == 0 ) cout << "No amorphous ions have been found" << endl;
									cout << "Cov1 = " << endl;
									for(unsigned int i1=0;i1<3;i1++){
										for(unsigned int i2=0;i2<3;i2++) cout << AverageCov_temp[index_clust_ass][i1][i2] << " ";
										cout << endl;
									}
									cout << "Mean1 = ";
									for(unsigned int i2=0;i2<3;i2++) cout << AverageMu_temp[index_clust_ass][i2] << " ";
									cout << endl;
									cout << "Cov2 = " << endl;
									for(unsigned int i1=0;i1<3;i1++){
										for(unsigned int i2=0;i2<3;i2++) cout << Adjusted_GMM.getCov()[i*9+i1*3+i2] << " ";
										cout << endl;
									}
									cout << "Mean1 = ";
									for(unsigned int i2=0;i2<3;i2++) cout << Adjusted_GMM.getMu()[i*3+i2] << " ";
									cout << endl;
									cout << "AveCov = " << endl;
									for(unsigned int i1=0;i1<3;i1++){
										for(unsigned int i2=0;i2<3;i2++) cout << AverageCov[current_nbGBAnalyzed-1][i1][i2] << " ";
										cout << endl;
									}
									cout << "Mean1 = ";
									for(unsigned int i2=0;i2<3;i2++) cout << AverageMu[current_nbGBAnalyzed-1][i2] << " ";
									cout << endl;


								// Store normals and mean positions for this GB and average of normals
								MeanAndCovarClust_Final.push_back(vector<double>());
								for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back(MeanAndCovarClust_temp[0][6*index_clust_ass+dim]);
								for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back(MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]);
								double x1 = MeanAndCovarClust_temp[0][6*index_clust_ass+0];
								double y1 = MeanAndCovarClust_temp[0][6*index_clust_ass+1];
								double z1 = MeanAndCovarClust_temp[0][6*index_clust_ass+2];
								double x2 = MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+0];
								double y2 = MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+1];
								double z2 = MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+2];
								if( ( fabs(x1) > 1e-1 && fabs(x2) > 1e-1 && x1*x2 < 0. ) ||  ( fabs(y1) > 1e-1 && fabs(y2) > 1e-1 && y1*y2 < 0. ) || ( fabs(z1) > 1e-1 && fabs(z2) > 1e-1 && z1*z2 < 0. ) ){
									double norm = 0.;
									for(unsigned int dim=0;dim<3;dim++) norm += pow((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]-MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/2.,2.);
									norm = sqrt(norm);
									for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]-MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/(2.*norm));
								}else{
									double norm = 0.;
									for(unsigned int dim=0;dim<3;dim++) norm += pow((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]+MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/2.,2.);
									for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]+MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/(2.*norm));
								}
								// Store ions ID for this GB
								for(unsigned int nb=0;nb<GBIons_temp[index_clust_ass].size();nb++) GBIons_Final[(nbGBAnalyzed-1)*3].push_back(GBIons_temp[index_clust_ass][nb]);
								for(unsigned int nb=0;nb<Final_Ids[i].size();nb++) GBIons_Final[(nbGBAnalyzed-1)*3+1].push_back(Final_Ids[i][nb]);

								// Stored ids of grain
								GBId_arr_Final.push_back(GBId_arr[g]);

								// delete this cluster
								for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_temp[0].erase(MeanAndCovarClust_temp[0].begin()+6*index_clust_ass);
								GBIons_temp.erase(GBIons_temp.begin()+index_clust_ass);
							}
						}
						delete[] covar;
						delete[] eigval;
						delete[] eigvec;
					}
				} // end nclust loop


				// free memory
				if( indexes_Previous.size() != 0 ){
					delete[] ClusterIndexAssociated;
					for(unsigned int dimc=0;dimc<indexes_Previous.size();dimc++){
						Prev_Final_Ids[dimc].clear();
						vector<unsigned int>().swap(Prev_Final_Ids[dimc]);
					}
					delete[] Prev_Final_Ids;
				}
			} // end gbid loop
			for(unsigned int gbid=0;gbid<3;gbid++) delete[] current_coords[gbid];
			current_coords.clear();
			vector<double*>().swap(current_coords);
			if( ClusterIdentifiedToDelete ){
				delete[] fitting_ok_all;
				for(unsigned int dimc=0;dimc<ClusterIdentifiedSize;dimc++){
					Means_Identified[dimc].clear();
					vector<double>().swap(Means_Identified[dimc]);
					for(unsigned int dimk=0;dimk<Covars_Identified[dimc].size();dimk++){
						Covars_Identified[dimc][dimk].clear();
						vector<double>().swap(Covars_Identified[dimc][dimk]);
					}
					vector<vector<double>>().swap(Covars_Identified[dimc]);
				}
				delete[] Means_Identified;
				delete[] Covars_Identified;
			}
		} // end g loop
							cout << 6 << endl;
		
		double *MeanStrain_1 = new double[8];
		double *MeanStrain_2 = new double[8];
		double *MeanStrain_am = new double[8];
		long double *MeanStress_1 = new long double[6];
		long double *MeanStress_2 = new long double[6];
		long double *MeanStress_am = new long double[6];
		double am_thick;

		// print GBs already in previous
		for(unsigned int i=0;i<nbGBAnalyzed_prev;i++){

			AverageStressAndStrain(MySystem, MeanStrain_1, MeanStrain_2, MeanStrain_am, MeanStress_1, MeanStress_2, MeanStress_am, MeanAndCovarClust_Final_prev, GBIons_Final_prev, AtStrain_ind, size_AtStrain, Stress_ind, size_Stress, AtVol_ind, size_AtVol, i, am_thick);

			string path_1 = "./GB_";
			string path_2 = "/";
			string fullpath = path_1+to_string(IndexesInPreviousVariables[i])+path_2;
			const char *env = fullpath.c_str();
        		DIR *dir;
        		if( (dir = opendir(env) ) == nullptr ) cout << "Warning the directory does not exist" << endl;

			string datafilename = fullpath+OutputFilename;
			std::ofstream f_out_res;
			f_out_res.open(datafilename.c_str(), std::ios_base::app);
			if( f_out_res.fail() ) cout << "Issue when printing" << endl;
			f_out_res << timesteps[t] << " " << MeanAndCovarClust_Final_prev[i][12] << " " << MeanAndCovarClust_Final_prev[i][13] << " " << MeanAndCovarClust_Final_prev[i][14] << " " << GBIdOverTimestep[IndexesInPreviousVariables[i]] - (GBIdOverTimestep[IndexesInPreviousVariables[i]]-floor(GBIdOverTimestep[IndexesInPreviousVariables[i]])) << " " << (int) 10.*(GBIdOverTimestep[IndexesInPreviousVariables[i]]-floor(GBIdOverTimestep[IndexesInPreviousVariables[i]])) << " " << am_thick ;
			for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_1[ii]; 
			for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_2[ii]; 
			for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_am[ii]; 
			for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_1[ii]; 
			for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_2[ii]; 
			for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_am[ii]; 
		        f_out_res << " " << GBIons_Final_prev[3*i].size() + GBIons_Final_prev[3*i+1].size() + GBIons_Final_prev[3*i+2].size();
			f_out_res << endl;
			f_out_res.close();

			string name = "GB_";
			string ext = ".cfg";
			string fullname = fullpath+name+timesteps[t]+ext;
			printAtomicFile(MySystem, fullname, timesteps[t], GBIons_Final_prev, i);
		}

		// PRINT NEW GB	
		cout << "Printing new GBs" << endl;
		for(unsigned int i=0;i<nbGBAnalyzed;i++){

			AverageStressAndStrain(MySystem, MeanStrain_1, MeanStrain_2, MeanStrain_am, MeanStress_1, MeanStress_2, MeanStress_am, MeanAndCovarClust_Final, GBIons_Final, AtStrain_ind, size_AtStrain, Stress_ind, size_Stress, AtVol_ind, size_AtVol, i, am_thick);

			string path_1 = "./GB_";
			string path_2 = "/";
			string fullpath = path_1+to_string(GBIdOverTimestep.size()+i)+path_2;
			const char *env = fullpath.c_str();
        		DIR *dir;
        		if( (dir = opendir(env) ) != nullptr ) cout << "Warning the directory already exist" << endl;
			else std::filesystem::create_directory(fullpath);

			string datafilename = fullpath+OutputFilename;
			std::ofstream f_out_res;
			f_out_res.open(datafilename.c_str());
			if( f_out_res.fail() ) cout << "Issue when printing" << endl;
			f_out_res << "timestep(0) nx(1) ny(2) nz(3) G1(4) G2(5) AmThick(6) Eps_xx_1(7) Eps_yy_1(8) Eps_zz_1(9) Eps_xy_1(10) Eps_xz_1(11) Eps_yz_1(12) ShearInv_1(13) SphInv_1(14) Eps_xx_2(15) Eps_yy_2(16) Eps_zz_2(17) Eps_xy_2(18) Eps_xz_2(19) Eps_yz_2(20) ShearInv_2(21) SphInv_2(22) Eps_xx_am(23) Eps_yy_am(24) Eps_zz_am(25) Eps_xy_am(26) Eps_xz_am(27) Eps_yz_am(28) ShearInv_am(29) SphInv_am(30) sigma_xx_1(31) sigma_yy_1(32) sigma_zz_1(33) sigma_xy_1(34) sigma_xz_1(35) sigma_yz_1(36) sigma_xx_2(37) sigma_yy_2(38) sigma_zz_2(39) sigma_xy_2(40) sigma_xz_2(41) sigma_yz_2(42) sigma_xx_am(43) sigma_yy_am(44) sigma_zz_am(45) sigma_xy_am(46) sigma_xz_am(47) sigma_yz_am(48) nbAt49)" << endl;
			f_out_res << timesteps[t] << " " << MeanAndCovarClust_Final[i][12] << " " << MeanAndCovarClust_Final[i][13] << " " << MeanAndCovarClust_Final[i][14] << " " << GBId_arr_Final[i] - (GBId_arr_Final[i]-floor(GBId_arr_Final[i])) << " " << (int) 10.*(GBId_arr_Final[i]-floor(GBId_arr_Final[i])) << " " << am_thick ;
			for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_1[ii]; 
			for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_2[ii]; 
			for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_am[ii]; 
			for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_1[ii]; 
			for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_2[ii]; 
			for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_am[ii]; 
		        f_out_res << " " << GBIons_Final[3*i].size() + GBIons_Final[3*i+1].size() + GBIons_Final[3*i+2].size();
			f_out_res << endl;
			f_out_res.close();

			string name = "GB_";
			string ext = ".cfg";
			string fullname = fullpath+name+timesteps[t]+ext;
			printAtomicFile(MySystem, fullname, timesteps[t], GBIons_Final, i);
		}
		cout << "Saving previous variables" << endl;
		// Save new GBs in Previous_Variables 
		for(unsigned int i=0;i<nbGBAnalyzed;i++){
			GBIdOverTimestep.push_back(GBId_arr_Final[i]);
			Previous_GBIons.push_back(vector<vector<unsigned int>>());
			Previous_NbClust.push_back(vector<unsigned int>());
			for(unsigned int gbid=0;gbid<2;gbid++){
				Previous_NbClust[GBIdOverTimestep.size()-1].push_back(NbClustToSaveInPrevious[i][gbid]);
				Previous_GBIons[GBIdOverTimestep.size()-1].push_back(vector<unsigned int>());
				for(unsigned int n=0;n<GBIons_Final[i*3+gbid].size();n++) Previous_GBIons[GBIdOverTimestep.size()-1][gbid].push_back(GBIons_Final[i*3+gbid][n]);
			}
		}
		cout << "Size of prevs : " << GBIdOverTimestep.size() << " " << Previous_GBIons.size() << endl;

		delete[] MeanStrain_1;
		delete[] MeanStrain_2;
		delete[] MeanStrain_am;
		delete[] MeanStress_1;
		delete[] MeanStress_2;
		delete[] MeanStress_am;

	} // end timestep loop
	// ### END TREAT ALL TIMESTEPS
	Dis.ExecutionTime();	
	return 0;
}
