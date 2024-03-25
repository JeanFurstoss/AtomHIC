// AtomHic library files
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

using namespace std;

int main(int argc, char *argv[])
{
	// Here we use a dump file containing (1) the structure of a given ion (either interface, amorph or crystal), (2) the GB Id, (3) the atomic volume, (4) strain and (5) stress tensor
	// the objective is to return a file for each GB Id with:
	// nx ny nz AmorphousThickness GBSurface Interface1Strain(AndStress)Tensor Interface2Strain(AndStress)Tensor AmorphousStrain(AndStress)Tensor Interface1Strain(AndStress)Invariant Interface2Strain(AndStress)Invariants AmorphousStrain(AndStress)Invariants 
	string InputFilename, OutputFilename;


	if( argc == 3 ){
		InputFilename = argv[1];
		OutputFilename = argv[2];
	}else{
		cerr << "Usage: ./AnalyzePolycrystal InputFilename OutputFilename" << endl;
		cerr << "TODO description" << endl;
		return EXIT_FAILURE;
	}
	AtomicSystem MySystem(InputFilename);

	const unsigned int nbAt = MySystem.getNbAtom();

	unsigned int Struct_GB = 5; //TODO add in argument of exe
	unsigned int Struct_Amorph = 0;
	
	// Get the different aux properties
	unsigned int size_Struct;
	unsigned Struct_ind = MySystem.getAuxIdAndSize("Struct",size_Struct);
	unsigned int size_GBId;
	unsigned GBId_ind = MySystem.getAuxIdAndSize("GBId",size_GBId);
	unsigned int size_AtVol;
	unsigned AtVol_ind = MySystem.getAuxIdAndSize("AtVol",size_AtVol);
	unsigned int size_Stress;
	unsigned Stress_ind = MySystem.getAuxIdAndSize("c_s",size_Stress);
	unsigned int size_AtStrain;
	unsigned AtStrain_ind = MySystem.getAuxIdAndSize("AtomicStrain",size_AtStrain);

	// Search which ions belongs to GB and differentiate the one in amorph and the one in interface
	vector<double> GBId_arr;
	vector<vector<unsigned int>> GBIons;
	cout << "Searching GB ions" << endl;
	// TODO Paralelize?
	for(unsigned int i=0;i<nbAt;i++) {
		bool IsGB = false;
		unsigned int GBOrAmorph;
		if( (unsigned int) MySystem.getAux(Struct_ind)[i*size_Struct] == Struct_GB ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
			GBOrAmorph = 0;       
			IsGB = true;
		}else if( (unsigned int) MySystem.getAux(Struct_ind)[i*size_Struct] == Struct_Amorph ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
			GBOrAmorph = 1;       
			IsGB = true;
		}
		if( IsGB ){
			bool IsAlreadyStored = false;
			for(unsigned int g=0;g<GBId_arr.size();g++){
				if( MySystem.getAux(GBId_ind)[i*size_GBId] == GBId_arr[g] ){
					GBIons[g*2+GBOrAmorph].push_back(i);
					IsAlreadyStored = true;
					break;
				}
			}
			if( !IsAlreadyStored ){
				GBId_arr.push_back(MySystem.getAux(GBId_ind)[i*size_GBId]);
				GBIons.push_back(vector<unsigned int>());
				GBIons.push_back(vector<unsigned int>());
				GBIons[(GBId_arr.size()-1)*2+GBOrAmorph].push_back(i); //not sure about the minus one (to verify)
			}
		}
	}
	cout << "Done" << endl;
	cout << "Compute GB plane normals" << endl;
	// Compute the normal plane to each GB using only interface ions (maybe compute here stress and strain for interfaces)
	// Pb: multiple planes for a given GBId (i.e. an other one formed by the periodic replica of the grain)
	vector<vector<double>> coords;
	double *Normals = new double[GBId_arr.size()*4];
	MathTools MT;

	double afit,bfit,cfit;
	for(unsigned int g=0;g<GBId_arr.size();g++){
		coords.clear();
		for(unsigned int i=0;i<GBIons[g*2].size();i++){
			coords.push_back(vector<double>());
			coords[i].push_back(MySystem.getWrappedPos(GBIons[g*2][i]).x);
			coords[i].push_back(MySystem.getWrappedPos(GBIons[g*2][i]).y);
			coords[i].push_back(MySystem.getWrappedPos(GBIons[g*2][i]).z);
		}
		MT.plane_fit(coords,afit,bfit,cfit);
		// Normalize
		double normFac = sqrt(pow(afit,2.)+pow(bfit,2.)+1);
		Normals[g*4] = afit/normFac;
		Normals[g*4+1] = bfit/normFac;
		Normals[g*4+2] = -1./normFac;
		Normals[g*4+3] = cfit/normFac;
		cout << "GBId : " << GBId_arr[g] << ", plane coeff are : " << Normals[g*4] << ", " << Normals[g*4+1] << ", " << Normals[g*4+2] << ", " << Normals[g*4+3] << ", nb of ion considered : " << GBIons[g*2].size() << endl; 
	}

	cout << "Done" << endl;

	// z = ax+by+c <=> apx+bpy+cpz+dp=0 and a,b,c normed
	// apx+bpy+cpax+cpbby+ccp+dp=0
	//Points[count].push_back(((double) x));
	//Points[count].push_back(((double) y));
	//Points[count].push_back((a*(double) x)+(b*(double) y)+c+(((double) (rand()%100) - 50.)/10));


	delete[] Normals;
	return 0;
}
