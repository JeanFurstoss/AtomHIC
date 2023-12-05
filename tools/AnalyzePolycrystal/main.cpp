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
	string InputFilename, OutputFilename;
	double SiO_d;
	double SecFac = 2.;
	if( argc == 3 ){
		InputFilename = argv[1];
		OutputFilename = argv[2];
	}else{
		cerr << "Usage: ./IdentifyGB InputFilename OutputFilename" << endl;
		cerr << "the GB field returned here consists of G1Id.G2Id => when the number of grain is higher than 9 there is doubt and this program shoudl be modified" << endl;
		return EXIT_FAILURE;
	}
	AtomicSystem MySystem(InputFilename);

	const unsigned int nbAt = MySystem.getNbAtom();
	const unsigned int nbNMax = MySystem.getNbMaxN();

	double *GBAux = new double[nbAt];

	vector<double> Struct_GB;
	Struct_GB.push_back(0); // amorph
	Struct_GB.push_back(5); // interface

	unsigned int size_struct;
	unsigned struct_ind = MySystem.getAuxIdAndSize("struct",size_struct);
	unsigned int size_gt;
	unsigned gt_ind = MySystem.getAuxIdAndSize("grainID",size_gt);

	unsigned int nbGrain = 8; // TODO generalize	
	
	MySystem.searchNeighbours(10.);
	
	// Paralelize
	for(unsigned int i=0;i<nbAt;i++) {
		MathTools MT;
		bool IsGB = false;
		for(unsigned int s=0;s<Struct_GB.size();s++){
			if( MySystem.getAux(struct_ind)[i*size_struct+1] == Struct_GB[s] ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
				IsGB = true;
				break;
			}
		}
		if( IsGB ){
			vector<vector<double>> Dist;
			for(unsigned int d=0;d<nbGrain;d++) Dist.push_back(vector<double>());
			vector<unsigned int> GrainId;
			bool IsNGB, IsGrainAlreadyStored;
			unsigned int id;
			for(unsigned int j=0;j<MySystem.getNeighbours(i*(nbNMax+1));j++){
				IsNGB = false;
				id = MySystem.getNeighbours(i*(nbNMax+1)+j+1);
				for(unsigned int s=0;s<Struct_GB.size();s++){
					if( MySystem.getAux(struct_ind)[id*size_struct+1] == Struct_GB[s] ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
						IsNGB = true;
						break;
					}
				}
				if( !IsNGB ){
					IsGrainAlreadyStored = false;
					for(unsigned int k=0;k<GrainId.size();k++){
						if( (unsigned int) MySystem.getAux(gt_ind)[id] == GrainId[k] ){
							IsGrainAlreadyStored = true;
							break;
						}
					}
					if( !IsGrainAlreadyStored ){
						GrainId.push_back((unsigned int) MySystem.getAux(gt_ind)[id]);
					}
					double xp = MySystem.getWrappedPos(id).x+MySystem.getCLNeighbours(i*nbNMax*3+j*3)*MySystem.getH1()[0]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+1)*MySystem.getH2()[0]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+2)*MySystem.getH3()[0]-MySystem.getWrappedPos(i).x;
					double yp = MySystem.getWrappedPos(id).y+MySystem.getCLNeighbours(i*nbNMax*3+j*3)*MySystem.getH1()[1]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+1)*MySystem.getH2()[1]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+2)*MySystem.getH3()[1]-MySystem.getWrappedPos(i).y;
					double zp = MySystem.getWrappedPos(id).z+MySystem.getCLNeighbours(i*nbNMax*3+j*3)*MySystem.getH1()[2]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+1)*MySystem.getH2()[2]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+2)*MySystem.getH3()[2]-MySystem.getWrappedPos(i).z;
					Dist[(unsigned int) MySystem.getAux(gt_ind)[id] - 1].push_back(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)); 
				}
			}
			vector<double> FinalArray;
			for(unsigned int n=0;n<GrainId.size();n++){
				FinalArray.push_back(MT.min(Dist[GrainId[n]]));
				FinalArray.push_back(GrainId[n]);
			}
			MT.sort(FinalArray,0,2,FinalArray);
			if( FinalArray[1] < FinalArray[3] ) GBAux[i] = FinalArray[1]+(FinalArray[3]/10.);
			else GBAux[i] = FinalArray[3]+(FinalArray[1]/10.);
		} else GBAux[i] = 0.;
	}

	MySystem.setAux(GBAux,"GBId");	
	MySystem.printSystem_aux(OutputFilename,"GBId grainID struct");
	cout << "New system printed" << endl;
	delete[] GBAux;
	return 0;
}
