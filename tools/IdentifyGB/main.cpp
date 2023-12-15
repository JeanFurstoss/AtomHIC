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
	double rcut_large;
	double SecFac = 2.;
	if( argc == 4 ){
		InputFilename = argv[1];
		OutputFilename = argv[2];
		istringstream iss_rc(argv[3]);
		iss_rc >> rcut_large;
	}else{
		cerr << "Usage: ./IdentifyGB InputFilename OutputFilename rcut_large" << endl;
		cerr << "the GB field returned here consists of G1Id.G2Id => when the number of grain is higher than 9 there is doubt and this program shoudl be modified" << endl;
		return EXIT_FAILURE;
	}
	AtomicSystem MySystem(InputFilename);

	const unsigned int nbAt = MySystem.getNbAtom();

	double *GBAux = new double[nbAt];
	for(unsigned int i=0;i<nbAt;i++) GBAux[i] = 0.;

	vector<double> Struct_GB;
	unsigned int nbStruct = 7;
	Struct_GB.push_back(0); // amorph
	Struct_GB.push_back(5); // interface

	unsigned int size_struct;
	unsigned struct_ind = MySystem.getAuxIdAndSize("Struct",size_struct);
	unsigned int size_gt;
	unsigned gt_ind = MySystem.getAuxIdAndSize("grainID",size_gt);

	unsigned int nbGrain = 8; // TODO generalize
	double rc_id_edge = 5.; // cutoff radius for search the edge of the grains
	MathTools MT;
	double PCIonsFrac = 0.75;
	double EdgeIonsFrac = 0.25;

	// In a first step, we search the edge atoms (i.e. the atoms at the frontiers of the grain in order to stringly restrict the following neighbour research for which the cutoff radius have to be very long=
	// this identification of edge atoms will also help to identify atoms misidentified as crystal (with very few neighbours)
	Crystal *MyCrystal = new Crystal("Forsterite");
	AtomicSystem *SystemG = new AtomicSystem();
	vector<unsigned int> IdAtGrain;
	vector<unsigned int> EdgeGrainId;
	vector<unsigned int> IsolGrainId;
	for(unsigned int g=0;g<nbGrain;g++){
		IdAtGrain.clear();
		// Search if the ion belongs to the considered grain
		for(unsigned int i=0;i<nbAt;i++){
			if( MySystem.getAux(gt_ind)[i*size_gt] == g+1 ){
				bool IsGB = false;
				for(unsigned int s=0;s<Struct_GB.size();s++){
					if( MySystem.getAux(struct_ind)[i*size_struct] == Struct_GB[s] ){
						IsGB = true;
						EdgeGrainId.push_back(i);
						break;
					}
				}
				if( !IsGB ) IdAtGrain.push_back(i);
			}
		}
		// Create an atomic system with this isolated grain
		Atom *AtomsG = new Atom[IdAtGrain.size()]; // TODO dont forget to delete
		unsigned int *NbNeight = new unsigned int[IdAtGrain.size()];
		for(unsigned int i=0;i<IdAtGrain.size();i++) AtomsG[i] = MySystem.getAtom(IdAtGrain[i]);
		SystemG->AtomListConstructor(AtomsG,IdAtGrain.size(),MyCrystal,MySystem.getH1(),MySystem.getH2(),MySystem.getH3());
		// Compute the neighbours and search if a given ion have medium neighbour number (=> edge ion) or low neighbour number (=> isolated and maybe mislabelled ion) compared to the maximum number of neighbour which is assumed to be the perfect crystal
		cout << "Grain " << g+1 << endl;
		unsigned int nbNMaxG = SystemG->searchNeighbours(rc_id_edge);
		for(unsigned int i=0;i<IdAtGrain.size();i++) NbNeight[i] = SystemG->getNeighbours(i*(nbNMaxG+1));
		unsigned int NMax = MT.max(NbNeight,IdAtGrain.size());
		for(unsigned int i=0;i<IdAtGrain.size();i++){
			if( NbNeight[i] > NMax*PCIonsFrac ){
				GBAux[IdAtGrain[i]] = 0.;
			}else if( NbNeight[i] > NMax*EdgeIonsFrac ){
				EdgeGrainId.push_back(IdAtGrain[i]);
				GBAux[IdAtGrain[i]] = 0.;
			}else{ // Isolated ions
				EdgeGrainId.push_back(IdAtGrain[i]);
				IsolGrainId.push_back(IdAtGrain[i]);
			}
		}
		delete[] AtomsG;
		delete[] NbNeight;
	}

	SystemG->~AtomicSystem();

	// search for isolated ions the structure the most represented in its neighbour and change it own struct according to
	cout << "Treating restricted system" << endl;
	cout << "Homogeneize interfaces" << endl;
	const unsigned int nbNMax = MySystem.searchNeighbours(rc_id_edge);
	unsigned int *count_struct = new unsigned int[nbStruct];
	for(unsigned int i=0;i<IsolGrainId.size();i++){
		for(unsigned int k=0;k<nbStruct;k++) count_struct[k] = 0;
		for(unsigned int j=0;j<MySystem.getNeighbours(IsolGrainId[i]*(nbNMax+1));j++) count_struct[(unsigned int) MySystem.getAux(struct_ind)[MySystem.getNeighbours(IsolGrainId[i]*(nbNMax+1)+j+1)*size_struct]] += 1;
		unsigned int ind_max = MT.max_ind(count_struct,nbStruct);
		MySystem.getAux(struct_ind)[IsolGrainId[i]*size_struct] = ind_max;
	}

	MySystem.deleteNeighList(); // free memory

	// Create new system with only edge ions and set struct and grainID auxiliary
	vector<unsigned int> ResSys_ToSearch;
	vector<unsigned int> ResSys_ForSearch;
	Atom *RestrictedAtoms = new Atom[EdgeGrainId.size()];
	double *AuxStruct = new double[EdgeGrainId.size()];
	double *AuxGT = new double[EdgeGrainId.size()];
	double *AuxRS = new double[EdgeGrainId.size()];
	for(unsigned int i=0;i<EdgeGrainId.size();i++){
		RestrictedAtoms[i] = MySystem.getAtom(EdgeGrainId[i]);
		AuxStruct[i] = MySystem.getAux(struct_ind)[EdgeGrainId[i]*size_struct];
		AuxGT[i] = MySystem.getAux(gt_ind)[EdgeGrainId[i]*size_gt];
		bool IsGB = false;
		for(unsigned int s=0;s<Struct_GB.size();s++){
			if( AuxStruct[i] == Struct_GB[s] ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
				IsGB = true;
				ResSys_ToSearch.push_back(i);
				AuxRS[i] = 0.;
				break;
			}
		}
		if( !IsGB ){
			ResSys_ForSearch.push_back(i);
			AuxRS[i] = 1.;
			GBAux[EdgeGrainId[i]] = 0.;
		}
	}
	AtomicSystem RestrictedSystem(RestrictedAtoms,EdgeGrainId.size(),MyCrystal,MySystem.getH1(),MySystem.getH2(),MySystem.getH3());
	RestrictedSystem.setAux(AuxStruct,"Struct");
	RestrictedSystem.setAux(AuxGT,"grainID");
	RestrictedSystem.setAux(AuxRS,"Neigh");
	RestrictedSystem.printSystem_aux("ResSys.cfg","Struct grainID Neigh");

	delete[] AuxStruct;	
	delete[] AuxGT;	
	delete[] AuxRS;

	unsigned int size_struct_res;
	unsigned struct_ind_res = RestrictedSystem.getAuxIdAndSize("Struct",size_struct_res);
	unsigned int size_gt_res;
	unsigned gt_ind_res = RestrictedSystem.getAuxIdAndSize("grainID",size_gt_res);
	cout << "Compute GB identifier" << endl;
	const unsigned long int nbNMax_res = RestrictedSystem.searchNeighbours_restricted(rcut_large,ResSys_ToSearch,ResSys_ForSearch);
	const unsigned long int one_l = 1;	
	const unsigned long int two_l = 2;	
	const unsigned long int three_l = 3;	
	// Paralelize
	#pragma omp parallel for
	for(unsigned long int i=0;i<ResSys_ToSearch.size();i++) {
		double xpos = RestrictedSystem.getWrappedPos(ResSys_ToSearch[i]).x;
		double ypos = RestrictedSystem.getWrappedPos(ResSys_ToSearch[i]).y;
		double zpos = RestrictedSystem.getWrappedPos(ResSys_ToSearch[i]).z;
		vector<vector<double>> Dist;
		for(unsigned int d=0;d<nbGrain;d++) Dist.push_back(vector<double>());
		vector<unsigned int> GrainId;
		bool IsNGB, IsGrainAlreadyStored;
		unsigned int id;
		for(unsigned long int j=0;j<RestrictedSystem.getNeighbours(i*(nbNMax_res+one_l));j++){
			id = RestrictedSystem.getNeighbours(i*(nbNMax_res+one_l)+j+one_l);
			IsGrainAlreadyStored = false;
			for(unsigned int k=0;k<GrainId.size();k++){
				if( (unsigned int) RestrictedSystem.getAux(gt_ind_res)[id] == GrainId[k] ){
					IsGrainAlreadyStored = true;
					break;
				}
			}
			if( !IsGrainAlreadyStored ) GrainId.push_back((unsigned int) RestrictedSystem.getAux(gt_ind_res)[id]);
			double xp = RestrictedSystem.getWrappedPos(id).x+RestrictedSystem.getCLNeighbours(i*nbNMax_res*three_l+j*three_l)*RestrictedSystem.getH1()[0]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*three_l+j*three_l+one_l)*RestrictedSystem.getH2()[0]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*three_l+j*three_l+two_l)*RestrictedSystem.getH3()[0]-xpos;
			double yp = RestrictedSystem.getWrappedPos(id).y+RestrictedSystem.getCLNeighbours(i*nbNMax_res*three_l+j*three_l)*RestrictedSystem.getH1()[1]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*three_l+j*three_l+one_l)*RestrictedSystem.getH2()[1]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*three_l+j*three_l+two_l)*RestrictedSystem.getH3()[1]-ypos;
			double zp = RestrictedSystem.getWrappedPos(id).z+RestrictedSystem.getCLNeighbours(i*nbNMax_res*three_l+j*three_l)*RestrictedSystem.getH1()[2]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*three_l+j*three_l+one_l)*RestrictedSystem.getH2()[2]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*three_l+j*three_l+two_l)*RestrictedSystem.getH3()[2]-zpos;
			Dist[(unsigned int) RestrictedSystem.getAux(gt_ind_res)[id] - 1].push_back(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)); 
		}
		if( GrainId.size() < 2 ){
			cout << "Issue ! Increase cutoff radius !" << endl;
			GBAux[EdgeGrainId[ResSys_ToSearch[i]]] = 0.;
		}else{
			vector<double> FinalArray;
			for(unsigned int n=0;n<GrainId.size();n++){
				FinalArray.push_back(MT.min_vec(Dist[GrainId[n]-1]));
				FinalArray.push_back(GrainId[n]);
			}
			MT.sort(FinalArray,0,2,FinalArray);
			if( FinalArray[1] < FinalArray[3] ) GBAux[EdgeGrainId[ResSys_ToSearch[i]]] = FinalArray[1]+(FinalArray[3]/10.);
			else GBAux[EdgeGrainId[ResSys_ToSearch[i]]] = FinalArray[3]+(FinalArray[1]/10.);
		}
	}
	//for(unsigned int i=0;i<EdgeGrainId.size();i++) {
	//	bool IsGB = false;
	//	for(unsigned int s=0;s<Struct_GB.size();s++){
	//		if( RestrictedSystem.getAux(struct_ind_res)[i*size_struct_res] == Struct_GB[s] ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
	//			IsGB = true;
	//			break;
	//		}
	//	}
	//	if( IsGB ){
	//		double xpos = RestrictedSystem.getWrappedPos(i).x;
	//		double ypos = RestrictedSystem.getWrappedPos(i).y;
	//		double zpos = RestrictedSystem.getWrappedPos(i).z;
	//		vector<vector<double>> Dist;
	//		for(unsigned int d=0;d<nbGrain;d++) Dist.push_back(vector<double>());
	//		vector<unsigned int> GrainId;
	//		bool IsNGB, IsGrainAlreadyStored;
	//		unsigned int id;
	//		for(unsigned int j=0;j<RestrictedSystem.getNeighbours(i*(nbNMax_res+1));j++){
	//			IsNGB = false;
	//			id = RestrictedSystem.getNeighbours(i*(nbNMax_res+1)+j+1);
	//			for(unsigned int s=0;s<Struct_GB.size();s++){
	//				if( RestrictedSystem.getAux(struct_ind_res)[id*size_struct_res] == Struct_GB[s] ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
	//					IsNGB = true;
	//					break;
	//				}
	//			}
	//			if( !IsNGB ){
	//				IsGrainAlreadyStored = false;
	//				for(unsigned int k=0;k<GrainId.size();k++){
	//					if( (unsigned int) RestrictedSystem.getAux(gt_ind_res)[id] == GrainId[k] ){
	//						IsGrainAlreadyStored = true;
	//						break;
	//					}
	//				}
	//				if( !IsGrainAlreadyStored ) GrainId.push_back((unsigned int) RestrictedSystem.getAux(gt_ind_res)[id]);
	//				double xp = RestrictedSystem.getWrappedPos(id).x+RestrictedSystem.getCLNeighbours(i*nbNMax_res*3+j*3)*RestrictedSystem.getH1()[0]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*3+j*3+1)*RestrictedSystem.getH2()[0]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*3+j*3+2)*RestrictedSystem.getH3()[0]-xpos;
	//				double yp = RestrictedSystem.getWrappedPos(id).y+RestrictedSystem.getCLNeighbours(i*nbNMax_res*3+j*3)*RestrictedSystem.getH1()[1]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*3+j*3+1)*RestrictedSystem.getH2()[1]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*3+j*3+2)*RestrictedSystem.getH3()[1]-ypos;
	//				double zp = RestrictedSystem.getWrappedPos(id).z+RestrictedSystem.getCLNeighbours(i*nbNMax_res*3+j*3)*RestrictedSystem.getH1()[2]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*3+j*3+1)*RestrictedSystem.getH2()[2]+RestrictedSystem.getCLNeighbours(i*nbNMax_res*3+j*3+2)*RestrictedSystem.getH3()[2]-zpos;
	//				Dist[(unsigned int) RestrictedSystem.getAux(gt_ind_res)[id] - 1].push_back(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)); 
	//			}
	//		}
	//		if( GrainId.size() < 2 ){
	//			cout << "Issue ! Increase cutoff radius !" << endl;
	//			GBAux[EdgeGrainId[i]] = 0.;
	//		}else{
	//			vector<double> FinalArray;
	//			for(unsigned int n=0;n<GrainId.size();n++){
	//				FinalArray.push_back(MT.min_vec(Dist[GrainId[n]-1]));
	//				FinalArray.push_back(GrainId[n]);
	//			}
	//			MT.sort(FinalArray,0,2,FinalArray);
	//			if( FinalArray[1] < FinalArray[3] ) GBAux[EdgeGrainId[i]] = FinalArray[1]+(FinalArray[3]/10.);
	//			else GBAux[EdgeGrainId[i]] = FinalArray[3]+(FinalArray[1]/10.);
	//		}
	//	} else GBAux[EdgeGrainId[i]] = 0.;
	//}

	cout << "Done !" << endl;
	cout << "Setting GBidentifier and printing system" << endl;

	MySystem.setAux(GBAux,"GBId");	
	MySystem.printSystem_aux(OutputFilename,"GBId grainID Struct");
	cout << "New system printed" << endl;
	
	delete[] GBAux;
	delete[] count_struct;
	delete[] RestrictedAtoms;
	//delete[] AuxStruct;
	//delete[] AuxGT;

	return 0;
}
