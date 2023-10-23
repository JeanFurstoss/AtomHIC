// AtomHic library files
#include <AtomicSystem.h>
#include <Bicrystal.h>
#include <ComputeAuxiliary.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include "MyStructs.h"

using namespace std;

int main(int argc, char *argv[])
{
	string InputFilename, OutputFilename, CrystalType, outfilename, NormalDir;
	double SiO_d;
	double SecFac = 2.;
	if( argc == 3 ){
		InputFilename = argv[1];
		OutputFilename = argv[2];
	}else{
		cerr << "Usage: ./SearchLinkedTetrahedra InputFilename Dist_SiO_Tetrahedra xlo xhi ylo yhi zlo zhi CrystalType NormalDirection OutputFilename(optional)" << endl;
		return EXIT_FAILURE;
	}
	AtomicSystem MySystem(InputFilename);
	
	const unsigned int nbAt = MySystem.getNbAtom();

	unsigned int nbO(0), nbSi(0), nbMg(0);
	for(unsigned int i=0;i<nbAt;i++) {
		if( MySystem.getAtom(i).type_uint == 1 ) nbMg += 1;
		else if( MySystem.getAtom(i).type_uint == 2 ) nbSi += 1;
		else if( MySystem.getAtom(i).type_uint == 3 ) nbO += 1;
	}

	int NbExtraMg, NbExtraO;
	NbExtraMg = nbMg - nbSi*2;
	NbExtraO = nbO - nbSi*4;

	if( NbExtraMg < 0 || NbExtraO < 0 ) cout << "Issue" << endl;
	else cout << "Nb of extra Mg : " << NbExtraMg << ", nb extra O : " << NbExtraO << endl;

	MySystem.searchNeighbours(5.);
	const unsigned int nbNMax = MySystem.getNbMaxN();
	double dmin, xpos, ypos, zpos, xp, yp, zp, dtemp;
	unsigned int id;
	vector<double> MgList, OList;
	for(unsigned int i=0;i<nbAt;i++) {
		xpos = MySystem.getWrappedPos(i).x;
		ypos = MySystem.getWrappedPos(i).y;
		zpos = MySystem.getWrappedPos(i).z;
		if( MySystem.getAtom(i).type_uint == 1 ){
			dmin = 3.;
			for(unsigned int j=0;j<MySystem.getNeighbours(i*(nbNMax+1));j++){
				id = MySystem.getNeighbours(i*(nbNMax+1)+j+1);
				xp = MySystem.getWrappedPos(id).x+MySystem.getCLNeighbours(i*nbNMax*3+j*3)*MySystem.getH1()[0]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+1)*MySystem.getH2()[0]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+2)*MySystem.getH3()[0]-xpos;
				yp = MySystem.getWrappedPos(id).y+MySystem.getCLNeighbours(i*nbNMax*3+j*3)*MySystem.getH1()[1]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+1)*MySystem.getH2()[1]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+2)*MySystem.getH3()[1]-ypos;
				zp = MySystem.getWrappedPos(id).z+MySystem.getCLNeighbours(i*nbNMax*3+j*3)*MySystem.getH1()[2]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+1)*MySystem.getH2()[2]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+2)*MySystem.getH3()[2]-zpos;
				dtemp = sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)); 
				if( dtemp < dmin ) dmin = dtemp;
			}
			if( dmin < 2.5 ){
				MgList.push_back(i);
				MgList.push_back(dmin);
			}
		}else if( MySystem.getAtom(i).type_uint == 3 ){
			dmin = 3.;
			for(unsigned int j=0;j<MySystem.getNeighbours(i*(nbNMax+1));j++){
				id = MySystem.getNeighbours(i*(nbNMax+1)+j+1);
				xp = MySystem.getWrappedPos(id).x+MySystem.getCLNeighbours(i*nbNMax*3+j*3)*MySystem.getH1()[0]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+1)*MySystem.getH2()[0]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+2)*MySystem.getH3()[0]-xpos;
				yp = MySystem.getWrappedPos(id).y+MySystem.getCLNeighbours(i*nbNMax*3+j*3)*MySystem.getH1()[1]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+1)*MySystem.getH2()[1]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+2)*MySystem.getH3()[1]-ypos;
				zp = MySystem.getWrappedPos(id).z+MySystem.getCLNeighbours(i*nbNMax*3+j*3)*MySystem.getH1()[2]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+1)*MySystem.getH2()[2]+MySystem.getCLNeighbours(i*nbNMax*3+j*3+2)*MySystem.getH3()[2]-zpos;
				dtemp = sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)); 
				if( dtemp < dmin ) dmin = dtemp;
			}
			if( dmin < 2.5 ){
				OList.push_back(i);
				OList.push_back(dmin);
			}
		}
	}

	MathTools MT;
	MT.sort(MgList,1,2,MgList);
	MT.sort(OList,1,2,OList);
	vector<unsigned int> RmO(NbExtraO);
	vector<unsigned int> RmMg(NbExtraMg);
	for(unsigned int i=0;i<NbExtraO;i++) RmO[i] = OList[i*2];
	for(unsigned int i=0;i<NbExtraMg;i++) RmMg[i] = MgList[i*2];

	Atom *NewAtoms = new Atom[nbAt-NbExtraO-NbExtraMg];
	bool Store;
	unsigned int count = 0;
	for(unsigned int i=0;i<nbAt;i++){
		if( MySystem.getAtom(i).type_uint == 1 ){
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
				count++;
			}
		}else if( MySystem.getAtom(i).type_uint == 2 ){
			NewAtoms[count] = MySystem.getAtom(i);
			count++;
		}else if( MySystem.getAtom(i).type_uint == 3 ){
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
				count++;
			}
		}
	}
	cout << "Number of atom in the new system : " << count << endl;

	Crystal *MyCrystal = new Crystal("Forsterite");
	
	AtomicSystem NewSys(NewAtoms,nbAt-NbExtraO-NbExtraMg,MyCrystal,MySystem.getH1(),MySystem.getH2(),MySystem.getH3());
	
	NewSys.printSystem(OutputFilename);
	cout << "New system printed" << endl;
	return 0;
}
