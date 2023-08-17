// AtomHic library files
#include <AtomicSystem.h>
#include <Bicrystal.h>
#include <ComputeAuxiliary.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"

using namespace std;

int main(int argc, char *argv[])
{
	string InputFilename, OutputFilename, CrystalType, outfilename, NormalDir;
	double SiO_d;
	double SecFac = 2.;
	if( argc == 5 || argc == 6 ){
		InputFilename = argv[1];
		istringstream iss_d(argv[2]);
		CrystalType = argv[3];
		NormalDir = argv[4];
		iss_d >> SiO_d;
	}else{
		cerr << "Usage: ./SearchLinkedTetrahedra InputFilename Dist_SiO_Tetrahedra xlo xhi ylo yhi zlo zhi CrystalType NormalDirection OutputFilename(optional)" << endl;
		return EXIT_FAILURE;
	}
	if( argc == 6 ){
		OutputFilename = argv[5];
	}
	Bicrystal MySystem1(InputFilename, NormalDir, CrystalType);

	double z_lo = MySystem1.getGBPos1() - (MySystem1.getGBwidth1()/2.);
	double z_hi = MySystem1.getGBPos1() + (MySystem1.getGBwidth1()/2.);

	AtomicSystem MySystem(InputFilename);
	SiO_d *= SecFac;
	MySystem.searchNeighbours(SiO_d);
	SiO_d /= SecFac;
	const unsigned int nbNMax = MySystem.getNbMaxN();
	const unsigned int nbAt = MySystem.getNbAtom();
	vector<unsigned int> At_index;
	At_index = MySystem.selectAtomInBox(-1e9,1e9,-1e9,1e9,z_lo,z_hi);
	unsigned int ind_O, ind_Si;
	for(unsigned int i=0;i<MySystem.getNbAtomType();i++){
		if( MySystem.getAtomType(i) == "Si" ) ind_Si = i+1;
		else if( MySystem.getAtomType(i) == "O" ) ind_O = i+1;
	}
	unsigned int id, nbO, NbLinked, nbSi;
	double xpos, ypos, zpos, xp, yp, zp, dist;
	unsigned int *nbT = new unsigned int[nbAt];
	vector<unsigned int> nbT_O;
	nbO = 0;
	NbLinked = 0;
	nbSi = 0;
	bool inside;
	for(unsigned int n=0;n<nbAt;n++){
		inside = false;
		for(unsigned int i=0;i<At_index.size();i++){
			if( n == At_index[i] ){
				if( MySystem.getAtom(At_index[i]).type_uint == ind_O ){
					nbO += 1;
					nbT_O.push_back(0);
					xpos = MySystem.getWrappedPos(At_index[i]).x;
        	        		ypos = MySystem.getWrappedPos(At_index[i]).y;
        	        		zpos = MySystem.getWrappedPos(At_index[i]).z;
					for(unsigned int j=0;j<MySystem.getNeighbours(At_index[i]*(nbNMax+1));j++){
						id = MySystem.getNeighbours(At_index[i]*(nbNMax+1)+j+1);
						if( MySystem.getAtom(id).type_uint == ind_Si ){
		        	                        xp = MySystem.getWrappedPos(id).x+MySystem.getCLNeighbours(At_index[i]*nbNMax*3+j*3)*MySystem.getH1()[0]+MySystem.getCLNeighbours(At_index[i]*nbNMax*3+j*3+1)*MySystem.getH2()[0]+MySystem.getCLNeighbours(At_index[i]*nbNMax*3+j*3+2)*MySystem.getH3()[0]-xpos;
		        	                        yp = MySystem.getWrappedPos(id).y+MySystem.getCLNeighbours(At_index[i]*nbNMax*3+j*3)*MySystem.getH1()[1]+MySystem.getCLNeighbours(At_index[i]*nbNMax*3+j*3+1)*MySystem.getH2()[1]+MySystem.getCLNeighbours(At_index[i]*nbNMax*3+j*3+2)*MySystem.getH3()[1]-ypos;
        			                        zp = MySystem.getWrappedPos(id).z+MySystem.getCLNeighbours(At_index[i]*nbNMax*3+j*3)*MySystem.getH1()[2]+MySystem.getCLNeighbours(At_index[i]*nbNMax*3+j*3+1)*MySystem.getH2()[2]+MySystem.getCLNeighbours(At_index[i]*nbNMax*3+j*3+2)*MySystem.getH3()[2]-zpos;
							dist = sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.));
							if( dist < SiO_d ) nbT_O[nbO-1] += 1;
						}
					}
					if( nbT_O[nbO-1] == 2 ){
						NbLinked += 1;
						nbT[At_index[i]] = 2;
					}else nbT[At_index[i]] = 1;
				}
				else if( MySystem.getAtom(At_index[i]).type_uint == ind_Si ){
					nbSi += 1;
					nbT[At_index[i]] = 1;
				}
				else nbT[At_index[i]] = 1;
				inside = true;
				At_index.erase(At_index.begin()+i);
				break;
			}
		}
		if( !inside ) nbT[n] = 0;
	}

	if( argc == 6 ){
		//MySystem.setAux(nbT,"LinkedTetrahedra");
		//MySystem.printSystem_aux(OutputFilename,"LinkedTetrahedra");
		ofstream writefile(OutputFilename);
		writefile << MySystem1.getGBPos1() << " " << MySystem1.getGBwidth1() << " " << MySystem1.getExcessVol() << " " << NbLinked << " " << nbSi;
		writefile.close();
	}

	cout <<	"Number of linked tetrahedra : " << NbLinked << " over a total number of " << nbSi << " tetrahedra" << endl;

	return 0;
}
