#ifndef ATOMHIC_ATOMICSYSTEM_H
#define ATOMHIC_ATOMICSYSTEM_H

#include "AtomHicExport.h"
#include "Crystal.h"
#include "MyStructs.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

/** Define an atomic system (i.e. a defined box containing atoms)
 */
class ATOMHIC_EXPORT AtomicSystem {
private:
	unsigned int nbAtom;
	unsigned int nbAtomType;
	unsigned int MaxAtomType=15;
	std::string *AtomType; // the different atom species present in the system
	unsigned int *AtomType_uint;
	double *AtomMass;
	Atom *AtomList; // List of the atom belonging to the system
	Position *WrappedPos;
	bool IsWrappedPos = false;
	unsigned int nbMaxN; // Maximum number of neighbours (computed as function of cuttof radius)
	unsigned int *Neighbours; // List of neigbors, structured as 2D array [ nbAtom x (nbMaxNeighbours+1) ] where the first value of each line corresponds to the number of neighbours for the ith atom
	int *CLNeighbours; // List of coefficient (NclX,Ncly,Nclz) applied to the atom position to be a neighbour [ nbAtom x (nbMaxNeighbours*3) ] where the first value of each line corresponds to the number of neighbours for the ith atom
	bool IsNeighbours = false;
	double *H1; // cell vectors of the system
	double *H2;
	double *H3;
	double *G1; // inverse cell vectors of the system
	double *G2;
	double *G3;
	double timestep; // timestep related to the atomic system (cfg lammps file)
	bool IsCharge;
	bool IsTilted;
	bool IsCrystalDefined;
	Crystal* MyCrystal;	
	std::vector<double*> Aux; // auxiliary atom properties
	bool IsSetAux = false;
	std::vector<std::string> Aux_name; // auxiliary atom properties
public:
	AtomicSystem();
	AtomicSystem(const std::string& filename); // construct AtomicSystem by reading file
	// getters
	int getNbAtom(){ return this->nbAtom; }
	Atom getAtom(const unsigned int Id){ return this->AtomList[Id]; }
	double* getAux(const unsigned int AuxId){ return this->Aux[AuxId]; }
	unsigned int* getNeighbours(){ return this->Neighbours; }
	unsigned int getNeighbours(const unsigned int Id){ return this->Neighbours[Id]; }
	int* getCLNeighbours(){ return this->CLNeighbours; }
	int getCLNeighbours(const unsigned int Id){ return this->CLNeighbours[Id]; }
	unsigned int getNbMaxN(){ return this->nbMaxN; }
	bool getIsNeighbours(){ return this->IsNeighbours; }
	double* getH1(){ return this->H1; }
	double* getH2(){ return this->H2; }
	double* getH3(){ return this->H3; }
	// methods
	void read_lmp_file(const std::string& filename);
	void read_cfg_file(const std::string& filename);
	void printSystem(const std::string& filename); // print a file containing atom type, charge, masse and position
	void print_lmp(const std::string& filename);
	void print_cfg(const std::string& filename);
	void printSystem_aux(const std::string& filename, const std::string& AuxId);
	void searchNeighbours(const double& rc);
	void computeWrap();
	void setAux(const double* aux, const std::string& AuxName);
	~AtomicSystem();
};

#endif
