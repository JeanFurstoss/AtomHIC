#ifndef ATOMICSYSTEM_H
#define ATOMICSYSTEM_H

#include "AtomHicExport.h"
#include "MathTools.h"
#include "MyStructs.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "Crystal.h"

class ATOMHIC_EXPORT AtomicSystem {
protected:
	unsigned int nbAtom;
	unsigned int nbAtomType;
	unsigned int MaxAtomType=15;
	std::string *AtomType; // the different atom species present in the system
	double *AtomMass;
	double *AtomCharge;
	Atom *AtomList; // List of the atom belonging to the system
	bool IsAtomListMine = true;
	Position *WrappedPos;
	bool IsWrappedPos = false;
	unsigned int nbMaxN; // Maximum number of neighbours (computed as function of cuttof radius)
	unsigned int *Neighbours; // List of neigbors, structured as 2D array [ nbAtom x (nbMaxNeighbours+1) ] where the first value of each line corresponds to the number of neighbours for the ith atom, i.e. Neighbours[i*(nbMaxN+1)] = nb of neighbour of atom i, Neighbours[i*(nbMaxN+1)+j+1] = id of the jth neighbour of atom i
	int *CLNeighbours; // List of coefficient (NclX,Ncly,Nclz) applied to the atom position to be a neighbour [ nbAtom x (nbMaxNeighbours*3) ] where the first value of each line corresponds to the number of neighbours for the ith atom, i.e. CLNeighbours[i*nbMaxN*3+j*3] = integer applied to H1 to the jth neighbour of atom i for being its neighbour, CLNeighbours[i*nbMaxN*3+j*3+1] = integer applied to H2 to the jth neighbour of atom i for being its neighbour, CLNeighbours[i*nbMaxN*3+j*3+2] = integer applied to H3 to the jth neighbour of atom i for being its neighbour
	bool IsNeighbours = false;
	std::vector<int> *NotSepTag; // array containing the atom which should not be separated: NotSepTag[i][0] = -1 => the atom i is stored because of the NotSep list, 0 => this atom is not related to the NotSep list, n => the atom has stored n neighbour because of the NotSep list, in the latter case: NotSepTag[i][j] return the index of the (j-1)th atom which has been stored with the ith one due to the NotSepList
	bool IsNotSepTag = false;
	double *H1; // cell vectors of the system
	double *H2;
	double *H3;
	bool IsCellVecMine = true;
	double *G1; // inverse cell vectors of the system
	double *G2;
	double *G3;
	bool IsG = false;
	double timestep; // timestep related to the atomic system (cfg lammps file)
	bool IsCharge;
	bool IsTilted;
	bool IsCrystalDefined = false;
	bool IsCrystalMine = false;
	bool FilenameConstructed = false;
	Crystal *_MyCrystal;	
	std::vector<double*> Aux; // auxiliary atom properties
	std::vector<unsigned int> Aux_size; // size of auxiliary atom properties
	bool IsSetAux = false;
	std::vector<std::string> Aux_name; // auxiliary atom properties
	MathTools *MT;
	std::vector<double*> density_prof; // density_prof[i][j*2] = density of auxiliary property i at sampled point j, density_prof[i][j*2+1] = coordinate (in density_name[i*2+1] direction) of the sample point j
	std::vector<int> density_nbPts; // density_nbPts[i] = number of sampled point
	std::vector<std::string*> density_name; // auxiliary atom properties => density_name[i*2] = auxiliary property used to compute density, density_name[i*2+1] = direction along which the density has been computed
	std::string File_Heading; // head of lmp printed file
	// Parameters read from Fixed_Parameter.dat file
	std::string FixedParam_Filename = "Fixed_Parameters.dat";
	double r_cut_n;
	double current_rc_neigh;
	int l_sph_st;
public:
	AtomicSystem(){};
	AtomicSystem(Crystal *_MyCrystal, double xhi, double yhi, double zhi, std::vector<int> cl_box); // construct atomic system from crystal and cell size
	AtomicSystem(const std::string& filename); // construct AtomicSystem by reading file
	AtomicSystem(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3); // construct AtomicSystem giving AtomList and cell vectors 
	void FilenameConstructor(const std::string& filename);
	// getters
	unsigned int getNbAtom(){ return this->nbAtom; }
	Atom getAtom(const unsigned int Id){ return this->AtomList[Id]; }
	double* getAux(const unsigned int AuxId){ return this->Aux[AuxId]; }
	double get_current_rc(){ return this->current_rc_neigh; }
	Position getWrappedPos(const unsigned int AuxId){ if( !IsWrappedPos ) computeWrap(); return this->WrappedPos[AuxId]; }
	unsigned int* getNeighbours(){ return this->Neighbours; }
	unsigned int getNeighbours(const unsigned int Id){ return this->Neighbours[Id]; }
	int* getCLNeighbours(){ return this->CLNeighbours; }
	std::vector<int>* getNotSepTag(){ return this->NotSepTag; }
	int getCLNeighbours(const unsigned int Id){ return this->CLNeighbours[Id]; }
	unsigned int getNbMaxN(){ return this->nbMaxN; }
	bool getIsNeighbours(){ return this->IsNeighbours; }
	bool getIsCrystalDefined(){ return this->IsCrystalDefined; }
	Crystal* getCrystal(){ return this->_MyCrystal; }
	double* getH1(){ return this->H1; }
	double* getH2(){ return this->H2; }
	double* getH3(){ return this->H3; }
	double get_rcut(){ return this->r_cut_n; }
	int get_lsph(){ return this->l_sph_st; }
	// setters
	void setAux(const double* aux, const std::string& AuxName);
	void setAux_vec(const double* aux, const unsigned int size, const std::string& AuxName);
	void modifyAux_vec(const double* aux, const std::string& AuxName);
	void setAux(const unsigned int* aux, const std::string& AuxName);
	void setAux(const int* aux, const std::string& AuxName);
	void setCrystal(Crystal* MyCrystal);
	void setCrystal(const std::string& CrystalName);
	void set_File_Heading(const std::string& Heading){ this->File_Heading = Heading; }
	// methods
	void read_params_atsys();
	void computeInverseCellVec();
	void read_lmp_file(const std::string& filename);
	void read_cfg_file(const std::string& filename);
	void printSystem(const std::string& filename); // print a file containing atom type, charge, masse and position
	void print_lmp(const std::string& filename);
	void print_cfg(const std::string& filename);
	void printSystem_aux(const std::string& filename, const std::string& AuxId);
	void searchNeighbours(const double& rc);
	void computeWrap();
	unsigned int Compute1dDensity(std::string auxname, std::string dir, double sigma, unsigned int nbPts); // compute and store the 1D density profile of a given auxiliary property, the metho return the index of the given density
	void Print1dDensity(std::string filename, std::string auxname);
	std::vector<unsigned int> selectAtomInBox(const double x_lo,const double x_hi,const double y_lo,const double y_hi,const double z_lo,const double z_hi);
	~AtomicSystem();
};

#endif
