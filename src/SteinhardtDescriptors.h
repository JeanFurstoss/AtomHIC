#ifndef STEINHARDTDESCRIPTORS_H
#define STEINHARDTDESCRIPTORS_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "Descriptors.h"

class ATOMHIC_EXPORT SteinhardtDescriptors : public Descriptors {
private:
	// Properties of Steinhardt descriptors
	double rc; // cutoff radius for neighbor research
	int l_sph; // harmonic degree
	std::string SteinhardtStyle; // multi, mono
	std::string AverageStyle; // multi mono
	std::string mode; // OneL Full
	// variables for calculation
	unsigned int nbNMax; // maximum number of neighbours of the AtomicSystem (used to access to the different arrays)
	unsigned int *Malpha;// array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	double *Calpha; // normalization factor, only used for OneL mode (BondOriParam) 
	std::complex<double> *Qlm; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m (a bit different for OneL mode see function)
	unsigned int lsph2; // (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1; // l_sph*2+1;

	// Variables for calculation and printing
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog;
	unsigned int count_t;
public:
	// constructors
	SteinhardtDescriptors(AtomicSystem *_MySystem);
	SteinhardtDescriptors(AtomicSystem *_MySystem, std::vector<std::string> _Properties);
	// methods
	void readProperties(std::vector<std::string> _Properties);
	void ComputeDescriptors();
	void setProperties();
	void readFixedParams();
	void InitializeArrays();
	void ComputeBondOriParam();
	void ComputeSteinhardtParameters_Mono();
	void ComputeSteinhardtParameters_Multi();
	void AverageSteinhardtParameters_Mono();
	void AverageSteinhardtParameters_Multi();
	void printDescriptorsPropToDatabase(std::ofstream &writefile);
	// getters
	int get_l_sph(){ return l_sph; }
	double get_rc(){ return rc; }
	std::string getMode(){ return mode; }
	// destructor
	~SteinhardtDescriptors();
	
};

#endif
