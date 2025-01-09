#ifndef DESCRIPTORS_H
#define DESCRIPTORS_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>

class ATOMHIC_EXPORT Descriptors {
protected:
	std::string name = "none";
	std::vector<std::string> Properties; // Descriptor properties could be set from DescriptorProperties.ath file when reading labelled data or from MachineLearning database reading or when computed by the setProperties() method
	AtomicSystem *_MySystem;
	MathTools *MT;
	unsigned int dim; // dimension of descriptors
	std::vector<unsigned int> nbDat; // nbDat[f] number of data with filter f
	unsigned int nbDatMax; // maximum number of data in all filter
	unsigned int nbDatTot;
	unsigned int *FilterIndex; // FilterIndex[f*nbDatMax+j] = i index in _Descriptor array of jth descriptor having filter f
	double *_Descriptors; // [i*dim+d] component d of descriptor i
	double *min_val; // [f*dim+d] component d of the minimum value of descriptors with filter f
	double *max_val; // [f*dim+d] component d of the maximum value of descriptors with filter f
	bool AreExtremums = false;

	bool AreDescriptorsMine = true;
	// label associated to each data, used for supervised learning (e.g. labeling of a GMM for structural analysis)
	unsigned int *Labels_uint; // [f*nbDatMax+i] = l with l the unsigned int label of descriptor i with filter f
	std::vector<std::string> Labels; // Labels[l] = name of the label l
	unsigned int *LabelsSize; // [f*Labels.size()+l] = number of descriptor having label l in filter f
	// Neighbours of descriptors
	std::vector<unsigned int*> Neighbours; // List of neigbors, structured as 2D array for each filter f such that Neighbours[f][i*(nbMaxN[f]+1)] = nb of neighbour of descriptor i, Neighbours[f][i*(nbMaxN[f]+1)+j+1] = id of the jth neighbour of descriptor i
	std::vector<unsigned int> nbMaxN;
	std::vector<double> rc_neighbours;
	std::vector<std::string> Neigh_FilterValue;
	bool IsNeighbours = false;

	// used to filtering the data (none by default but could be by element for instance i.e. the descriptor list is separate for the different types)
	unsigned int nbFilter = 0; //
	std::string FilteringType = "none";
	std::vector<std::string> FilterValue; // FilterValue[f] = name of the f filter
	unsigned int nbMaxFilter = 100; // only used for warning
	
	std::string FixedParam_Filename = "Fixed_Parameters.dat";

public:
	// constructors
	Descriptors(AtomicSystem *_MySystem);
	Descriptors(AtomicSystem *_MySystem, std::vector<std::string> _Properties);
	Descriptors(AtomicSystem *_MySystem, std::string DescriptorName, std::string _FilteringType); // Decriptor name could be an auxiliary property of the AtomicSystem or Position of atoms
	Descriptors(const std::string& FilenameOrDir); // constructor from file reading 
	Descriptors(const std::string& FilenameOrDir, const std::string& DescriptorName); // constructor from file reading 
	Descriptors(double *Descriptors, unsigned int nbdat, unsigned int dim);
	// methods
	void printDescriptorsPropToDatabase(std::ofstream &writefile);
	void setProperties(){};
	void ComputeDescriptors(){};
	void readFixedParams();
	void ConstructFilterIndexArray(AtomicSystem *_MySystem);
	void ConstructDescriptorFromAtomicPosition();
	void searchNeighbours(double rc, std::string filter_val, std::string distFunc="Euclidian");
	void ComputeExtremums();
	void SquaredDistance(unsigned int id_1, unsigned int id_2, unsigned int filter, double &squared_dist);// For the moment Euclidian distance (to add other distance and rename)
	void readProperties(std::vector<std::string> _Properties);
	// getters
	double *getDescriptors(){ return _Descriptors; }
	unsigned int *getLabels_uint(){ return Labels_uint; }
	unsigned int getLabelsSize(const unsigned int &f, unsigned int &l){ return LabelsSize[f*Labels.size()+l]; }
	unsigned int getNbLabels(){ return Labels.size(); }
	std::string getLabels(const unsigned int &l){ return Labels[l]; }
	unsigned int getDim(){ return dim; }
	unsigned int getNbDat(const unsigned int &f){ return nbDat[f]; }
	unsigned int getNbDatMax(){ return nbDatMax; }
	unsigned int getNbFilter(){ return nbFilter; }
	unsigned int getFilterIndex(const unsigned int &f){ return FilterIndex[f]; }
	std::string getFilteringType(){ return FilteringType; }
	std::string getFilterValue(const unsigned int &f){ return FilterValue[f]; }
	unsigned int* getNeighbours(unsigned int f){ return Neighbours[f]; }
	unsigned int getNbMaxNAndFilter(std::string filter_name, unsigned int &f);
	bool getIsNeighbours(std::string filter_name);
	double get_current_rc(std::string filter_name);
	// destructor
	~Descriptors();
	
};

#endif
