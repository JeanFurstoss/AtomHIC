//**********************************************************************************
//*   Descriptors.h                                                                *
//**********************************************************************************
//* This file contains the declaration of the Descriptors class                    *
//* This class is not an abstract class because it could be constructed by reading *
//* the descriptor values from files (for dump files it uses the AtomicSystem      *
//* methods, and for other files it has its own readers)			   *
//* The descriptors are filtered (based on given properties) which can make the    *
//* access to the different arrays (e.g. descriptor values, neighbours) a bit      *
//* tricky (particularly true for MachineLearningModels). The FilterIndex array    *
//* should be used for this point.						   *
//* This class mainly contains:							   *
//*	- the properties of the descriptor (e.g. number of dimension, filtering    *
//* type) and the way of printing them 						   *
//*	- the readers and the methods for filtering the data			   *
//*	- distance function(s)							   *
//*	- a N dimensional neighbour research without periodic BC		   *
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
//*	- implement other type of distance functions (e.g. Mahalanobis)            *
//*	- test the neighbour research at high dimensions                           *
//*	- implement other types of descriptors (for the moment we only have        *
//* Steinhardt) such as SOAP, ACE..                                                *
//*	- construct subset of descriptors (i.e. classed by filter and/or label) to *
//* give them more efficiently (or more transparently) to ML models                *
//*	- add possibility to label descriptors from auxiliary property of dump file*
//**********************************************************************************

#ifndef DESCRIPTORS_H
#define DESCRIPTORS_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include <Eigen/Dense>

using namespace Eigen;

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
	unsigned int *FilterIndex = nullptr; // FilterIndex[f*nbDatMax+j] = i index in _Descriptor array of jth descriptor having filter f
	double *_Descriptors; // [i*dim+d] component d of descriptor i
	// TEST subarrays
	std::vector<MatrixXd*> DescriptorsSubarray; // [s](i,d) component d of descriptor i in descriptor subarray s
	std::vector<unsigned int> SubarraySize; // [s] number of descriptors in subarray s
	std::vector<std::string> SubarrayDescription; // [s] number of descriptors in subarray s
	std::vector<MatrixXd*> TestDataset; // [s](i,d) component d of descriptor i in descriptor test dataset s
	std::vector<unsigned int> TestDatasetSize; // [s] number of descriptors in subarray s
	std::vector<std::string> TestDatasetDescription; // [s] number of descriptors in subarray s
	std::vector<unsigned int*> CorresIndexSubarray; // [s][j] = i index in _Descriptors of jth descriptors in subarray s
	std::vector<unsigned int*> CorresIndexTestDataset; // [s][j] = i index in _Descriptors of jth descriptors in TestDataset s
	double _RatioTestDataset = 0.;
	bool subarray_defined = false;
	bool IsTestDataset = false;
	bool subarray_filter;
	bool subarray_label;
	unsigned int *DescriptorLabel = nullptr; // [i] = Label (and Filter, respectively) of ith descriptor in _Descriptor array
	unsigned int *DescriptorFilter = nullptr;

	// END TEST subarrays
	double *min_val = nullptr; // [f*dim+d] component d of the minimum value of descriptors with filter f
	double *max_val = nullptr; // [f*dim+d] component d of the maximum value of descriptors with filter f
	bool AreExtremums = false;

	bool AreDescriptorsMine = true;
	// label associated to each data, used for supervised learning (e.g. labeling of a GMM for structural analysis)
	unsigned int *Labels_uint = nullptr; // [f*nbDatMax+i] = l with l the unsigned int label of descriptor i with filter f
	std::vector<std::string> Labels; // Labels[l] = name of the label l
	unsigned int *LabelsSize = nullptr; // [f*Labels.size()+l] = number of descriptor having label l in filter f
	unsigned int *LabelIndex = nullptr; // [f*Labels.size()*nbDatMaxLabel+l*nbDatMaxLabel+j] = i index in _Descriptor array of j_th descriptor having label l in filter f
	unsigned int nbLabels = 1; // number of labels
	unsigned int nbDatMaxLabel; // maximum number of data in all labels (used to acces LabelIndex array)
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
	void ConstructLabelIndexArray();
	void ConstructDescriptorFromAtomicPosition();
	void searchNeighbours(double rc, std::string filter_val, std::string distFunc="Euclidian");
	void ComputeExtremums();
	void SquaredDistance(unsigned int id_1, unsigned int id_2, unsigned int filter, double &squared_dist);// For the moment Euclidian distance (to add other distance and rename)
	void readProperties(std::vector<std::string> _Properties);
	// TEST subarrays
	void constructSubarrays(bool filter=false, bool label=false);
	void constructSubarrays(double &RatioTestDataset, bool filter=false, bool label=true);
	MatrixXd *getSubarray(unsigned int &subarray_nbdat, std::string filter_name="none", std::string label_name="none");
	MatrixXd *getTestDataset(unsigned int &subarray_nbdat, std::string filter_name="none", std::string label_name="none");
	unsigned int *getCorresIndexSubarray(unsigned int &subarray_nbdat, std::string filter_name="none", std::string label_name="none");
	unsigned int *getCorresIndexTestDataset(unsigned int &subarray_nbdat, std::string filter_name="none", std::string label_name="none");
	unsigned int getLabels_uint(unsigned int &i){ return DescriptorLabel[i]; }
	unsigned int current_filter(const std::string &filter_name);
	unsigned int current_label(const std::string &label_name);
	// END TEST subarrays
	// getters
	double *getDescriptors(){ return _Descriptors; }
	unsigned int *getLabels_uint(){ return Labels_uint; }
	unsigned int getLabelsSize(const unsigned int &f, unsigned int &l){ return LabelsSize[f*Labels.size()+l]; }
	unsigned int getNbLabels(){ return Labels.size(); }
	std::string getLabels(const unsigned int &l){ return Labels[l]; }
	unsigned int getNbDatMaxLabel(){ return nbDatMaxLabel; }
	unsigned int getLabelIndex(const unsigned int &f){ return LabelIndex[f]; }
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
