//**********************************************************************************
//*   KMeans.h                                                                     *
//**********************************************************************************
//* This file contains the declaration of the KMeans class (herited from           *
//* MachineLearningModel) which employ the KMeans Algorythm to identify clusters   *
//* of Descriptors with a provided number of cluster.                              *
//* This class mainly contains:							   *
//*	- the KMeans++ initialization						   *
//*	- the KMeans algorithm							   *
//*	- the comptation of the weights, centroids and variances of each cluster   *
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
//*	- implement the Classify method (and remove the Data2Cluster array and     *
//* method)						   			   *
//*	- think about how to classify unknown descriptors                          *
//*	- use the distance function of Descriptors 				   *
//*	- 						                           *
//**********************************************************************************

#ifndef KMEANS_H
#define KMEANS_H

#include "AtomHicExport.h"
#include <string>
#include "MachineLearningModel.h"
#include "KMeansTools.h"
#include <Eigen/Dense>

class ATOMHIC_EXPORT KMeans : public MachineLearningModel {
private:
	// parameters of the KMeans for searching optimal number of cluster
	unsigned int optimal_nbClust;
	std::vector<std::vector<double>> optimal_centroids;

	// parameters of the KMeans ML model (filtered( and labelled)) the first index of the following vectors represent the filter value
	std::vector<unsigned int> nbClust;
	std::vector<std::vector<std::vector<double>>> centroids;

	std::vector<std::vector<unsigned int>> ClusterLabel;
	std::vector<long double*> AveLabelProb;
	std::vector<bool> IsLabelled;

	std::vector<KMeansTools*> MyKMeansTools;
	std::vector<bool> IsKMeansTools;

	// Variables for the labelling of the KMeans model
	double tolLabelSize = .7; // tolerance for warning message for repartition of descriptor in labels
	double tol2ndProb = .05; // tolerance for warning message for 2nd prob after labeling with the highest prob
	
	// FixedParameters
	double fac_elbow = 0.1; // reduction factor for considering that its is a real elbow
	unsigned int nb_sil_increase = 1;
	std::string after_elbow_choice = "Max";

public:
	// constructors
	KMeans();
	// base methods of MachineLearningModels
	void setDescriptors(Descriptors *D);
	void Classify();
	void TrainModel(unsigned int &nbClust, std::string filter_value="none"); // fit a KMeans with a given number of cluster and filter value 
	void TrainModel(unsigned int &nbClust_min, unsigned int &nbClust_max, const bool &softly_label=true); // train model (i.e. search optimal nb_clust and label the KMeans if the descriptors are labeled) for all filters
	void TrainModel(unsigned int &nbClust_min, unsigned int &nbClust_max, const std::string &filter_value, const bool &softly_label=true); // train model (i.e. search optimal nb_clust and label the KMeans if the descriptors are labeled) for a given filter
	void TrainModel(std::vector<unsigned int> &nbClust_min, std::vector<unsigned int> &nbClust_max, std::vector<std::string> &filter_value);// train model (i.e. search optimal nb_clust and label the KMeans if the descriptors are labeled) for different filters each having different values of nClustmin and max
	void fitOptimalKMeans(unsigned int &_nbClust_min, unsigned int &_nbClust_max, std::string namefile="");
	void Labelling(std::string filter_value);
	void Labelling();
private:
	void Classify(std::string filter_name);
	void ChangeFilterIndex();
	void readFixedParams();
	void appendOptimalKMeans(const std::string &filter_value="none");
	void NormalizeWeights(std::string filter_value);
	void setKMeansTools(std::string filter_value);
	void ComputeAveLabelProb(std::string filter_value);
	void PrintLabelling(std::string filter_value);

public:
	void PrintModelParams(std::string filename);
	//void PrintModelParams(std::string filename, std::vector<std::string> label_order); // specifiy the order of the label to print
	void PrintToDatabase(const std::string &name_of_database);
	void ReadModelParamFromDatabase(const std::string &name_of_database);
	void ReadProperties(std::vector<std::string> Properties);
	// getters
	unsigned int getNbClust(unsigned int &filter_value){ return nbClust[filter_value]; }
	// destructor
	~KMeans();
	
};

#endif
