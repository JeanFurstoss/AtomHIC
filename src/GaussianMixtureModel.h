//**********************************************************************************
//*   GaussianMixtureModel.h                                                       *
//**********************************************************************************
//* This file contains the declaration of the GaussianMixtureModel class (herited  *
//* from MachineLearningModel) which employ the Expectation Maximization algorythm *
//* to fit a mixture of gaussian distribution on descriptors.			   *
//* This class mainly contains:							   *
//*	- initialization methods (can be set from FixedParameters) 		   *
//* 	- a method for labeling the fitted GMM if the descriptors are labelled     *
//* (used for the Steinhardt Gaussian Mixture Analysis (SGMA) presented in         *
//* (Furstoss et al., 2025, Comp. Phys. Comm.)					   *
//* 	- readers and printers of labeled GMM in the Machine Learning database     *
//*	- methods for determining the optimal number of clusters		   *
//*	- classification of unknown descriptors based on Maximum Likelihood 	   *
//* Classifier			   						   *
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
//*	-                                                                          *
//**********************************************************************************

#ifndef GAUSSIANMIXTUREMODEL_H
#define GAUSSIANMIXTUREMODEL_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "MachineLearningModel.h"
#include "KMeans.h"
#include "GMMTools.h"

class ATOMHIC_EXPORT GaussianMixtureModel : public MachineLearningModel {
private:
	// parameters of the GMM for searching optimal number of cluster
	unsigned int optimal_nbClust;
	std::vector<double> optimal_weights;
	std::vector<std::vector<double>> optimal_mu;
	std::vector<std::vector<double>> optimal_V;

	// parameters of the GMM ML model (filtered( and labelled)) the first index of the following vectors represent the filter value
	std::vector<unsigned int> nbClust;
	std::vector<std::vector<double>> weights;
	std::vector<std::vector<std::vector<double>>> mu;
	std::vector<std::vector<std::vector<double>>> V;
	
	std::vector<std::vector<unsigned int>> ClusterLabel;
	std::vector<long double*> AveLabelProb;
	std::vector<bool> IsLabelled;

	std::vector<GMMTools*> MyGMMTools;
	std::vector<bool> IsGMMTools;

	// Variables for the labelling of the GMM
	double tolLabelSize = .7; // tolerance for warning message for repartition of descriptor in labels
	double tol2ndProb = .05; // tolerance for warning message for 2nd prob after labeling with the highest prob

	// FixedParameters
	double fac_elbow = 0.1; // reduction factor for considering that its is a real elbow
	unsigned int nb_bic_increase = 1;
	std::string after_elbow_choice = "Max";

public:
	// constructors
	GaussianMixtureModel();
	// base methods of MachineLearningModels
	void setDescriptors(Descriptors *D);
	void Classify();
	void TrainModel(unsigned int &nbClust, std::string filter_value="none"); // fit a GMM with a given number of cluster and filter value 
	void TrainModel(unsigned int &nbClust_min, unsigned int &nbClust_max, const bool &softly_label=true); // train model (i.e. search optimal nb_clust and label the GMM if the descriptors are labeled) for all filters
	void TrainModel(unsigned int &nbClust_min, unsigned int &nbClust_max, const std::string &filter_value, const bool &softly_label=true); // train model (i.e. search optimal nb_clust and label the GMM if the descriptors are labeled) for a given filter
	void TrainModel(std::vector<unsigned int> &nbClust_min, std::vector<unsigned int> &nbClust_max, std::vector<std::string> &filter_value);// train model (i.e. search optimal nb_clust and label the GMM if the descriptors are labeled) for different filters each having different values of nClustmin and max
	void fitOptimalGMM(unsigned int &_nbClust_min, unsigned int &_nbClust_max, std::string namefile="");
	void Labelling(std::string filter_value);
	void Labelling();
	// specific methods of GMM
private:
	void Classify(std::string filter_name);
	void ChangeFilterIndex();
	void readFixedParams();
	void appendOptimalGMM(const std::string &filter_value="none");
	unsigned int getCurrentFIndex(std::string filter_value);
	void NormalizeWeights(std::string filter_value);
	void setGMMTools(std::string filter_value);
	void ComputeAveLabelProb(std::string filter_value);
	void PrintLabelling(std::string filter_value);

public:
	void PrintModelParams(std::string filename);
	void PrintModelParams(std::string filename, std::vector<std::string> label_order); // specifiy the order of the label to print
	void PrintToDatabase(const std::string &name_of_database);
	void ReadModelParamFromDatabase(const std::string &name_of_database);
	void ReadProperties(std::vector<std::string> Properties);
	// getters
	std::vector<std::string> getAvailableDatabases();
	unsigned int getNbClust(unsigned int &filter_value){ return nbClust[filter_value]; }
	//long double *getCov(){ return V; }
	//double *getMu(){ return mu; }
	// destructor
	~GaussianMixtureModel();
};

#endif
