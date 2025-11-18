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

#ifndef KMEANSTOOLS_H
#define KMEANSTOOLS_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

class ATOMHIC_EXPORT KMeansTools {
private:
	unsigned int _nbClust; // number of cluster in the GMMs
	// variables
	double _inertia;
	double _optimal_inertia;
	MatrixXd _centroids; // _centroids(k,d) = d component of the kth cluster of the kmeans model
	MatrixXd _centroids_old;
	MatrixXd _optimal_centroids;
	MatrixXd _V;

	unsigned int *_Data2Cluster = nullptr; // [i] = k cluster id of datapoint i
	unsigned int *_nbDat2Cluster = nullptr; // nbDat2Cluster[k] = number of points belonging to cluster k
	MathTools *MT;

	MatrixXd *_dataMat; // _dataMat(i,d) = d component of the ith descriptor
	unsigned int _nbDat;
	unsigned int _dim;
	unsigned int _dim2;

	// FixedParameters
	std::string FixedParam_Filename = "FixedParameters.ath";
	double tol_KMeans = 1e-5;
	unsigned int MaxIter_KMeans = 500;
	unsigned int nbInit = 10; // number of random initialization (we keep at the end the one with the highest likelihood)
	
	bool FixedSeed = false;
	unsigned int seed;
private:
	void InitializeKMeansVariables();
	void InitializeDataVariables();

public:
	// constructors
	KMeansTools(unsigned int &nbClust, MatrixXd *dataMat, unsigned int &nbDat, unsigned int &dim);
	KMeansTools(unsigned int &nbClust, MatrixXd *dataMat, unsigned int &nbDat, unsigned int &dim, std::vector<std::vector<double>> &centroids);
	KMeansTools(unsigned int &nbClust, unsigned int &dim, std::vector<std::vector<double>> &centroids);

	// methods
	void AffectData2Cluster();
	void KMeansPPInitialization();
	void ComputeFullVariances();
	void fit();
	void setSeed(unsigned int &seed);
	double ComputeSilhouette();

	// getters
	unsigned int getNbDat2Cluster(unsigned int &k){ return _nbDat2Cluster[k]; }
	MatrixXd getCentroids(){ return _optimal_centroids; }
	MatrixXd getV(){ return _V; }
	unsigned int *getData2Cluster(){ return _Data2Cluster; }

	// setters
	void setDataMat(MatrixXd *dataMat, unsigned int &nbDat);
	void setKMeansVariables(std::vector<std::vector<double>> &centroids);

	// read of fixed parameters
	void readFixedParams();
	void ReadProperties(std::vector<std::string> Properties);

	// destructor
	~KMeansTools();
};

#endif
