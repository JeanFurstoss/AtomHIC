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
	MatrixXd _centroids;
	MatrixXd _centroids_old;
	MatrixXd _optimal_centroids;
	MatrixXd _V;

	unsigned int *_Data2Cluster; // [i] = k cluster id of datapoint i
	unsigned int *_nbDat2Cluster; // nbDat2Cluster[k] = number of points belonging to cluster k
	MathTools *MT;

	MatrixXd *_dataMat;
	unsigned int _nbDat;
	unsigned int _dim;
	unsigned int _dim2;

	// FixedParameters
	std::string FixedParam_Filename = "Fixed_Parameters.dat";
	double tol_KMeans;
	unsigned int MaxIter_KMeans;
	unsigned int nbInit; // number of random initialization (we keep at the end the one with the highest likelihood)


public:
	// constructors
	KMeansTools(unsigned int &nbClust, MatrixXd *dataMat, unsigned int &nbDat, unsigned int &dim);
	// methods
	//double SquareEuclidianDistance2Centroid(const unsigned int &i, const unsigned int &ClusterIndex);
	void AffectData2Cluster();
	void KMeansPPInitialization();
	void ComputeFullVariances();
	void fit();

	// getters
	unsigned int getNbDat2Cluster(unsigned int &k){ return _nbDat2Cluster[k]; }
	//double *getCentroids(unsigned int &k){ return _optimal_centroids[k]; }
	MatrixXd getCentroids(){ return _optimal_centroids; }
	MatrixXd getV(){ return _V; }
	//long double *getV(unsigned int &k){ return _V[k]; }
	unsigned int *getData2Cluster(){ return _Data2Cluster; }

	// read of fixed parameters
	void readFixedParams();
	void ReadProperties(std::vector<std::string> Properties);

	// destructor
	~KMeansTools();
};

#endif
