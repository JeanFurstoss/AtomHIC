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

#ifndef GMMTOOLS_H
#define GMMTOOLS_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "KMeansTools.h"
#include <vector>
#include <chrono>
#include <Eigen/Dense>

using namespace Eigen;

class ATOMHIC_EXPORT GMMTools {
private:
	unsigned int _nbClust; // number of cluster in the GMMs
	
	// variables for EM 
	double _LogLikelihood, _optimal_LogLikelihood;
	VectorXd _weights; // weights[k] weight of cluster k
	MatrixXd _mu; // mu[k][d] d component of the mean of cluster k
	MatrixXd _V; // V[k][d1*dim+d2] d1,d2 component of variance of cluster k
	std::vector<LDLT<MatrixXd>> _V_ldlt;
	std::vector<long double> _det_V; // [k] determinant of variance of cluster k
	
	VectorXd _nbDat2Cluster; // [k] determinant of variance of cluster k
	MatrixXd _WeightedClusterLogProb; // [k][i] logarithm of the weighted gaussian probability of datapoint i in cluster k
	MatrixXd _Resp; // [k][i] logarithm of responsibilities for datapoint i in cluster k
	VectorXd _WeightedLogProb; // [i] logarithm of the sum of weighted gaussian probability of datapoint i for all clusters
	double logFacNorm; // log normalization the factor for gaussian prob (log(2pi^(-dim/2)))

	// optimal variables
	VectorXd _optimal_weights; // weights[k] weight of cluster k
	MatrixXd _optimal_mu; // mu[k][d] d component of the mean of cluster k
	MatrixXd _optimal_V; // V[k][d1*dim+d2] d1,d2 component of variance of cluster k
	std::vector<LDLT<MatrixXd>> _optimal_V_ldlt;
	std::vector<long double> _optimal_det_V; // [k] determinant of variance of cluster k
	
	MatrixXd regcovar;

	MatrixXd *_dataMat;

	unsigned int _nbDat;
	unsigned int _dim;

	KMeansTools *_MyKMeansTools;
	bool IsKMeans = false;
	bool FixedSeed = false;
	unsigned int seed;

	// FixedParameters
	std::string FixedParam_Filename = "FixedParameters.ath";
	double tol_Lkh_EM = 1e-4;
	unsigned int MaxIter_EM = 100;
	unsigned int nbInit = 10; // number of times the GMM is fitted (we keep at the end the one with the highest likelihood)
	std::string InitMethod = "KMEANSPP"; // number of random initialization (we keep at the end the one with the highest likelihood)
	std::vector<std::string> current_Properties;

	// times
	std::chrono::high_resolution_clock::time_point time_beg, time_end;
	double time_cluster_prob = 0.;
	double time_weighted_cluster_prob = 0.;
	double time_weighted_prob = 0.;
	double time_resp = 0.;
	double time_means = 0.;
	double time_variances = 0.;
	double time_nbdat = 0.;
	double time_kmeans = 0.;

private:
	void InitializeGMMVariables();
	void InitializeDataVariables();

public:
	// constructors
	GMMTools(unsigned int &nbClust, MatrixXd *dataMat, unsigned int &nbDat, unsigned int &dim);
	GMMTools(unsigned int &nbClust, MatrixXd *dataMat, unsigned int &nbDat, unsigned int &dim, std::vector<double> &weights, std::vector<std::vector<double>> &mu, std::vector<std::vector<double>> &V);
	GMMTools(unsigned int &nbClust, unsigned int &dim, std::vector<double> &weights, std::vector<std::vector<double>> &mu, std::vector<std::vector<double>> &V);
	
	// initialization methods
	void Initialize();
	void RandomInit();
	void InitFromKMeans();
	void setSeed(unsigned int &seed);
	
	// gaussian related methods
	void computeClusterProb();
	void ComputeWeightedClusterLogProb();
	void ComputeWeightedLogProb();
	double AverageWeightedLogProb();
	void ComputeResp();
	void ComputeWeights();
	void ComputeMeans();
	void ComputeVariances();
	void ComputeNbDat2Cluster();
	double ComputeBIC();
	void ComputeMLC(MatrixXd &MLC);

	// fitting methods
	void fit();
	double expectation(); // return LogLikelihood
	void maximization();
	void SaveOptimalRun();
	void SetOptimalRun();

	// getters
	double getWeights(unsigned int &k){ return _optimal_weights[k]; }
	double getMu(unsigned int &k, unsigned int &d){ return _optimal_mu(k,d); }
	double getV(unsigned int &k, unsigned int &d1, unsigned int &d2){ return _optimal_V(k*_dim+d1,d2); }

	// setters
	void setDataMat(MatrixXd *dataMat, unsigned int &nbDat);
	void setGMMVariables(std::vector<double> &weights, std::vector<std::vector<double>> &mu, std::vector<std::vector<double>> &V);

	// read of fixed parameters
	void readFixedParams();
	void ReadProperties(std::vector<std::string> Properties);

	~GMMTools();
};

#endif
