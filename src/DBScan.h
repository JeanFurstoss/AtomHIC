//**********************************************************************************
//*   DBScan.h                                                                     *
//**********************************************************************************
//* This file contains the declaration of the DBScan class (herited from           *
//* MachineLearningModel) which employ the Density Based Scan Algorythm to         *
//* identify clusters of Descriptors.                                               *
//* The returned Classificator contains the cluster id and the status of the given *
//* Descriptor.
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
//*	- think about how to separate the TrainModel and the Classify methods      *
//*	- implement a way for minPts automatic calculation or using the read value *
//* from FixedParameters                                                           *
//*	- think about how to classify unknown descriptors                          *
//*	- 						                           *
//**********************************************************************************

#ifndef DBSCAN_H
#define DBSCAN_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "MachineLearningModel.h"

class ATOMHIC_EXPORT DBScan : public MachineLearningModel {
private:
	bool IsDescriptor = false;
	unsigned int *nbClust; // [f]
	double *centroids; // centroids[k*dim*nbFilter+d*nbFilter+f] d component of the centroid of cluster k with filter f
	AtomicSystem *_MySystem;

	// Fixed parameters
	double eps = 3.;
	unsigned int minPts = 5;

public:
	// constructors
	DBScan();
	// methods
	void readFixedParams();
	void setDescriptors(Descriptors *D);
	void TrainModel(std::string filter_name);
	void Classify(){};
	void ReadProperties(std::vector<std::string> Properties);
	//void ReadModelParamFromDatabase();
	// destructor
	~DBScan();
	
};

#endif
