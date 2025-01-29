//**********************************************************************************
//*   MachineLearningModel.h                                                       *
//**********************************************************************************
//* This file contains the declaration of the MachineLearningModel class which is  *
//* basis class of the ML models of AtomHic 			  	           *
//* This class has:								   *
//*	- some descriptors (and should care about the filter of descriptors)	   *
//*	- classificators (or estimators)					   *
//*	- readers and writters of the model in the ML database of AtomHic	   *
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
//*	- 						                           *
//**********************************************************************************

#ifndef MACHINELEARNINGMODEL_H
#define MACHINELEARNINGMODEL_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "Descriptors.h"

class ATOMHIC_EXPORT MachineLearningModel {
protected:
	std::string name;
	Descriptors *_MyDescriptors;
	bool IsDescriptor = false;
	MathTools *MT;

	double *buffer_vec_1_dim;
	double *buffer_vec_2_dim;
	
	unsigned int dim;
	unsigned int dim2;
	unsigned int nbDatMax;
	unsigned int *nbDat; // [f] number of data with filter f
	unsigned int nbFilter; // number of filter (e.g. if filtered by element, the number of element)
	// Variables for labelling the model
	bool IsLabelled = false;
	unsigned int nbLabel;
	
	// Variables for reading the database parameters
	std::vector<std::string> Labels; // Labels[l] = name of the label l
	std::string FilteringType;
	std::vector<std::string> FilterValue; // FilterValue[f] = name of the filter f
	std::string DescriptorName; // to be used for initializing descriptors
	std::vector<std::string> DescriptorProperties; // to be used for initializing descriptors
	unsigned int *FilterIndexToModify;
	bool IsFilterIndexModified = false;
	bool IsRead = false;

	// Variables for classification
	double *Classificator; // depend on the ML model but generally for classification, [n*2] = index of the label with highest probability for descriptors n (index in _Descriptor array), [n*2+1] = probability
	bool IsClassified = false;	
	
	std::string FixedParam_Filename = "Fixed_Parameters.dat";

public:
	// constructors
	MachineLearningModel();
	// methods
	void setDescriptors(Descriptors *D);
	void TrainModel();
	std::string getMLDatabasePath();
	std::string getDatabasePath(const std::string &name_of_database);
	void ChangeFilterIndex();
	//void ReadModelParamFromDatabase(); // TODO implement base fonction here
	void Classify();
	void PrintClassifiedData(std::string filename);
	void readFixedParams();
	// getters
	double *getClassificator(){ return Classificator; }
	std::string getFilteringType(){ return FilteringType; }
	std::vector<std::string> getDescriptorProperties(){ return DescriptorProperties; }
	// destructor
	~MachineLearningModel();
	
};

#endif
