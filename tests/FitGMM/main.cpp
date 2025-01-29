//**********************************************************************************
//*   FitGMM/main.cpp                                                              *
//**********************************************************************************
//* This file contains the implementation of the FitGMM test.                      *
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

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "Descriptors.h"
#include "GaussianMixtureModel.h"

using namespace std;

int main(int argc, char *argv[])
{
	Descriptors MyDescriptors("ForTest");
	GaussianMixtureModel GMM;
	vector<string> Properties;
	Properties.push_back("GMM_TOL_LKH_EM 1e-5");
	Properties.push_back("GMM_MAX_ITER_EM 500");
	Properties.push_back("GMM_NB_MAX_CLUSTER 50");
	Properties.push_back("GMM_ELBOW_FACTOR 0.1");
	Properties.push_back("GMM_NB_BIC_INCREASE_FOR_MIN 2");
	Properties.push_back("GMM_NB_INIT 1");
	Properties.push_back("GMM_INIT_METHOD KMEANS");
	GMM.ReadProperties(Properties);
	vector<string> KMProp;
	KMProp.push_back("KMEANS_NB_MAX_CLUSTER 100");
	KMProp.push_back("KMEANS_TOL 1e-5");
	KMProp.push_back("KMEANS_MAX_ITER 1000");
	KMProp.push_back("KMEANS_NB_INIT 10");
	GMM.SetKMeansProperties(KMProp);
	GMM.setDescriptors(&MyDescriptors);
	unsigned int nmin= 3;
	unsigned int nmax= 15;
	GMM.fitOptimalGMM(nmin,nmax);
	GMM.Labelling();
	vector<string> label_order;
	string lab="lab";
	for(unsigned int l=0;l<5;l++){
		string num = to_string(l);
		string result = lab + num;
		label_order.push_back("Lab"+to_string(l+1));
	}
	GMM.PrintModelParams("output.dat", label_order);
}
