//**********************************************************************************
//*   SGMA_Alumine/main.cpp                                                        *
//**********************************************************************************
//* This file contains the implementation of the SGMA_Alumine test.                *
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
#include "SteinhardtDescriptors.h"
#include "GaussianMixtureModel.h"
#include "AtomicSystem.h"
#include <chrono>
#include <omp.h>

using namespace std;

int main(int argc, char *argv[])
{
		GaussianMixtureModel GMM;
		GMM.ReadModelParamFromDatabase("AluminaPhases_Furstoss_CompPhysCom");
		AtomicSystem MySystem("Input.cfg");
		string DescriptorName = "Steinhardt";
		string ftype="element";
		Descriptors MyDescriptors(&MySystem,DescriptorName,ftype);
		GMM.setDescriptors(&MyDescriptors);
		GMM.Classify();
		vector<string> label_order;
		label_order.push_back("Liquid");
		label_order.push_back("Gamma");
		label_order.push_back("FreeSurface");
		label_order.push_back("Alpha");
		GMM.setLabelOrder(label_order);
		MySystem.setAux_vec(GMM.getClassificator(),2,"Struct");
		MySystem.printSystem_aux("Output.cfg","Struct");
}
