//**********************************************************************************
//*   DBScanClustering/main.cpp                                                    *
//**********************************************************************************
//* This file contains the implementation of the DBScanClustering executable.      *
//* It allows to make a cluster analysis of atomic position using the density	   *
//* based clustering algorythm							   *
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
//*	- add possibility of filtering ions					   *
//*	- generalize to the clustering of descriptors and not only on position	   *
//**********************************************************************************

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "Descriptors.h"
#include "DBScan.h"
#include "AtomicSystem.h"

using namespace std;

int main(int argc, char *argv[])
{
	AtomicSystem MySystem("Dump.cfg");
	Descriptors MyDes(&MySystem,"S","element");
	DBScan MyDB;
	vector<string> Properties;
	Properties.push_back("DBSCAN_EPS AUTO");
	Properties.push_back("DBSCAN_MINPTS AUTO");
	MyDB.ReadProperties(Properties);
	MyDB.setDescriptors(&MyDes);
	MyDB.TrainModel();
	MySystem.setAux_vec(MyDB.getClassificator(),2,"ClusterId");
	MySystem.printSystem_aux("output.cfg","ClusterId S");
	return 0;
}
