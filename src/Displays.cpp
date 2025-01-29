//**********************************************************************************
//*   Displays.cpp                                                                   *
//**********************************************************************************
//* This file contains the implementation of the Displays class                      *
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


#include "Displays.h"
#include "AtomHicConfig.h"
#include <iostream>
#include <stdlib.h>
#include <omp.h>

using namespace std; 
using namespace std::chrono;

Displays::Displays(){
	time_beg = high_resolution_clock::now();	
}

void Displays::Logo(){

    cout << R"(
       █████╗ ████████╗ ██████╗ ███╗   ███╗██╗  ██╗██╗ ██████╗
      ██╔══██╗╚══██╔══╝██╔═══██╗████╗ ████║██║  ██║██║██╔════╝
      ███████║   ██║   ██║   ██║██╔████╔██║███████║██║██║
      ██╔══██║   ██║   ██║   ██║██║╚██╔╝██║██╔══██║██║██║
      ██║  ██║   ██║   ╚██████╔╝██║ ╚═╝ ██║██║  ██║██║╚██████╗
      ╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝ ╚═════╝
    )";
    cout << "                  (C) 2025 Jean Furstoss " << endl;

    cout << "              =====================================" << endl;
    cout << endl;
    cout << "Calculation running using " << omp_get_max_threads() << " threads" << endl;
    cout << endl;
}

void Displays::ExecutionTime(){
	std::chrono::high_resolution_clock::time_point end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - time_beg);
	cout << "Execution time : " << duration.count()/1.e3 << " s" << endl;
}
