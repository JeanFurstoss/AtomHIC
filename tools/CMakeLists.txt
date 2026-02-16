#**********************************************************************************
#*   CMakeLists.txt                                                               *
#**********************************************************************************
#* This file contains the building instructions of an AtomHic executable          *
#**********************************************************************************
#* (C) Jan 2025 - Jean Furstoss                                                   *
#*     Université de Poitiers, Institut PPRIME                                    *
#*     UPR CNRS 3346, 86360 Chasseuneuil-du-Poitou, France                        *
#*     jean.furstoss@univ-poitiers.fr                                             *
#* Last modification: J. Furstoss - 28 Janv 2025                                  *
#**********************************************************************************
#* This program is free software: you can redistribute it and/or modify           *
#* it under the terms of the GNU General Public License as published by           *
#* the Free Software Foundation, either version 3 of the License, or              *
#* (at your option) any later version.                                            *
#*                                                                                *
#* This program is distributed in the hope that it will be useful,                *
#* but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *
#* GNU General Public License for more details.                                   *
#*                                                                                *
#* You should have received a copy of the GNU General Public License              *
#* along with this program.  If not, see <http://www.gnu.org/licenses/>.          *
#**********************************************************************************

# Generic tools
add_subdirectory(AnalyzeCSL)
add_subdirectory(CreateOrientedPlane)
add_subdirectory(CreateGB)
add_subdirectory(CreateFacetGB)
add_subdirectory(ComputeDescriptors)
add_subdirectory(BondOriParam)
add_subdirectory(AtomicStrain)
add_subdirectory(Compute_D2Min)
add_subdirectory(SaveNonCSCrystalBondOriParam)
add_subdirectory(FitAndSaveGMM)
add_subdirectory(FitAndSaveKMeans)
add_subdirectory(GMMClassification)
add_subdirectory(KMeansClassification)
add_subdirectory(DBScanClustering)
add_subdirectory(GenerateYAML4ACE)
# add other subdirectories of new generic tool here

# other tools
add_subdirectory(scripts)
add_subdirectory(compare)
add_subdirectory(MyTools)
