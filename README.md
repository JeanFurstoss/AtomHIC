# AtomHIC: Atomic systems Handling In Computational simulations 

       █████╗ ████████╗ ██████╗ ███╗   ███╗██╗  ██╗██╗ ██████╗
      ██╔══██╗╚══██╔══╝██╔═══██╗████╗ ████║██║  ██║██║██╔════╝
      ███████║   ██║   ██║   ██║██╔████╔██║███████║██║██║
      ██╔══██║   ██║   ██║   ██║██║╚██╔╝██║██╔══██║██║██║
      ██║  ██║   ██║   ╚██████╔╝██║ ╚═╝ ██║██║  ██║██║╚██████╗
      ╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝ ╚═════╝
		     (C) 2025 Jean Furstoss

## Description

AtomHIC is a C++ library which can be used to generate or post-treat atomic systems for atomistic simulations (e.g. molecular dynamics, ab-initio)

Its compilation generates several executables which can directly used. The principal executables allows for instance to:  
	* create bicrystalline orthogonal simulation cell containing a given grain boundary (by giving the five macroscopic degree of freedom)  
	* compute bond orientational parameters even for complex crystals allowing to highlight defective atomic environments  
	* employ machine-learning techniques to classify local atomic environments  

Executables tailored to specific needs can also be created using the different features of the library such as:  
	* readers and printers of atomic systems  
	* crystallographic and bicrystallographic tools adapted to non-cubic and multi-element crystals  
	* managing of crystal and machine-learning database  
	* efficient N-dimensional neighbor research  
	* computation of atomic descriptors  
	* machine-learning models  

This library is intended to be further developed in the future by, among others, adding new descriptors and machine-learning models. Contributors are then welcomed and could contact jean.furstoss@univ-poitiers.fr for helping in any development needs.

## Table of contents

1. [Prerequisite](#prerequisites)
2. [Compilation and execution](#compilation-and-execution)
3. [Main features](#main-features)
4. [Citations](#citations)
5. [Code structure](#code-structure)
6. [Bugs and further developments](#bugs-and-further-developments)
7. [License](#license)

---

## Prerequisites

- Linux distribution (not tested on windows for now)
- [CMake](https://cmake.org/) (version >= 2.8.9)
- C++ compiler, [gcc](https://gcc.gnu.org/) (version >= 8)
- [OpenMP](https://www.openmp.org/) (optional)

## Compilation and execution

Download and compile the AtomHIC repository:

```bash
git clone https://github.com/JeanFurstoss/AtomHIC.git 
cd AtomHIC
mkdir build
cd build
cmake ..
make
```

At this point you can run the test cases to be sure that all works properly by typing 
```bash
ctest
```
in the build directory and verify that 100% of tests pass

All executables are build in the /build/bin/ directory  
To use them just type path_to_bin/NameOfExecutable (a short description of the each executable and the required arguments will then be printed)  
Some part of the code are parallelized using openmp, to benefit for this, use: 
```bash
export OMP_NUM_THREADS="number of thread you want to use"
```
, before using the wanted executable  

## Main features

**1. Creating bicrystalline systems**

The executable CreateGB allows to create an atomic system containing a grain boundary (GB) by giving the five macroscopic degrees of freedom (misorientation axis, angle and GB plane).  
To create the GB, the program uses the Crystal database present in (/data/Crystal/). Thus, if it is not already present in the database, the crystal has to be defined following the example file present in /data/ExampleFiles/Crystal.ath. This file can easily be produced from some crystal databases such as the [AMS](https://rruff.geo.arizona.edu/AMS/amcsd.php) one by downloading the .cif file and converting it to .lmp using [atomsk](https://atomsk.univ-lille.fr/).  

Once the crystal is present in the database, one can create any GB from the executable. For instance, 
```bash
./CreateGB 1 0 0 60.8 0 1 1 Forsterite 1
```

will create the closest rational forsterite (Mg2SiO4) GB to a 60° misoriented around [100] axis and with a (011) GB plane such that:

![figureGB](/doc/GB_And_BondOriParam/figureGB.png)

**2. Compute bond orientational parameter**

Computing order parameter for complex crystal structure (e.g. non-centrosymmetric such as forsterite) is not a trivial task.   
It necessitates firstly to compute some reference bond orientational parameters for the different crystallographic sites of the crystal.  
This can be done using the SaveNonCSCrystalBondOriParam executable which will identify the different sites and store the values needed for the computation of the order parameter in the crystal database.  
Once these values are stored, the order parameter can be easily computed using the BondOriParam executable. More detail on the computation of the order parameter can be found [here](https://www.sciencedirect.com/science/article/pii/S0927025623007620)  
This methodology applied to a bicrystalline forsterite system leads to:
 
![figure_order](/doc/GB_And_BondOriParam/figure_order.png)

**3. Machine-learning techniques for structural analysis**

Identifying specific local atomic environments in atomic systems is a major challenge in material physics. AtomHIC provides tools for performing this task using a supervised classification technique called Steinhardt Gaussian Mixture Analysis (SGMA).  
More details on this technique can be found [here](https://www.sciencedirect.com/science/article/pii/S001046552400403X) and an extensive documentation for applying this technique using AtomHIC can be found in /doc/GaussianMixtureAnalysis.txt  
Briefly, this technique lies with the labeling of a gaussian mixture model fitted on the Steinhardt parameters atomic descriptors such as schematized by:  

![figure_GMM](/doc/FitGMM/FigureFitGMM.png)

## Code structure

The library is structure following:

* /src  
	* % implementation of the different classes of the library  
* /tests  
	* % the different tests of the library (run using ctest)  
* /cmake  
	* % cmake files for the compilation of the library  
* /doc  
	* % some documentation  
* /data  
	* /ExampleFiles  
		* % Contains example files  
	* /Crystals  
		* % The crystal database  
	* /Masses  
		* % Masses database  
	* /FixedParameters  
		* % The parameters used by the different executables (no need to recompile when changing the values)  
	* /MachineLearningModels  
		* % Machine-learning databases  
		* /GaussianMixtureModel  
			* % Existing database for SGMA method  
* /tools  
	* % directories containing sources of improved executables (CreateGB, BondOriParam, ..)  
	* /PythonToolsGMM  
		* % some python tools related to gaussian mixture models  
	* /MyTools  
		* % directory used to develop specific executables  
	* /compare  
		* % tools for ctest  
	* /scripts  
		* % tools for ctest  

The file /doc/CodeStructure.txt gives an overview of the different classes used in the library.  
To develop tailored executables, one can use the /data/ExampleFiles/Tool.cpp to see the main functions that can be used from the library or refer to the comments in the source code.  

## Citations

Please cite following papers:

* [Furstoss, J., Hirel, P., Carrez, P., Gouriet, K., Meko-Fotso, V., & Cordier, P. (2024). Structures and energies of twist grain boundaries in Mg2SiO4 forsterite. Computational Materials Science, 233, 112768.](https://www.sciencedirect.com/science/article/pii/S0927025623007620)
* [Furstoss, J., Salazar, C. R., Carrez, P., Hirel, P., & Lam, J. (2025). All-around local structure classification with supervised learning: The example of crystal phases and dislocations in complex oxides. Computer Physics Communications, 109480.](https://www.sciencedirect.com/science/article/pii/S001046552400403X)


## Bugs and further developments

For any bugs or unexpected behavior of the program, you can either try to fix the issues by yourself using a debugger such as [gdb](https://www.sourceware.org/gdb/) or don't hesitate to contact jfurstoss@orange.fr  

If any feature you need is missing in the code, don't hesitate to contact jfurstoss@orange.fr for helping on the development of the missing feature.  

## License 

(C) Jean Furstoss 2025  
This program is distributed under the GNU/GPL  
(General Public License) version 3 or any later version.  
A copy of this license can be found in the file LICENSE  
that is provided with this program, or on the Web at:  
<http://www.gnu.org/licenses/gpl.html>
