//**********************************************************************************
//*   MathTools.h		                                                   *
//**********************************************************************************
//* This file contains the declaration of the MathTools class which contains       *
//* numerous mathematical tools (e.g. linear algebra, array sorting, min max..)    *
//* Almost every class of AtomHic has Its own MathTools for performing simple tasks*
//* Other features can easily be implemented here.				   *
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

#ifndef MATHTOOLS_H
#define MATHTOOLS_H
#include "AtomHicExport.h"
#include "MyStructs.h"
#include <vector>
#include <complex>

class ATOMHIC_EXPORT MathTools {
protected:
	double *buffer_vec_1;
	double *buffer_vec_2;
	double *buffer_vec_3;
	double *buffer_vec_4;
	unsigned int *buffer_vec_uint;
	int *buffer_vec_int;
	double *buffer_mat_1;
	double *buffer_mat_2;
	unsigned int *buffer_mat_uint;
	std::vector<std::vector<double>> buffer_vec_vec_1, buffer_vec_vec_2;
public:
	MathTools();
	// those functions return the max/min of a given array
	unsigned int max(const unsigned int arr[], unsigned int size);
	unsigned int min(const unsigned int arr[], unsigned int size);
	int max(const int arr[], unsigned int size);
	int min(const int arr[], unsigned int size);
	double max(const double arr[], unsigned int size);
	double min(const double arr[], unsigned int size);
	double min_p(const double* arr, unsigned int size);
	double max_p(const double* arr, unsigned int size);
	unsigned int min_p(const unsigned int* arr, unsigned int size);
	unsigned int max_p(const unsigned int* arr, unsigned int size);
	// those functions return the indice of the max/min of the vector
	unsigned int max(std::vector<double> arr);
	unsigned int max(std::vector<long double> arr);
	unsigned int min(std::vector<double> arr);
	unsigned int min_abs(std::vector<double> arr);
	unsigned int max_ind(const unsigned int arr[], unsigned int size);
	unsigned int min_p_ind(const double* arr, unsigned int &size);
	unsigned int max_p_ind(const double* arr, unsigned int &size);
	unsigned int min_p_ind(const double* arr, unsigned int size, unsigned int dim, unsigned int col);
	unsigned int max_p_ind(const double* arr, unsigned int size, unsigned int dim, unsigned int col);
	unsigned int max_p_ind(const long double* arr, unsigned int &size);
	// those functions return the max/min of a double vector
	double max_vec(const std::vector<double> arr);
	unsigned int max_vec(const std::vector<unsigned int> arr);
	double min_vec(const std::vector<double> arr);
	double min_vec_abs(const std::vector<double> arr);
	double dotProd(const double *vector1, const double *vector2);
	void mixedProd(const double *vector1, const double *vector2, double *prod);
	void reduce_vec(const int *vec1, int *vec2);
	double gaussian(double x, double mu, double sigma);
	double gaussian_prefac(double x, double mu, double sigma, double prefac);
	double det(const double *mat);
	double det(const std::vector<std::vector<double>> &mat);
	void invert3x3(const double *mat, double *inv);
	double invert3x3(const std::vector<std::vector<double>> &mat, std::vector<std::vector<double>> &inv);
	void plane_fit(const std::vector<std::vector<double>> & data, double &a, double &b, double &c);
	void gaussian_fit(const std::vector<double> data, double &mu, double &sigma, double &prefac);
	void Vec2rotMat(const double *vec, const double &theta, double *rotMat);
	void MatDotVec(const double *mat, const double *vec, double *prod);
	void MatDotRawVec(const double *mat, const double *vec, double *prod);
	void MatDotRawVec_vec(const std::vector<std::vector<double>> mat, const std::vector<double> vec, std::vector<double> &prod);
	void OuterVecProduct(const std::vector<double> vec1, const std::vector<double> vec2, std::vector<std::vector<double>> &prod);
	void MatDotVec_vec(const std::vector<std::vector<double>> mat, const std::vector<double> vec, std::vector<double> &prod);
	void VecDotMat(const double *vec, const double *mat, double *prod);
	void MatDotAt(const double *mat, const Atom &At, Atom &At_prod);
	void MatDotMat(const double *mat1, const double *mat2, double *prod);
	void MatDotTMat(const double *mat1, const double *mat2, double *prod);
	void MatDotMatVec(const std::vector<std::vector<double>> mat1, const std::vector<std::vector<double>> mat2, std::vector<std::vector<double>> &prod);
	void MatDotMatVec(const long double *mat1, const long double *mat2, long double *prod, unsigned int &dim1_1, unsigned int &dim1_2, unsigned int &dim2_1, unsigned int &dim2_2);
	void printMat(const double *mat);
	void printMatVec(const std::vector<std::vector<double>> mat);
	void printVec(const double *vec);
	void crossProd(const double *vec1, const double *vec2, double *Prod);
	void dia_sym_mtx(const double *Mat, double *Prod); // compute the symmetric over the second diagonal
	void get_right_hand(const double *Mat, double *Prod); // get right handed lattice of Mat 
	unsigned int find_integer_vector(const double *vec, double tol, unsigned int sigma, int *int_vec, bool &IsFound);
	unsigned int gcd(const unsigned int nb1, const unsigned int nb2); // greatest common dividor
	unsigned int gcd_mult(const unsigned int *arr, const unsigned int size); // greatest common dividor
	void LLL(const double *Mat, double *Prod); //lattice reduction
	void Gram_Schmidt(const double *Mat, double *Prod); // orthogonalized basis
	void sort(const std::vector<double> vec, const unsigned int col, const unsigned int NbCol, std::vector<double> &sorted); // sort a vector with respect to a given column (col), giving the number of column in the vector
	void sort_abs(const std::vector<double> vec, const unsigned int col, const unsigned int NbCol, std::vector<double> &sorted); // sort a vector with respect to a given column (col), giving the number of column in the vector
	void computeTiltTrans(const double *xh, const double *yh, const double  *zh, double *TiltTrans); // compute a transformation matrix permitting to transform an inclined parallelepiped (with xh, yh and zh its cell parameters) into an orthogonal cell with dimension |xh| |yh| and |zh|
	void MultidimGaussian(const std::vector<std::vector<double>> data, std::vector<double> &mu, std::vector<std::vector<double>> &C); // compute the esperance and covariance of multidimensionnal data assuming normal law
	long double Prob_MultidimGaussian(const std::vector<std::vector<double>> C_inv, std::vector<double> mu, const long double det_C, const std::vector<double> X); // compute the probability associated to a multidim gaussian with a given esperance, inverse covariant matrix and determinant of covariant matrix 
	long double LogLikelihoodMultidimGaussian(const std::vector<std::vector<double>> C_inv, std::vector<double> mu, const long double det_C, const std::vector<std::vector<double>> data, double &BIC); // compute the log likelihood of a Gaussian distribution regarding a dataset and the Bayes Informed Criteria (BIC) 
	long double LogLikelihoodGMM(const std::vector<std::vector<std::vector<double>>> C_inv, std::vector<std::vector<double>> mu, const std::vector<long double> det_C, const std::vector<double> weight, const std::vector<std::vector<double>> data, double &BIC); // compute the log likelihood of a Gaussian Mixture Model regarding a dataset and the Bayes Informed Criteria (BIC) 
	double ExpectationMaximization_GMM(const std::vector<std::vector<std::vector<double>>> C_inv_0, std::vector<std::vector<double>> mu_0, const std::vector<long double> det_C_0, const std::vector<double> weight_0, std::vector<std::vector<std::vector<double>>> &C_inv, std::vector<std::vector<double>> &mu, std::vector<long double> &det_C, std::vector<double> &weight, const std::vector<std::vector<double>> data, double &BIC); // One iteration of the EM algorithm, return the difference of likelihood with updated parameters of GMM distrib TODO no one iter but all optimization to declare array only once
	void invMat_LU(const std::vector<std::vector<double>> mat, std::vector<std::vector<double>> &inv, long double &det); // invert a square matrix using the LU method and compute the determinant
	void invMat_LU(const long double *mat, long double *inv, long double &det, unsigned int &dim); // invert a square matrix using the LU method and compute the determinant
	void invMat_LU(long double *mat, long double *inv, unsigned int dim, unsigned int index, unsigned int nbFilter, unsigned int filter_value, long double &det); // invert a square matrix using the LU method and compute the determinant
	std::complex<double> spherical_harmonics(const unsigned int& l, int& m, double& theta, double& phi);
	void EigenDecomposition(double *Matrix, unsigned int dim, double *EigenValues, double *EigenVectors); // find eigenvalues and eigenvectors of Matrix and return them sorted from the highest to the lowest eigenvalue (power iterative method)
	void EigenDecomposition(std::vector<std::vector<double>> &Matrix, double *EigenValues, double *EigenVectors); // find eigenvalues and eigenvectors of Matrix and return them sorted from the highest to the lowest eigenvalue (power iterative method)
	void GenerateCombinations(std::vector<int>& combination, int dim, int currentIndex, const int minValue, const int maxValue, std::vector<std::vector<int>>& results);
	std::vector<std::vector<int>> GenerateNDCombinations(int dim, int minValue, int maxValue);
	~MathTools();
};

#endif
