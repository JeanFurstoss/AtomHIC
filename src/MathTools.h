#ifndef MATHTOOLS_H
#define MATHTOOLS_H
#include "AtomHicExport.h"
#include "MyStructs.h"
#include <vector>

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
	unsigned int min(std::vector<double> arr);
	unsigned int min_abs(std::vector<double> arr);
	// those functions return the max/min of a double vector
	double max_vec(const std::vector<double> arr);
	double min_vec(const std::vector<double> arr);
	double min_vec_abs(const std::vector<double> arr);
	double dotProd(const double *vector1, const double *vector2);
	void mixedProd(const double *vector1, const double *vector2, double *prod);
	void reduce_vec(const int *vec1, int *vec2);
	double gaussian(double x, double mu, double sigma);
	double gaussian_prefac(double x, double mu, double sigma, double prefac);
	double det(const double *mat);
	void invert3x3(const double *mat, double *inv);
	void gaussian_fit(const std::vector<double> data, double &mu, double &sigma, double &prefac);
	void Vec2rotMat(const double *vec, const double &theta, double *rotMat);
	void MatDotVec(const double *mat, const double *vec, double *prod);
	void MatDotVec_vec(const std::vector<std::vector<double>> mat, const std::vector<double> vec, std::vector<double> &prod);
	void VecDotMat(const double *vec, const double *mat, double *prod);
	void MatDotAt(const double *mat, const Atom &At, Atom &At_prod);
	void MatDotMat(const double *mat1, const double *mat2, double *prod);
	void MatDotMatVec(const std::vector<std::vector<double>> mat1, const std::vector<std::vector<double>> mat2, std::vector<std::vector<double>> &prod);
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
	double Prob_MultidimGaussian(const std::vector<std::vector<double>> C_inv, std::vector<double> mu, const double det_C, const std::vector<double> X); // compute the probability associated to a multidim gaussian with a given esperance, inverse covariant matrix and determinant of covariant matrix 
	void invMat_LU(const std::vector<std::vector<double>> mat, std::vector<std::vector<double>> &inv, double &det); // invert a square matrix using the LU method and compute the determinant
	~MathTools();
};

#endif
