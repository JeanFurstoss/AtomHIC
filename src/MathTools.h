#ifndef MATHTOOLS_H
#define MATHTOOLS_H
#include "AtomHicExport.h"
#include "MyStructs.h"
#include <vector>

class ATOMHIC_EXPORT MathTools {
protected:
	double *buffer_vec_1;
	double *buffer_vec_2;
	double *buffer_mat_1;
	double *buffer_mat_2;
public:
	MathTools();
	// those functions return the max/min of a given array
	unsigned int max(const unsigned int arr[], unsigned int size);
	unsigned int min(const unsigned int arr[], unsigned int size);
	int max(const int arr[], unsigned int size);
	int min(const int arr[], unsigned int size);
	double max(const double arr[], unsigned int size);
	double min(const double arr[], unsigned int size);
	// those functions return the indice of the max/min of the vector
	unsigned int max(std::vector<double> arr);
	unsigned int min(std::vector<double> arr);
	double dotProd(const double *vector1, const double *vector2);
	void mixedProd(const double *vector1, const double *vector2, double *prod);
	double gaussian(double x, double mu, double sigma);
	double gaussian_prefac(double x, double mu, double sigma, double prefac);
	void invert3x3(const double *mat, double *inv);
	void gaussian_fit(const std::vector<double> data, double &mu, double &sigma, double &prefac);
	void Vec2rotMat(const double *vec, const double &theta, double *rotMat);
	void MatDotVec(const double *mat, const double *vec, double *prod);
	void VecDotMat(const double *vec, const double *mat, double *prod);
	void MatDotAt(const double *mat, const Atom &At, Atom &At_prod);
	void MatDotMat(const double *mat1, const double *mat2, double *prod);
	void printMat(const double *mat);
	void printVec(const double *vec);
	void sort(const std::vector<double> vec, const unsigned int col, const unsigned int NbCol, std::vector<double> &sorted);
	~MathTools();
};

#endif
