#ifndef MATHTOOLS_H
#define MATHTOOLS_H
#include "AtomHicExport.h"
#include <vector>

class ATOMHIC_EXPORT MathTools {
public:
	MathTools(){};
	// those functions return the max/min of a given array
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
	~MathTools(){};
};

#endif
