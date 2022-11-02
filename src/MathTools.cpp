#include <iostream>
#include <vector>
#include "MathTools.h"
#include <cmath>

using namespace std;

int MathTools::max(const int arr[], unsigned int size){
	int max = arr[0];
	for(unsigned int i=0;i<size;i++) if( arr[i] > max ) max = arr[i];
	return max;
}

int MathTools::min(const int arr[], unsigned int size){
	int min = arr[0];
	for(unsigned int i=0;i<size;i++) if( arr[i] < min ) min = arr[i];
	return min;
}

double MathTools::max(const double arr[], unsigned int size){
	double max = arr[0];
	for(unsigned int i=0;i<size;i++) if( arr[i] > max ) max = arr[i];
	return max;
}

double MathTools::min(const double arr[], unsigned int size){
	double min = arr[0];
	for(unsigned int i=0;i<size;i++) if( arr[i] < min ) min = arr[i];
	return min;
}

unsigned int MathTools::max(vector<double> arr){
	unsigned int imax = 0;
	for(unsigned int i=1;i<arr.size();i++) if( arr[i] > arr[imax] ) imax = i;
	return imax;
}

unsigned int MathTools::min(vector<double> arr){
	unsigned int imin = 0;
	for(unsigned int i=1;i<arr.size();i++) if( arr[i] < arr[imin] ) imin = i;
	return imin;
}

double MathTools::dotProd(const double *vector1, const double *vector2){
	double dotProd = 0;
	for(unsigned int i=0;i<3;i++) dotProd += vector1[i]*vector2[i];
	return dotProd;
}

void MathTools::mixedProd(const double *vector1, const double *vector2, double *prod){
	prod[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
	prod[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
	prod[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
}

double MathTools::gaussian(double x, double mu, double sigma){
        return (1./(sigma*sqrt(M_PI*2.)))*exp(-pow(x-mu, 2.)/(2.*pow(sigma,2.)));
}

double MathTools::gaussian_prefac(double x, double mu, double sigma, double prefac){
        return (prefac/(sigma*sqrt(M_PI*2.)))*exp(-pow(x-mu, 2.)/(2.*pow(sigma,2.)));
}

void MathTools::gaussian_fit(const vector<double> data, double &mu, double &sigma, double &prefac){
	unsigned int nbPts = data.size()/2;
	unsigned int maxIt = 100;
	double eps = 1e-6; // for numerical derivative
	double *jacobian = new double[nbPts*3];
	double *JTJ = new double[9];
	double *invJTJ = new double[9];
	double *residual = new double[nbPts];
	double *JTR = new double[3];
	double delta_mu, delta_sigma, delta_prefac, det;
	double tol = 1e-5;
	double NormRes = 0.;
	unsigned int iter = 1;
	// compute initial residual
	for(unsigned int i=0;i<nbPts;i++){
		residual[i] = data[i*2] - this->gaussian_prefac(data[i*2+1], mu, sigma, prefac);
		NormRes += pow(residual[i],2.);
	}
	NormRes = sqrt(NormRes);
	for(unsigned int n=0;n<maxIt;n++){
		for(unsigned int i=0;i<nbPts;i++){
			// compute Jacobian
			jacobian[i*3] = (this->gaussian_prefac(data[i*2+1], mu+eps, sigma, prefac)-this->gaussian_prefac(data[i*2+1], mu-eps, sigma, prefac))/(2.*eps);
			jacobian[i*3+1] = (this->gaussian_prefac(data[i*2+1], mu, sigma+eps, prefac)-this->gaussian_prefac(data[i*2+1], mu, sigma-eps, prefac))/(2.*eps);
			jacobian[i*3+2] = (this->gaussian_prefac(data[i*2+1], mu, sigma, prefac+eps)-this->gaussian_prefac(data[i*2+1], mu, sigma, prefac-eps))/(2.*eps);
		}
		// Compute jacobian*T(jacobian)
		for(unsigned int i=0;i<3;i++){
			for(unsigned int j=0;j<3;j++){
				JTJ[i*3+j] = 0.;
				for(unsigned int k=0;k<nbPts;k++) JTJ[i*3+j] += jacobian[k*3+j]*jacobian[k*3+i];
			}
		}
		// Invert JTJ
		this->invert3x3(JTJ,invJTJ);
		// Compute Jacobian times transposed residual
		for(unsigned int i=0;i<3;i++){
			JTR[i] = 0.;	
			for(unsigned int k=0;k<nbPts;k++) JTR[i] += jacobian[k*3+i]*residual[k];
		}
		// compute deltas
		delta_mu = (invJTJ[0]*JTR[0]) + (invJTJ[3]*JTR[1]) + (invJTJ[6]*JTR[2]);
		delta_sigma = (invJTJ[1]*JTR[0]) + (invJTJ[4]*JTR[1]) + (invJTJ[7]*JTR[2]);
		delta_prefac = (invJTJ[2]*JTR[0]) + (invJTJ[5]*JTR[1]) + (invJTJ[8]*JTR[2]);
		mu += delta_mu;
		sigma += delta_sigma;
		prefac += delta_prefac;
		// recompute residual
		NormRes = 0.;
		for(unsigned int i=0;i<nbPts;i++){
			residual[i] = data[i*2] - this->gaussian_prefac(data[i*2+1], mu, sigma, prefac);
			NormRes += pow(residual[i],2.);
		}
		NormRes = sqrt(NormRes);
		if( sqrt(pow(delta_mu,2.)+pow(delta_sigma,2.)+pow(delta_prefac,2.)) < tol || NormRes < tol ) break;
		iter++;
	}
	if( iter > maxIt-1 ) cout << "Warning the optimal parameters of gaussian distribution have not been found" << endl;
	// desallocate memory
	delete[] jacobian;
	delete[] JTJ;
	delete[] invJTJ;
	delete[] residual;
	delete[] JTR;
}

void MathTools::invert3x3(const double *mat, double *inv){
	// compute determinant
	double det = 0.;
	det += mat[0]*(mat[4]*mat[8]-(mat[5]*mat[7]));
	det -= mat[1]*(mat[3]*mat[8]-(mat[5]*mat[6]));
	det += mat[2]*(mat[3]*mat[7]-(mat[4]*mat[6]));
	if( fabs(det) < 1e-50 ){
		cout << "Trying to inverse a null determinant matrix !" << endl;
		return;
	}
	inv[0] = (mat[4]*mat[8]-(mat[5]*mat[7]))/det;
	inv[1] = -(mat[3]*mat[8]-(mat[5]*mat[6]))/det;
	inv[2] = (mat[3]*mat[7]-(mat[4]*mat[6]))/det;
	inv[3] = -(mat[1]*mat[8]-(mat[2]*mat[7]))/det;
	inv[4] = (mat[0]*mat[8]-(mat[2]*mat[6]))/det;
	inv[5] = -(mat[0]*mat[7]-(mat[6]*mat[1]))/det;
	inv[6] = (mat[1]*mat[5]-(mat[2]*mat[4]))/det;
	inv[7] = -(mat[0]*mat[5]-(mat[2]*mat[3]))/det;
	inv[8] = (mat[0]*mat[4]-(mat[1]*mat[3]))/det;
}