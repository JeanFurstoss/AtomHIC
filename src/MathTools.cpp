#include <iostream>
#include <vector>
#include "MathTools.h"
#include <cmath>

using namespace std;

MathTools::MathTools(){
	this->buffer_vec_1 = new double[3];
	this->buffer_vec_2 = new double[3];
	this->buffer_mat_1 = new double[9];
	this->buffer_mat_2 = new double[9];
}
unsigned int MathTools::max(const unsigned int arr[], unsigned int size){
	unsigned int max = arr[0];
	for(unsigned int i=0;i<size;i++) if( arr[i] > max ) max = arr[i];
	return max;
}

unsigned int MathTools::min(const unsigned int arr[], unsigned int size){
	unsigned int min = arr[0];
	for(unsigned int i=0;i<size;i++) if( arr[i] < min ) min = arr[i];
	return min;
}
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

double MathTools::min_p(const double* arr, unsigned int size){
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

void MathTools::Vec2rotMat(const double *vec, const double &theta, double *rotMat){
	double c = cos(theta);
	double s = sin(theta);
	rotMat[0] = pow(vec[0],2.)*(1.-c)+c;
	rotMat[1] = vec[0]*vec[1]*(1.-c)-(vec[2]*s);
	rotMat[2] = vec[0]*vec[2]*(1.-c)+(vec[1]*s);
	rotMat[3] = vec[0]*vec[1]*(1.-c)+(vec[2]*s);
	rotMat[4] = pow(vec[1],2.)*(1.-c)+c;
	rotMat[5] = vec[1]*vec[2]*(1.-c)-(vec[0]*s);
	rotMat[6] = vec[0]*vec[2]*(1.-c)-(vec[1]*s);
	rotMat[7] = vec[1]*vec[2]*(1.-c)+(vec[0]*s);
	rotMat[8] = pow(vec[2],2.)*(1.-c)+c;
}

void MathTools::MatDotVec(const double *mat, const double *vec, double *prod){
	for(unsigned int i=0;i<3;i++){
		this->buffer_vec_1[i] = vec[i];
		prod[i] = 0.;
	}
	for(unsigned int i=0;i<3;i++){
		for(unsigned int j=0;j<3;j++) prod[i] += mat[i*3+j]*this->buffer_vec_1[j];
	}
}

void MathTools::VecDotMat(const double *vec, const double *mat, double *prod){
	for(unsigned int i=0;i<3;i++){
		this->buffer_vec_1[i] = vec[i];
		prod[i] = 0.;
	}
	for(unsigned int i=0;i<3;i++){
		for(unsigned int j=0;j<3;j++) prod[i] += mat[j*3+i]*this->buffer_vec_1[j];
	}
}

void MathTools::MatDotMat(const double *mat1, const double *mat2, double *prod){
	for(unsigned int i=0;i<9;i++){
		this->buffer_mat_1[i] = mat1[i];
		this->buffer_mat_2[i] = mat2[i];
	}
	for(unsigned int i=0;i<3;i++){
		for(unsigned int j=0;j<3;j++){
			prod[i*3+j] = 0.;
			for(unsigned int i_m=0;i_m<3;i_m++) prod[i*3+j] += this->buffer_mat_1[i_m*3+j]*this->buffer_mat_2[i*3+i_m];
		}
	}
}

void MathTools::printMat(const double *mat){
	for(unsigned int i=0;i<3;i++){
		for(unsigned int j=0;j<3;j++){
			cout << mat[i*3+j] << " ";
		}
		cout << endl;
	}
}
void MathTools::printVec(const double *vec){
	for(unsigned int i=0;i<3;i++) cout << vec[i] << " " ;
	cout << endl;
}

void MathTools::MatDotAt(const double *mat, const Atom &At, Atom &At_prod){
	this->buffer_vec_2[0] = At.pos.x;
	this->buffer_vec_2[1] = At.pos.y;
	this->buffer_vec_2[2] = At.pos.z;
	MatDotVec(mat,buffer_vec_2,buffer_vec_2);
	At_prod.pos.x = this->buffer_vec_2[0];
	At_prod.pos.y = this->buffer_vec_2[1];
	At_prod.pos.z = this->buffer_vec_2[2];
}

void MathTools::sort(const vector<double> vec, const unsigned int col, const unsigned int NbCol, vector<double> &sorted){
	vector<double> vec_j = vec;
	unsigned int init_size = vec.size()/NbCol;
	unsigned int jmin;
	for(unsigned int i=0;i<init_size;i++){
		jmin = 0;
		for(unsigned int j=0;j<vec_j.size()/NbCol;j++) if( vec_j[j*NbCol+col] < vec_j[jmin*NbCol+col] ) jmin = j;
		for(unsigned int n=0;n<NbCol;n++){
			sorted[i*NbCol+n] = vec_j[jmin*NbCol+n];
		}
		for(unsigned int n=0;n<NbCol;n++){
			vec_j.erase(vec_j.begin()+jmin*NbCol);
		}
	}
}

MathTools::~MathTools(){
	delete[] buffer_mat_1;
	delete[] buffer_mat_2;
	delete[] buffer_vec_1;
	delete[] buffer_vec_2;
}
