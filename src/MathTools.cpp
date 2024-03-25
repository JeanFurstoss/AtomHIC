#include <iostream>
#include <vector>
#include "MathTools.h"
#include <cmath>

using namespace std;

MathTools::MathTools(){
	this->buffer_vec_1 = new double[3];
	this->buffer_vec_2 = new double[3];
	this->buffer_vec_3 = new double[3];
	this->buffer_vec_4 = new double[3];
	this->buffer_vec_int = new int[3];
	this->buffer_vec_uint = new unsigned int[3];
	this->buffer_mat_1 = new double[9];
	this->buffer_mat_2 = new double[9];
	this->buffer_mat_uint = new unsigned int[9];
}

double MathTools::min_vec(const vector<double> arr){
	if( arr.size() > 0 ){
		double min = arr[0];
		for(unsigned int i=0;i<arr.size()-1;i++) if( arr[i+1] < min ) min = arr[i+1];
		return min;
	}else{
		return 0.;
	}
}

double MathTools::min_vec_abs(const vector<double> arr){
	if( arr.size() > 0 ){
		double min = arr[0];
		for(unsigned int i=0;i<arr.size()-1;i++) if( fabs(arr[i+1]) < fabs(min) ) min = arr[i+1];
		return min;
	}else{
		return 0.;
	}
}

double MathTools::max_vec(const vector<double> arr){
	if( arr.size() > 0 ){
		double max = arr[0];
		for(unsigned int i=0;i<arr.size()-1;i++) if( arr[i+1] > max ) max = arr[i+1];
		return max;
	}else{
		return 0.;
	}
}

unsigned int MathTools::max(const unsigned int arr[], unsigned int size){
	unsigned int max = arr[0];
	for(unsigned int i=0;i<size;i++) if( arr[i] > max ) max = arr[i];
	return max;
}

unsigned int MathTools::max_ind(const unsigned int arr[], unsigned int size){
	unsigned int max = 0;
	for(unsigned int i=0;i<size;i++) if( arr[i] > arr[max] ) max = i;
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

double MathTools::max_p(const double* arr, unsigned int size){
	double max = arr[0];
	for(unsigned int i=0;i<size;i++) if( arr[i] > max ) max = arr[i];
	return max;
}

unsigned int MathTools::min_p(const unsigned int* arr, unsigned int size){
	unsigned min = arr[0];
	for(unsigned int i=0;i<size;i++) if( arr[i] < min ) min = arr[i];
	return min;
}

unsigned int MathTools::max_p(const unsigned int* arr, unsigned int size){
	unsigned max = arr[0];
	for(unsigned int i=0;i<size;i++) if( arr[i] > max ) max = arr[i];
	return max;
}

unsigned int MathTools::max(vector<double> arr){
	unsigned int imax = 0;
	for(unsigned int i=1;i<arr.size();i++) if( arr[i] > arr[imax] ) imax = i;
	return imax;
}

unsigned int MathTools::max(vector<long double> arr){
	unsigned int imax = 0;
	for(unsigned int i=1;i<arr.size();i++) if( arr[i] > arr[imax] ) imax = i;
	return imax;
}

unsigned int MathTools::min(vector<double> arr){
	unsigned int imin = 0;
	for(unsigned int i=1;i<arr.size();i++) if( arr[i] < arr[imin] ) imin = i;
	return imin;
}

unsigned int MathTools::min_abs(vector<double> arr){
	unsigned int imin = 0;
	for(unsigned int i=1;i<arr.size();i++) if( fabs(arr[i]) < fabs(arr[imin]) ) imin = i;
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

void MathTools::plane_fit(const vector<vector<double>> & data, double &a, double &b, double &c){
	unsigned int nbPts = data.size();
	double *A = new double[9];
	double *B = new double[3];
	double *sol = new double[3];
	double *invA = new double[9];
	for(unsigned int i=0;i<8;i++) A[i] = 0.;
	A[8] = nbPts;
	for(unsigned int i=0;i<3;i++) B[i] = 0.;
	for(unsigned int i=0;i<nbPts;i++){
		A[0] += data[i][0]*data[i][0];
		A[1] += data[i][0]*data[i][1];
		A[2] += data[i][0];
		A[3] += data[i][0]*data[i][1];
		A[4] += data[i][1]*data[i][1];
		A[5] += data[i][1];
		A[6] += data[i][0];
		A[7] += data[i][1];
		B[0] += data[i][0]*data[i][2];
		B[1] += data[i][1]*data[i][2];
		B[2] += data[i][2];
	}
	this->invert3x3(A,invA);
	this->VecDotMat(B,invA,sol);
	a = sol[0];
	b = sol[1];
	c = sol[2];

	delete[] sol;
	delete[] invA;
	delete[] A;
	delete[] B;
}

double MathTools::det(const double *mat){
	double det = 0.;
	det += mat[0]*(mat[4]*mat[8]-(mat[5]*mat[7]));
	det -= mat[1]*(mat[3]*mat[8]-(mat[5]*mat[6]));
	det += mat[2]*(mat[3]*mat[7]-(mat[4]*mat[6]));
	return det;
}
void MathTools::invert3x3(const double *mat, double *inv){
	// compute determinant
	double deter = det(mat);
	if( fabs(deter) < 1e-50 ){
		cout << "Trying to inverse a null determinant matrix !" << endl;
		return;
	}
	for(unsigned int i=0;i<9;i++) buffer_mat_1[i] = mat[i];

	//inv[0] = (mat[4]*mat[8]-(mat[5]*mat[7]))/deter; // to see if it has not change other results
	//inv[1] = -(mat[3]*mat[8]-(mat[5]*mat[6]))/deter;
	//inv[2] = (mat[3]*mat[7]-(mat[4]*mat[6]))/deter;
	//inv[3] = -(mat[1]*mat[8]-(mat[2]*mat[7]))/deter;
	//inv[4] = (mat[0]*mat[8]-(mat[2]*mat[6]))/deter;
	//inv[5] = -(mat[0]*mat[7]-(mat[6]*mat[1]))/deter;
	//inv[6] = (mat[1]*mat[5]-(mat[2]*mat[4]))/deter;
	//inv[7] = -(mat[0]*mat[5]-(mat[2]*mat[3]))/deter;
	//inv[8] = (mat[0]*mat[4]-(mat[1]*mat[3]))/deter;
	inv[0] = (buffer_mat_1[4]*buffer_mat_1[8]-(buffer_mat_1[5]*buffer_mat_1[7]))/deter;
	inv[3] = -(buffer_mat_1[3]*buffer_mat_1[8]-(buffer_mat_1[5]*buffer_mat_1[6]))/deter;
	inv[6] = (buffer_mat_1[3]*buffer_mat_1[7]-(buffer_mat_1[4]*buffer_mat_1[6]))/deter;
	inv[1] = -(buffer_mat_1[1]*buffer_mat_1[8]-(buffer_mat_1[2]*buffer_mat_1[7]))/deter;
	inv[4] = (buffer_mat_1[0]*buffer_mat_1[8]-(buffer_mat_1[2]*buffer_mat_1[6]))/deter;
	inv[7] = -(buffer_mat_1[0]*buffer_mat_1[7]-(buffer_mat_1[6]*buffer_mat_1[1]))/deter;
	inv[2] = (buffer_mat_1[1]*buffer_mat_1[5]-(buffer_mat_1[2]*buffer_mat_1[4]))/deter;
	inv[5] = -(buffer_mat_1[0]*buffer_mat_1[5]-(buffer_mat_1[2]*buffer_mat_1[3]))/deter;
	inv[8] = (buffer_mat_1[0]*buffer_mat_1[4]-(buffer_mat_1[1]*buffer_mat_1[3]))/deter;
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
		for(unsigned int j=0;j<3;j++) prod[i] += mat[j*3+i]*this->buffer_vec_1[j];
	}
}

void MathTools::MatDotRawVec(const double *mat, const double *vec, double *prod){
	for(unsigned int i=0;i<3;i++){
		this->buffer_vec_1[i] = vec[i];
		prod[i] = 0.;
	}
	for(unsigned int i=0;i<3;i++){
		for(unsigned int j=0;j<3;j++) prod[i] += mat[i*3+j]*this->buffer_vec_1[j];
	}
}


void MathTools::MatDotVec_vec(const vector<vector<double>> mat, const vector<double> vec, vector<double> &prod){
	unsigned int dim=vec.size();
	if( mat[0].size() != dim ) cerr << "The dimensions of matrix and vector to multiply are not consistent (i.e. " << mat.size() << "x" << mat[0].size() << " * " << dim << ")" << endl;
	vector<double> temp_vec(dim);
	for(unsigned int i=0;i<dim;i++){
		temp_vec[i] = vec[i];
		prod[i] = 0.;
	}
	for(unsigned int i=0;i<mat.size();i++){
		for(unsigned int j=0;j<dim;j++)	prod[i] += mat[i][j]*temp_vec[j];
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
			for(unsigned int i_m=0;i_m<3;i_m++) prod[i*3+j] += this->buffer_mat_1[i*3+i_m]*this->buffer_mat_2[i_m*3+j];
		}
	}
}

void MathTools::MatDotTMat(const double *mat1, const double *mat2, double *prod){
	for(unsigned int i=0;i<9;i++){
		this->buffer_mat_1[i] = mat1[i];
		this->buffer_mat_2[i] = mat2[i];
	}
	
	for(unsigned int d1=0;d1<3;d1++){
		for(unsigned int d2=d1;d2<3;d2++){
			prod[d1*3+d2] = 0.;
			for(unsigned int d=0;d<3;d++) prod[d1*3+d2] += this->buffer_mat_1[d*3+d1]*this->buffer_mat_2[d*3+d2];
		}
	}
}

void MathTools::MatDotMatVec(const vector<vector<double>> mat1, const vector<vector<double>> mat2, vector<vector<double>> &prod){
	unsigned int dim1_1 = mat1.size();
	unsigned int dim1_2 = mat1[0].size();
	unsigned int dim2_1 = mat2.size();
	unsigned int dim2_2 = mat2[0].size();
	if( dim1_2 != dim2_1 ){
		cerr << "Warning the dimensions of the matrix to multiply are not consistent (matrix 1: " << dim1_1 << "x" << dim1_2 << ", matrix 2: " << dim2_1 << "x" << dim2_2 << ")" << endl;
	}else{
		buffer_vec_vec_1.clear();
		buffer_vec_vec_2.clear();
		for(unsigned int i=0;i<mat1.size();i++){
			buffer_vec_vec_1.push_back(vector<double>());
			for(unsigned int j=0;j<mat1[i].size();j++) buffer_vec_vec_1[i].push_back(mat1[i][j]);
		}
		for(unsigned int i=0;i<mat2.size();i++){
			buffer_vec_vec_2.push_back(vector<double>());
			for(unsigned int j=0;j<mat2[i].size();j++) buffer_vec_vec_2[i].push_back(mat2[i][j]);
		}
		prod.clear();
		for(unsigned int i=0;i<dim1_1;i++){
			prod.push_back(vector<double>());
			for(unsigned int j=0;j<dim2_2;j++){
				prod[i].push_back(0.);
				for(unsigned int j_m=0;j_m<dim1_2;j_m++) prod[i][j] += buffer_vec_vec_1[i][j_m]*buffer_vec_vec_2[j_m][j];
			}
		}
	}

}

void MathTools::printMatVec(const vector<vector<double>> mat){
	for(unsigned int i=0;i<mat.size();i++){
		for(unsigned int j=0;j<mat.size();j++) cout << mat[i][j] << " ";
		cout << endl;
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
	MatDotRawVec(mat,buffer_vec_2,buffer_vec_2);
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

void MathTools::sort_abs(const vector<double> vec, const unsigned int col, const unsigned int NbCol, vector<double> &sorted){
	vector<double> vec_j = vec;
	unsigned int init_size = vec.size()/NbCol;
	unsigned int jmin;
	for(unsigned int i=0;i<init_size;i++){
		jmin = 0;
		for(unsigned int j=0;j<vec_j.size()/NbCol;j++) if( fabs(vec_j[j*NbCol+col]) < fabs(vec_j[jmin*NbCol+col]) ) jmin = j;
		for(unsigned int n=0;n<NbCol;n++){
			sorted[i*NbCol+n] = vec_j[jmin*NbCol+n];
		}
		for(unsigned int n=0;n<NbCol;n++){
			vec_j.erase(vec_j.begin()+jmin*NbCol);
		}
	}
}


void MathTools::reduce_vec(const int *vec1, int *vec2){
	if( vec1[0] == 0 && vec1[1] == 0 && vec1[2] == 0 ) return;
	for(unsigned int i=0;i<3;i++){
		this->buffer_vec_uint[i] = abs(vec1[i]);
		this->buffer_vec_int[i] = vec1[i];
	}
	unsigned int ComDen = gcd_mult(this->buffer_vec_uint,3);
	for(unsigned int i=0;i<3;i++) vec2[i] = this->buffer_vec_int[i] / (int) ComDen;
}

unsigned int MathTools::find_integer_vector(const double *vec, double tol, unsigned int sigma, int *int_vec, bool &IsFound){
	bool Find = false;
	unsigned int commonDenom;
	double Res;
	for(unsigned int i=1;i<sigma+1;i++){
		if( fabs((double) i*vec[0]- (double) round((double) i*vec[0])) < tol && fabs((double) i*vec[1]-(double) round((double) i*vec[1])) < tol && fabs((double) i*vec[2]- (double) round((double) i*vec[2])) < tol ){
			commonDenom = i;
			Find = true;
			break;
		}
	}
	if( Find ){
		for(unsigned int i=0;i<3;i++) int_vec[i] = round(commonDenom*vec[i]);
		reduce_vec(int_vec,int_vec);
		IsFound = true;
		return commonDenom;
	}else{
		//cout << "We don't find rationnal vector of : ";
		//printVec(vec);
		//cout << "within lcd = " << sigma << endl;
		IsFound = false;
		return 0;
	}
}

unsigned int MathTools::gcd(const unsigned int nb1, const unsigned int nb2){
	unsigned int HiNb, LoNb;
	if( nb1 == 0 ) return nb2;
	else if( nb2 == 0 ) return nb1;
	else if( nb1 > nb2 ){
		HiNb = nb1;
		LoNb = nb2;
	}else if( nb2 > nb1 ){
		HiNb = nb2;
		LoNb = nb1;
	}else return nb1;
	unsigned int Rem;
	unsigned int count = 0;
	unsigned int MaxCount = 10000;
	Rem = HiNb%LoNb;
	while( Rem != 0 && count < MaxCount ){
		HiNb = LoNb;
		LoNb = Rem;
		Rem = HiNb%LoNb;
		count++;
	}
	if( Rem == 0 ) return LoNb;
	else{
		cout << "We don't have find greatest common divisor of " << nb1 << " and " << nb2 << endl;
		return 0;
	}
}

unsigned int MathTools::gcd_mult(const unsigned int *arr, const unsigned int size){
	unsigned int buffer1_gcd;
	unsigned int buffer2_gcd;
	buffer1_gcd = gcd(arr[0],arr[1]);
	for(unsigned int i=0;i<size-2;i++){
		buffer2_gcd = gcd(buffer1_gcd,arr[i+2]);
		buffer1_gcd = buffer2_gcd;
	}
	return buffer1_gcd;		
}

// lattice reduction
void MathTools::LLL(const double *Mat, double *Prod){
	for(unsigned int i=0;i<9;i++) buffer_mat_1[i] = Mat[i];
	Gram_Schmidt(buffer_mat_1,buffer_mat_2);
	int k = 1;
	unsigned int km1;
	double delta = 3./4.;
	double proj, norm, sp1, sp2;
	unsigned int count = 0;
	unsigned int MaxCount = 100;
	while( k<=2 && count < MaxCount ){
		for(int j=k-1;j>=0;j--){
			norm = 0.;
			proj = 0.;
			for(unsigned int i=0;i<3;i++){
				proj += buffer_mat_1[i*3+ (unsigned int) k]*buffer_mat_2[i*3+ (unsigned int) j];
				norm += pow(buffer_mat_2[i*3+ (unsigned int) j],2.);
			}
			proj /= norm;
			if( fabs(proj) > 1./2. ){
				for(unsigned int i=0;i<3;i++) buffer_mat_1[i*3+k] -= round(proj)*buffer_mat_1[i*3+ (unsigned int) j];
				Gram_Schmidt(buffer_mat_1,buffer_mat_2);
			}
		}
		proj = 0.;
		norm = 0.;
		sp1 = 0.;
		sp2 = 0.;
		if( k <= 0 ) km1 = k+2;
		else km1 = k-1;
		for(unsigned int i=0;i<3;i++){
			proj += buffer_mat_1[i*3+k]*buffer_mat_2[i*3+km1];
			norm += pow(buffer_mat_2[i*3+km1],2.);
			sp1 += buffer_mat_2[i*3+k]*buffer_mat_2[i*3+k];
			sp2 += buffer_mat_2[i*3+km1]*buffer_mat_2[i*3+km1];
		}
		proj /= norm;
		if( sp1 >= ( (delta-pow(proj,2.))*sp2 ) ) k += 1;
		else{
			// exchange col k and k-1 of buffer_mat_1
			for(unsigned int i=0;i<3;i++){
				buffer_vec_1[i] = buffer_mat_1[i*3+k];
				buffer_mat_1[i*3+k] = buffer_mat_1[i*3+km1];
				buffer_mat_1[i*3+km1] = buffer_vec_1[i];
			}
			Gram_Schmidt(buffer_mat_1,buffer_mat_2);
			if( k > 2 ) k -= 1;
			else k = 1;
		}
		count += 1;
	}
	for(unsigned int i=0;i<9;i++) Prod[i] = buffer_mat_1[i];
}

// orthogonalized basis
void MathTools::Gram_Schmidt(const double *Mat, double *Prod){
	double proj, norm;
	for(unsigned int i=0;i<9;i++) buffer_mat_1[i] = Mat[i];
	for(unsigned int i=0;i<3;i++){
		for(unsigned int j=0;j<3;j++) Prod[j*3+i] = buffer_mat_1[j*3+i];
		if( i != 0 ){
			for(unsigned int j=0;j<3;j++){
				for(unsigned int col=0;col<i;col++){
					proj = 0.;
					norm = 0.;
					for(unsigned int jproj=0;jproj<3;jproj++){
						proj += buffer_mat_1[jproj*3+i]*buffer_mat_1[jproj*3+col];
						norm += pow(buffer_mat_1[jproj*3+col],2.);
					}
					proj /= norm;
					Prod[j*3+i] -= proj*buffer_mat_1[j*3+col];
				}
			}
		}
	}
}
// compute the symmetric over the second diagonal
void MathTools::dia_sym_mtx(const double *Mat, double *Prod){
	for(unsigned int i=0;i<9;i++) buffer_mat_1[i] = Mat[i];
	for(unsigned int i=0;i<3;i++){
		for(unsigned int j=0;j<3;j++) Prod[i*3+j] = buffer_mat_1[(2-j)*3+(2-i)];
	}
}

// cross product of two vectors //TODO verify if the two following function work well
void MathTools::crossProd(const double *vec1, const double *vec2, double *Prod){
	for(unsigned int i=0;i<3;i++){
		buffer_vec_1[i] = vec1[i];
		buffer_vec_2[i] = vec2[i];
	}
	Prod[0] = (buffer_vec_1[1]*buffer_vec_2[2]) - (buffer_vec_1[2]*buffer_vec_2[1]);
	Prod[1] = (buffer_vec_1[2]*buffer_vec_2[0]) - (buffer_vec_1[0]*buffer_vec_2[2]);
	Prod[2] = (buffer_vec_1[0]*buffer_vec_2[1]) - (buffer_vec_1[1]*buffer_vec_2[0]);
}

// get right handed lattice of Mat 
void MathTools::get_right_hand(const double *Mat, double *Prod){
	for(unsigned int i=0;i<9;i++){
		buffer_mat_1[i] = Mat[i];
		Prod[i] = Mat[i];
	}
	for(unsigned int i=0;i<3;i++){
		buffer_vec_3[i] = buffer_mat_1[i*3];
		buffer_vec_4[i] = buffer_mat_1[i*3+1];
	}
	crossProd(buffer_vec_3,buffer_vec_4,buffer_vec_3);
	double sp = 0.;
	for(unsigned int i=0;i<3;i++) sp += buffer_vec_4[i]*buffer_mat_1[i*3+2];
	if( sp < 0 ){
		for(unsigned int i=0;i<3;i++) Prod[i*3+2] = -buffer_mat_1[i*3+2];
	}
}
// compute a transformation matrix permitting to transform an inclined parallelepiped (with xh, yh and zh its cell parameters) into an orthogonal cell with dimension |xh| |yh| and |zh|
void MathTools::computeTiltTrans(const double *xh, const double *yh, const double  *zh, double *TiltTrans){
	// box dimension
	double xbox = sqrt( pow(xh[0],2.) + pow(xh[1],2.) + pow(xh[2],2.));
	double ybox = sqrt( pow(yh[0],2.) + pow(yh[1],2.) + pow(yh[2],2.));
	double zbox = sqrt( pow(zh[0],2.) + pow(zh[1],2.) + pow(zh[2],2.));
	// compute corrections to apply to the cell vectors and atomic positions to have an orthogonal cell
	buffer_mat_1[0] = 1.;
	buffer_mat_1[1] = -yh[0]/yh[1]+((zh[0]-(yh[0]*zh[1]/yh[1]))*yh[2]/((zh[2]-(yh[2]*zh[1]/yh[1]))*yh[1]));
	buffer_mat_1[2] = -(zh[0]-(yh[0]*zh[1]/yh[1]))/(zh[2]-(yh[2]*zh[1]/yh[1]));
	buffer_mat_1[3] = 0.;
	buffer_mat_1[4] = ybox/yh[1]+(zh[1]*(ybox/yh[1])*yh[2]/((zh[2]-(yh[2]*zh[1]/yh[1]))*yh[1]));
	buffer_mat_1[5] = -zh[1]*(ybox/yh[1])/(zh[2]-(yh[2]*zh[1]/yh[1]));
	buffer_mat_1[6] = 0.;
	buffer_mat_1[7] = -yh[2]*zbox/(yh[1]*(zh[2]-(yh[2]*zh[1]/yh[1])));
	buffer_mat_1[8] = zbox/(zh[2]-(yh[2]*zh[1]/yh[1]));
	
	TiltTrans[0] = xbox/(xh[0]+buffer_mat_1[1]*xh[1]+buffer_mat_1[2]*xh[2]);
	TiltTrans[1] = buffer_mat_1[1]*xbox/(xh[0]+buffer_mat_1[1]*xh[1]+buffer_mat_1[2]*xh[2]);
	TiltTrans[2] = buffer_mat_1[2]*xbox/(xh[0]+buffer_mat_1[1]*xh[1]+buffer_mat_1[2]*xh[2]);
	TiltTrans[3] = -(xh[1]*buffer_mat_1[4]+buffer_mat_1[5]*xh[2])/(xh[0]+buffer_mat_1[1]*xh[1]+buffer_mat_1[2]*xh[2]);
	TiltTrans[4] = buffer_mat_1[4]-(xh[1]*buffer_mat_1[4]+buffer_mat_1[5]*xh[2])*buffer_mat_1[1]/(xh[0]+buffer_mat_1[1]*xh[1]+buffer_mat_1[2]*xh[2]);
	TiltTrans[5] = buffer_mat_1[5]-(xh[1]*buffer_mat_1[4]+buffer_mat_1[5]*xh[2])*buffer_mat_1[2]/(xh[0]+buffer_mat_1[1]*xh[1]+buffer_mat_1[2]*xh[2]);
	TiltTrans[6] = -(xh[2]*buffer_mat_1[8]+buffer_mat_1[7]*xh[1])/(xh[0]+buffer_mat_1[1]*xh[1]+buffer_mat_1[2]*xh[2]);
	TiltTrans[7] = buffer_mat_1[7]-(xh[2]*buffer_mat_1[8]+buffer_mat_1[7]*xh[1])*buffer_mat_1[1]/(xh[0]+buffer_mat_1[1]*xh[1]+buffer_mat_1[2]*xh[2]);
	TiltTrans[8] = buffer_mat_1[8]-(xh[2]*buffer_mat_1[8]+buffer_mat_1[7]*xh[1])*buffer_mat_1[2]/(xh[0]+buffer_mat_1[1]*xh[1]+buffer_mat_1[2]*xh[2]);
}

void MathTools::MultidimGaussian(const vector<vector<double>> data, vector<double> &mu, vector<vector<double>> &C){
	for(unsigned int i=0;i<C.size();i++) C[i].clear();
	C.clear();
	mu.clear();
	unsigned int nbdat = data.size();
	unsigned int dim = data[0].size();
	for(unsigned int i=0;i<dim;i++){
		mu.push_back(0.);
		C.push_back(vector<double>());
		for(unsigned int j=0;j<dim;j++)	C[C.size()-1].push_back(0.);
	}
	// compute esperance
	for(unsigned int i=0;i<dim;i++){
		for(unsigned int j=0;j<nbdat;j++) mu[i] += data[j][i];
		mu[i] /= nbdat;
	}
	// compute covariance matrix
	for(unsigned int i=0;i<dim;i++){
		for(unsigned int j=0;j<dim;j++){
			for(unsigned int k=0;k<nbdat;k++) C[i][j] += (data[k][i]-mu[i])*(data[k][j]-mu[j]);
			C[i][j] /= nbdat;
		}
	}
}

long double MathTools::Prob_MultidimGaussian(const vector<vector<double>> C_inv, vector<double> mu, const long double det_C, const vector<double> X){
	unsigned int dim=mu.size();
	vector<double> prod(dim);
	vector<double> XMinusMu(dim);
	double sp = 0.;

	for(unsigned int i=0;i<dim;i++) XMinusMu[i] = X[i]-mu[i];
	MatDotVec_vec(C_inv,XMinusMu,prod);
	for(unsigned int i=0;i<dim;i++) sp += XMinusMu[i]*prod[i];
	return ( 1./ ( pow(2.*M_PI, dim/2.) * sqrt(det_C) ) ) * exp( -.5*sp );
}

void MathTools::invMat_LU(const vector<vector<double>> mat, vector<vector<double>> &inv, long double &det){
	unsigned int dim = mat.size();
	if( dim != mat[0].size() ){
		cerr << "We cannot invert the matrix because it is not square matrix" << endl; 
	}
	// Initialisation of L and Lt
	vector<vector<double>> L; //Triangular lower matrix
	vector<vector<double>> U; //Triangular upper matrix => L*U=mat
	vector<vector<double>> L_inv; //Inverse of Triangular lower matrix
	vector<vector<double>> U_inv; //Inverse of Triangular upper matrix => L*U=mat
	vector<vector<double>> temp; // buffer
	vector<double> vec;
	vector<double> v(mat.size(), 0);
	
	for(unsigned int i(0);i<mat.size();i++) {
		L.push_back(v);
		U.push_back(v);
		L_inv.push_back(v);
		U_inv.push_back(v);
		temp.push_back(v);
		for(unsigned int j=0;j<mat[i].size();j++) temp[i][j] = mat[i][j];
	}
	double sum;
	for(unsigned int j=0;j<dim;j++){
		L[j][j] = 1.;
		for(unsigned int i=0;i<j+1;i++){
			sum = 0.;
			for(unsigned int k=0;k<i;k++) sum += U[k][j]*L[i][k];
			U[i][j] = mat[i][j] - sum;
		}
		for(unsigned int i=j;i<dim;i++){
			sum = 0.;
			for(unsigned int k=0;k<j;k++) sum += U[k][j]*L[i][k];
			L[i][j] = ( mat[i][j] - sum ) / U[j][j];
		}
	}
	// Compute the determinant (product of diagonal element of U)
	det = U[0][0];
	for(unsigned int i=1;i<dim;i++) det *= U[i][i];
	//cout << "deter " << det << endl;
	//printMatVec(U);
	// Compute L_inv
	for(unsigned int i=0;i<dim;i++){
		L_inv[i][i] = 1./L[i][i];
		for(unsigned int j=0;j<i;j++){
			for(unsigned int k=j;k<=i-1;k++) L_inv[i][j] -= L[i][k]*L_inv[k][j];
			L_inv[i][j] /= L[i][i];
		}
	}
	// Compute U_inv
	unsigned int i_uint;
	for(int i=dim-1;i>=0;i--){
		i_uint = (unsigned int) i;
		U_inv[i_uint][i_uint] = 1./U[i_uint][i_uint];
		for(unsigned int j=dim-1;j>i_uint;j--){
			for(unsigned int k=i_uint+1;k<=j;k++) U_inv[i_uint][j] -= U[i_uint][k]*U_inv[k][j];
			U_inv[i_uint][j] /= U[i_uint][i_uint];
		}
	}
	// mat-1 = L-1 * U-1
	MatDotMatVec(U_inv,L_inv,inv);
}

MathTools::~MathTools(){
	delete[] buffer_mat_1;
	delete[] buffer_mat_2;
	delete[] buffer_mat_uint;
	delete[] buffer_vec_1;
	delete[] buffer_vec_2;
	delete[] buffer_vec_3;
	delete[] buffer_vec_4;
	delete[] buffer_vec_uint;
	delete[] buffer_vec_int;
}
