// AtomHic library files
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "MathTools.h"
#include "Descriptors.h"
#include "GaussianMixtureModel.h"

using namespace std;

int main(int argc, char *argv[])
{
	string InputFilename = argv[1];
	Descriptors MyDescriptors(InputFilename);
	GaussianMixtureModel GMM;
	GMM.setDescriptors(&MyDescriptors);
	GMM.TrainModel();
	ifstream file(InputFilename, ios::in);
	vector<vector<double>> data;
	double x1,x2,x3;
	string line;
	if(file){
		while(getline(file,line)){
			istringstream text(line);
			text >> x1 >> x2 >> x3;
			data.push_back(vector<double>());
			data[data.size()-1].push_back(x1);
			data[data.size()-1].push_back(x2);
			data[data.size()-1].push_back(x3);
		}
	}
	double mu_1_1 = 1.1;
	double mu_1_2 = 1.9;
	double mu_1_3 = 3.1;
	double cov_1_11 = 0.1;
	double cov_1_22 = 0.22;
	double cov_1_33 = 0.32;
	double cov_1_12 = 0.023;
	double cov_1_13 = 0.034;
	double cov_1_23 = 0.042;
	vector<vector<double>> cov_1;
	vector<double> mu_1;
	for(unsigned int i=0;i<3;i++) cov_1.push_back(vector<double>());
	cov_1[0].push_back(cov_1_11);
	cov_1[0].push_back(cov_1_12);
	cov_1[0].push_back(cov_1_13);
	cov_1[1].push_back(cov_1_12);
	cov_1[1].push_back(cov_1_22);
	cov_1[1].push_back(cov_1_23);
	cov_1[2].push_back(cov_1_13);
	cov_1[2].push_back(cov_1_23);
	cov_1[2].push_back(cov_1_33);
	mu_1.push_back(mu_1_1);
	mu_1.push_back(mu_1_2);
	mu_1.push_back(mu_1_3);
	MathTools MT;
	vector<vector<double>> inv_cov_1;
	long double det_cov_1;
	MT.invMat_LU(cov_1,inv_cov_1,det_cov_1);

	double mu_2_1 = 0.1;
	double mu_2_2 = 2.9;
	double mu_2_3 = 4.1;
	double cov_2_11 = 0.5;
	double cov_2_22 = 0.12;
	double cov_2_33 = 0.22;
	double cov_2_12 = 0.083;
	double cov_2_13 = 0.014;
	double cov_2_23 = 0.092;
	vector<vector<double>> cov_2;
	vector<double> mu_2;
	for(unsigned int i=0;i<3;i++) cov_2.push_back(vector<double>());
	cov_2[0].push_back(cov_2_11);
	cov_2[0].push_back(cov_2_12);
	cov_2[0].push_back(cov_2_13);
	cov_2[1].push_back(cov_2_12);
	cov_2[1].push_back(cov_2_22);
	cov_2[1].push_back(cov_2_23);
	cov_2[2].push_back(cov_2_13);
	cov_2[2].push_back(cov_2_23);
	cov_2[2].push_back(cov_2_33);
	mu_2.push_back(mu_2_1);
	mu_2.push_back(mu_2_2);
	mu_2.push_back(mu_2_3);
	vector<vector<double>> inv_cov_2;
	long double det_cov_2;
	MT.invMat_LU(cov_2,inv_cov_2,det_cov_2);

	vector<vector<vector<double>>> inv_cov, inv_cov_0, cov;
	vector<vector<double>> mu, mu_0;
	vector<long double> det_cov, det_cov_0;
	vector<double> weight, weight_0;
	
	cov.push_back(cov_2);
	cov.push_back(cov_1);
	
	inv_cov.push_back(inv_cov_2);
	mu.push_back(mu_2);
	det_cov.push_back(det_cov_2);
	weight.push_back(0.6);
	
	inv_cov.push_back(inv_cov_1);
	mu.push_back(mu_1);
	det_cov.push_back(det_cov_1);
	weight.push_back(0.4);
	
	inv_cov_0.push_back(inv_cov_2);
	mu_0.push_back(mu_2);
	det_cov_0.push_back(det_cov_2);
	weight_0.push_back(0.6);
	
	inv_cov_0.push_back(inv_cov_1);
	mu_0.push_back(mu_1);
	det_cov_0.push_back(det_cov_1);
	weight_0.push_back(0.4);
	
	double bic;
	long double temp;
	double L_0 = MT.LogLikelihoodGMM(inv_cov,mu,det_cov,weight,data,bic);
	cout << "Iter " << 0 << ", L = " << L_0 << ", bic = " << bic << endl;
	//for(unsigned int k=0;k<inv_cov.size();k++){
	//	cout << "Cluster " << k << endl;
	//	cout << "weight = " << weight[k] << endl;
	//       	cout << "mu = ";	
	//	for(unsigned int i=0;i<3;i++) cout << mu[k][i] << " ";
	//	cout << endl;
	//	cout << "V_inv = " << endl; 
	//	for(unsigned int i=0;i<3;i++){
	//		for(unsigned int j=0;j<3;j++) cout << inv_cov[k][i][j] << " ";
	//		cout << endl;
	//	}
	//}
	//double L = MT.LogLikelihoodMultidimGaussian(inv_cov_1,mu_1,det_cov_1,data,bic);
	unsigned int nbMaxIterEM = 500;
	double LKH_tol = 1e-4;
	double eps, L;
	for(unsigned int it=0;it<nbMaxIterEM;it++){
		L = MT.ExpectationMaximization_GMM(inv_cov_0,mu_0,det_cov_0,weight_0,inv_cov,mu,det_cov,weight,data,bic);
		eps = L-L_0;
		for(unsigned int k=0;k<inv_cov.size();k++) MT.invMat_LU(inv_cov[k],cov[k],temp);
		//cout << "Iter " << it+1 << ", L = " << L << ", L_diff = " << eps << ", bic = " << bic << endl;
		//for(unsigned int k=0;k<inv_cov.size();k++){
		//	cout << "Cluster " << k << endl;
		//	cout << "weight = " << weight[k] << endl;
		//       	cout << "mu = ";	
		//	for(unsigned int i=0;i<3;i++) cout << mu[k][i] << " ";
		//	cout << endl;
		//	cout << "V_inv = " << endl; 
		//	for(unsigned int i=0;i<3;i++){
		//		for(unsigned int j=0;j<3;j++) cout << cov[k][i][j] << " ";
		//		cout << endl;
		//	}
		//}
		for(unsigned int d=0;d<inv_cov.size();d++){
			weight_0[d] = weight[d];
			//det_cov_0[d] = det_cov[d];
			for(unsigned int x1=0;x1<mu[d].size();x1++){
				mu_0[d][x1] = mu[d][x1];
				for(unsigned int x2=0;x2<mu[d].size();x2++) inv_cov_0[d][x1][x2] = inv_cov[d][x1][x2];
			}
		}
		if( eps < LKH_tol ) break;
		L_0 = L;
	}
		for(unsigned int k=0;k<inv_cov.size();k++){
			cout << "Cluster " << k << endl;
			cout << "weight = " << weight[k] << endl;
		       	cout << "mu = ";	
			for(unsigned int i=0;i<3;i++) cout << mu[k][i] << " ";
			cout << endl;
			cout << "V_inv = " << endl; 
			for(unsigned int i=0;i<3;i++){
				for(unsigned int j=0;j<3;j++) cout << cov[k][i][j] << " ";
				cout << endl;
			}
		}

	cout << bic << endl;

}
