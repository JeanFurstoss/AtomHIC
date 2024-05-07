// AtomHic library files
#include <filesystem>
#include <AtomicSystem.h>
#include <Bicrystal.h>
#include <Crystal.h>
#include <ComputeAuxiliary.h>
#include <MathTools.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include "MyStructs.h"

#include <GaussianUtils/TrainSet.h>
#include <GaussianMixtureModel/ExpectationMaximization.h>
#include <GaussianMixtureModel/GaussianMixtureModel.h>
#include <GaussianMixtureModel/GaussianMixtureModelFactory.h>
#include <Eigen/Eigenvalues>


using namespace std;

bool read_GTFile(string GTFilename, vector<unsigned int> &IdForGrains){
	ifstream file_i(GTFilename, ios::in);
	string line;
	if( file_i ){
		getline(file_i,line);
		while(file_i){
			getline(file_i,line);
			IdForGrains.push_back(0);
		}
	}else{
		cout << "cannot open grain tag file" << endl;
		return false;
	}
	ifstream file(GTFilename, ios::in);
	getline(file,line);
	unsigned int IdIon, vecRef;
	for(unsigned int i=0;i<IdForGrains.size();i++){
		getline(file,line);
		istringstream text(line);
		text >> IdIon >> vecRef;
		unsigned int idGrain = floor(vecRef/10.);
		unsigned int cellParam = vecRef-idGrain*10;
		IdForGrains[(idGrain-1)*4+cellParam-1] = IdIon;
	}
	return true;
			
}


bool ComputeTransVec(double *coords, unsigned int nbAt, double *H1, double *H2, double *H3, double *G1, double *G2, double *G3, double &transvec_x, double &transvec_y, double &transvec_z){
	double tolBound = 2.; // in A, for an ions to be consider in contact with box boundaries
	double tolRedX = 2./H1[0]; 
	double tolRedY = 2./H2[1]; 
	double tolRedZ = 2./H3[2]; 
	bool XCentered = true;
	bool YCentered = true;
	bool ZCentered = true;
	for(unsigned int i=0;i<nbAt;i++){
		// Compute reduced coordinates
		double xpos = coords[i*3]*G1[0]+coords[i*3+1]*G2[0]+coords[i*3+2]*G3[0];
		double ypos = coords[i*3]*G1[1]+coords[i*3+1]*G2[1]+coords[i*3+2]*G3[1];
		double zpos = coords[i*3]*G1[2]+coords[i*3+1]*G2[2]+coords[i*3+2]*G3[2];
		coords[i*3] = xpos;
		coords[i*3+1] = ypos;
		coords[i*3+2] = zpos;
		if( coords[i*3] < tolRedX || coords[i*3] > 1.-tolRedX ) XCentered = false;
		if( coords[i*3+1] < tolRedY || coords[i*3+1] > 1.-tolRedY ) YCentered = false;
		if( coords[i*3+2] < tolRedZ || coords[i*3+2] > 1.-tolRedZ ) ZCentered = false;
	}
	transvec_x = 0.;
	transvec_y = 0.;
	transvec_z = 0.;
	unsigned int max_count = 10;
	unsigned int count = 0;
	double *TempVec = new double[3];
	while( ( !XCentered || !YCentered || !ZCentered ) && ( count < max_count ) ){
		if( !XCentered ){
			XCentered = true;
			YCentered = true;
			ZCentered = true;
			transvec_x += 1./( (double) max_count );
			for(unsigned int i=0;i<nbAt;i++){
				// Change reduced coordinates
				coords[i*3] += 1./( (double) max_count );
				// Wrapped reduced coordinates
				if( coords[i*3] >= 1. || coords[i*3] < 0. ) coords[i*3] = coords[i*3]-floor(coords[i*3]);
				// Wrapped normal coordinates
				for(unsigned int k=0;k<3;k++) TempVec[k] = coords[i*3]*H1[k] + coords[i*3+1]*H2[k] + coords[i*3+2]*H3[k];
				// New reduced coordinates
				for(unsigned int k=0;k<3;k++) coords[i*3+k] = TempVec[0]*G1[k] + TempVec[1]*G2[k] + TempVec[2]*G3[k];
				if( coords[i*3] < tolRedX || coords[i*3] > 1.-tolRedX ) XCentered = false;
				if( coords[i*3+1] < tolRedY || coords[i*3+1] > 1.-tolRedY ) YCentered = false;
				if( coords[i*3+2] < tolRedZ || coords[i*3+2] > 1.-tolRedZ ) ZCentered = false;
			}
		}
		if( !YCentered ){
			XCentered = true;
			YCentered = true;
			ZCentered = true;
			transvec_y += 1./( (double) max_count );
			for(unsigned int i=0;i<nbAt;i++){
				// Change reduced coordinates
				coords[i*3+1] += 1./( (double) max_count );
				// Wrapped reduced coordinates
				if( coords[i*3+1] >= 1. || coords[i*3+1] < 0. ) coords[i*3+1] = coords[i*3+1]-floor(coords[i*3+1]);
				// Wrapped normal coordinates
				for(unsigned int k=0;k<3;k++) TempVec[k] = coords[i*3]*H1[k] + coords[i*3+1]*H2[k] + coords[i*3+2]*H3[k];
				// New reduced coordinates
				for(unsigned int k=0;k<3;k++) coords[i*3+k] = TempVec[0]*G1[k] + TempVec[1]*G2[k] + TempVec[2]*G3[k];
				if( coords[i*3+0] < tolRedX || coords[i*3+0] > 1.-tolRedX ) XCentered = false;
				if( coords[i*3+1] < tolRedY || coords[i*3+1] > 1.-tolRedY ) YCentered = false;
				if( coords[i*3+2] < tolRedZ || coords[i*3+2] > 1.-tolRedZ ) ZCentered = false;
			}
		}
		if( !ZCentered ){
			XCentered = true;
			YCentered = true;
			ZCentered = true;
			transvec_z +=1./( (double) max_count );
			for(unsigned int i=0;i<nbAt;i++){
				// Change reduced coordinates
				coords[i*3+2] += 1./( (double) max_count );
				// Wrapped reduced coordinates
				if( coords[i*3+2] >= 1. || coords[i*3+2] < 0. ) coords[i*3+2] = coords[i*3+2]-floor(coords[i*3+2]);
				// Wrapped normal coordinates
				for(unsigned int k=0;k<3;k++) TempVec[k] = coords[i*3]*H1[k] + coords[i*3+1]*H2[k] + coords[i*3+2]*H3[k];
				// New reduced coordinates
				for(unsigned int k=0;k<3;k++) coords[i*3+k] = TempVec[0]*G1[k] + TempVec[1]*G2[k] + TempVec[2]*G3[k];
				if( coords[i*3+0] < tolRedX || coords[i*3+0] > 1.-tolRedX ) XCentered = false;
				if( coords[i*3+1] < tolRedY || coords[i*3+1] > 1.-tolRedY ) YCentered = false;
				if( coords[i*3+2] < tolRedZ || coords[i*3+2] > 1.-tolRedZ ) ZCentered = false;
			}
		}
	count += 1;	
	} // end while
	if( count == max_count ){
		cout << "We didn't find translation" << endl;
		return false;
	}else{
		for(unsigned int k=0;k<3;k++) TempVec[k] = transvec_x*H1[k] + transvec_y*H2[k] + transvec_z*H3[k];
		transvec_x = TempVec[0];
		transvec_y = TempVec[1];
		transvec_z = TempVec[2];
		//cout << "Translation find with vector : " << transvec_x << " " << transvec_y << " " << transvec_z << endl; 
		return true;
	}
}


Eigen::VectorXd make_3d_vector(const double val1, const double val2,
                               const double val3) {
  Eigen::VectorXd result(3);
  result << val1, val2, val3;
  return result;
};

unsigned int getMax(Eigen::VectorXd vec){
        unsigned int max = 0;
        for(unsigned int i=0;i<vec.size();i++){
                if( vec(i) > vec(max) ) max = i;
        }
        return max;
}

gauss::gmm::GaussianMixtureModel fitOptimalGMM(
    const gauss::TrainSet &train_set,
    const vector<size_t> &N_clusters_to_try,
    const size_t &Iterations = 1000) {
  if (N_clusters_to_try.empty()) {
    throw gauss::Error("No clusters to try");
  }
  gauss::gmm::TrainInfo info;
  info.maxIterations = Iterations;

  double MaxDiff = 0.05;
  double old_likelihood, new_likelihood;
  vector<gauss::gmm::Cluster> old_model = vector<gauss::gmm::Cluster>(ExpectationMaximization(train_set, N_clusters_to_try.front(), info, &old_likelihood));
  for (size_t k = 1; k < N_clusters_to_try.size(); ++k) {
    vector<gauss::gmm::Cluster> new_model = vector<gauss::gmm::Cluster>(ExpectationMaximization(train_set, N_clusters_to_try[k], info, &new_likelihood));
    if( fabs(2.*(old_likelihood-new_likelihood)/(old_likelihood+new_likelihood)) < MaxDiff ) break;
    //cout << N_clusters_to_try[k] << " " << 2.*(old_likelihood-new_likelihood)/(old_likelihood+new_likelihood) << endl;
    old_likelihood = new_likelihood;
    old_model = move(new_model);
  }
  return gauss::gmm::GaussianMixtureModel (old_model);
}

gauss::gmm::GaussianMixtureModel fitOptimalGMM_am(
    const gauss::TrainSet train_set,
    const vector<size_t> &N_clusters_to_try,
    const size_t &Iterations = 1000) {
  if (N_clusters_to_try.empty()) {
    throw gauss::Error("No clusters to try");
  }
  gauss::gmm::TrainInfo info;
  info.maxIterations = Iterations;

  double MaxDiff = 0.05;
  double old_likelihood, new_likelihood;
  vector<gauss::gmm::Cluster> old_model =
      vector<gauss::gmm::Cluster>(ExpectationMaximization(
          train_set, N_clusters_to_try.front(), info, &old_likelihood));
  for (size_t k = 1; k < N_clusters_to_try.size(); ++k) {
    vector<gauss::gmm::Cluster> new_model =
        vector<gauss::gmm::Cluster>(ExpectationMaximization(
            train_set, N_clusters_to_try[k], info, &new_likelihood));
    if( fabs(2.*(old_likelihood-new_likelihood)/(old_likelihood+new_likelihood)) < MaxDiff ) break;
    //cout << N_clusters_to_try[k] << " " << 2.*(old_likelihood-new_likelihood)/(old_likelihood+new_likelihood) << endl;
    old_likelihood = new_likelihood;
    old_model = move(new_model);
  }
  return gauss::gmm::GaussianMixtureModel (old_model);
}

int main(int argc, char *argv[])
{
	// Here we use a dump file containing (1) the structure of a given ion (either interface, amorph or crystal), (2) the GB Id, (3) the atomic volume, (4) strain and (5) stress tensor
	// the objective is to return a file for each GB Id with:
	// nx ny nz AmorphousThickness GBSurface Interface1Strain(AndStress)Tensor Interface2Strain(AndStress)Tensor AmorphousStrain(AndStress)Tensor Interface1Strain(AndStress)Invariant Interface2Strain(AndStress)Invariants AmorphousStrain(AndStress)Invariants 
	string InputFilename, OutputFilename, GTFilename;
	string timestep;
	double facClust, tolspangle;
	//unsigned int timestep;

	if( argc == 5 ){
		InputFilename = argv[1];
		GTFilename = argv[2];
		OutputFilename = argv[3];
		timestep = argv[4];
		//istringstream iss_rc(argv[4]);
		//iss_rc >> timestep;
		//istringstream iss_c(argv[4]);
		//iss_c >> tolspangle;
	}else{
		cerr << "Usage: ./AnalyzePolycrystal InputFilename GrainTagsFilename OutputFilename Timestep" << endl;
		cerr << "TODO description" << endl;
		return EXIT_FAILURE;
	}

	tolspangle = 5;
	double tol_sp = cos(tolspangle*M_PI/180.);
	
	// Search the cell vector of each grains
	vector<unsigned int> IdIonsGrain;
	bool ok = read_GTFile(GTFilename, IdIonsGrain);
	if( ok ) cout << "Grain tag file read successfully" << endl;
	else EXIT_FAILURE;
	unsigned int nbGrains = IdIonsGrain.size()/4;
	AtomicSystem MySystem(InputFilename);
	const unsigned int nbAt = MySystem.getNbAtom();
	double *CellVectors = new double[9*nbGrains];
	for(unsigned int i=0;i<nbGrains;i++){
		for(unsigned int dim=0;dim<3;dim++){
			CellVectors[i*9+dim*3] = MySystem.getAtom(IdIonsGrain[i*4+dim+1]).pos.x - MySystem.getAtom(IdIonsGrain[i*4]).pos.x;
			CellVectors[i*9+dim*3+1] = MySystem.getAtom(IdIonsGrain[i*4+dim+1]).pos.y - MySystem.getAtom(IdIonsGrain[i*4]).pos.y;
			CellVectors[i*9+dim*3+2] = MySystem.getAtom(IdIonsGrain[i*4+dim+1]).pos.z - MySystem.getAtom(IdIonsGrain[i*4]).pos.z;
		}
		// verify that the three vectors are almost normal
		double na(0.),nb(0.),nc(0.),sp(0.);
		for(unsigned int dim=0;dim<3;dim++){
			na += pow(CellVectors[i*9+dim],2.);
			nb += pow(CellVectors[i*9+3+dim],2.);
			nc += pow(CellVectors[i*9+6+dim],2.);
		}
		na = sqrt(na);
		nb = sqrt(nb);
		nc = sqrt(nc);
		//cout << "na = " << na << endl;
		//cout << "nb = " << nb << endl;
		//cout << "nc = " << nc << endl;
		// a,b
		for(unsigned int k=0;k<3;k++) sp += CellVectors[i*9+k]*CellVectors[i*9+3+k];
		if( sp/(na*nb) > tol_sp ) cout << "Warning is seems that a and b are not normals in grain " << i+1 << endl;
		//cout << "a,b angle = " << acos(sp/(na*nb))*180./M_PI << endl;
		// a,c
		sp = 0.;
		for(unsigned int k=0;k<3;k++) sp += CellVectors[i*9+k]*CellVectors[i*9+6+k];
		if( sp/(na*nc) > tol_sp ) cout << "Warning is seems that a and c are not normals in grain " << i+1 << endl;
		//cout << "a,c angle = " << acos(sp/(na*nc))*180./M_PI << endl;
		// b,c
		sp = 0.;
		for(unsigned int k=0;k<3;k++) sp += CellVectors[i*9+3+k]*CellVectors[i*9+6+k];
		if( sp/(nb*nc) > tol_sp ) cout << "Warning is seems that b and c are not normals in grain " << i+1 << endl;
		//cout << "b,c angle = " << acos(sp/(nb*nc))*180./M_PI << endl;
	}
		
	string namef1 = "GrainCellVec_";
	string nameext = ".dat";
	std::ofstream f_out_cv(namef1+timestep+nameext);
	for(unsigned int i=0;i<nbGrains;i++){
		for(unsigned int dim1=0;dim1<3;dim1++){
			f_out_cv << CellVectors[i*9+dim1*3];
			for(unsigned int dim2=1;dim2<3;dim2++) f_out_cv << " " << CellVectors[i*9+dim1*3+dim2];
			f_out_cv << endl;
		}
	}
	f_out_cv.close();


	facClust = 0.2;



	unsigned int Struct_GB = 5; //TODO add in argument of exe
	unsigned int Struct_Amorph = 0;
	
	// Get the different aux properties
	unsigned int size_Struct;
	unsigned Struct_ind = MySystem.getAuxIdAndSize("Struct",size_Struct);
	unsigned int size_GBId;
	unsigned GBId_ind = MySystem.getAuxIdAndSize("GBId",size_GBId);
	unsigned int size_AtVol;
	unsigned AtVol_ind = MySystem.getAuxIdAndSize("AtVol",size_AtVol);
	unsigned int size_Stress;
	unsigned Stress_ind = MySystem.getAuxIdAndSize("c_s",size_Stress);
	unsigned int size_AtStrain;
	unsigned AtStrain_ind = MySystem.getAuxIdAndSize("AtomicStrain",size_AtStrain);

	// Search which ions belongs to GB and differentiate the one in amorph and the one in interface
	vector<double> GBId_arr;
	vector<vector<unsigned int>> GBIons;
	cout << "Searching GB ions" << endl;
	// TODO Paralelize?
	for(unsigned int i=0;i<nbAt;i++) {
		bool IsGB = false;
		unsigned int GBOrAmorph_1;
		unsigned int GBOrAmorph_2;
		if( (unsigned int) MySystem.getAux(Struct_ind)[i*size_Struct] == Struct_GB ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
			GBOrAmorph_1 = 0;       
			GBOrAmorph_2 = 0;       
			IsGB = true;
		}else if( (unsigned int) MySystem.getAux(Struct_ind)[i*size_Struct] == Struct_Amorph ){ // Warning here maybe safer to compare fabs() < eps because we compare doubles
			GBOrAmorph_1 = 2;       
			GBOrAmorph_2 = 1;       
			IsGB = true;
		}
		if( IsGB ){
			bool IsAlreadyStored = false;
			for(unsigned int g=0;g<GBId_arr.size();g++){
				if( MySystem.getAux(GBId_ind)[i*size_GBId] == GBId_arr[g] ){ // 1.6 case first row (or if amorph third row)
					GBIons[g*3+GBOrAmorph_1].push_back(i);
					IsAlreadyStored = true;
					break;
				}else if( fabs(MySystem.getAux(GBId_ind)[i*size_GBId] - ( 10.*(GBId_arr[g]-floor(GBId_arr[g])) + ((double) floor(GBId_arr[g])/10.) ) ) < 1e-5 ){ // 6.1 case second raw (or third raw in amorph case) # TODO warning this method should not work if a grain index is higher than 10
					GBIons[g*3+1+GBOrAmorph_2].push_back(i);
					IsAlreadyStored = true;
					break;
				}
			}
			if( !IsAlreadyStored ){
				GBId_arr.push_back(MySystem.getAux(GBId_ind)[i*size_GBId]);
				GBIons.push_back(vector<unsigned int>());
				GBIons.push_back(vector<unsigned int>());
				GBIons.push_back(vector<unsigned int>());
				GBIons[(GBId_arr.size()-1)*3+GBOrAmorph_1].push_back(i); //not sure about the minus one (to verify)
			}
		}
	}
	cout << "Done" << endl;
	// for each GBId, translate the ions to have the clusters not separated by box boundaries
	double tolBound = 2.; // in A, for an ions to be consider in contact with box boundaries
	double tolRedX = 2./MySystem.getH1()[0]; 
	double tolRedY = 2./MySystem.getH2()[1]; 
	double tolRedZ = 2./MySystem.getH3()[2]; 
	vector<Eigen::VectorXd> coords;
	vector<Eigen::VectorXd> coords_am;
	double *TransVec = new double[GBId_arr.size()*6];
	double *TempVec = new double[3];
	unsigned int nbMinIonsInGB = 100;
	std::vector<std::size_t> NbClustToTest;
	std::size_t NMax = 3;
  	for(unsigned int i=1;i<NMax+1;i++) NbClustToTest.push_back(i);
	const size_t clusters_size = 2;
	vector<vector<double>> MeanAndCovarClust_temp;
	vector<vector<unsigned int>> GBIons_temp;
	vector<vector<double>> MeanAndCovarClust_Final;
	vector<double> GBId_arr_Final;
	vector<vector<unsigned int>> GBIons_Final;
	//double facClust = 0.1;
	unsigned int nbGBAnalyzed = 0;
	unsigned int current_nbGBAnalyzed;
	MathTools MT;
	for(unsigned int g=0;g<GBId_arr.size();g++){
		cout << "Treating " << GBId_arr[g] << " GB" << endl;
		if( GBIons[g*3].size() < nbMinIonsInGB || GBIons[g*3+1].size() < nbMinIonsInGB ){
			cout << "skipping, not enough ions" << endl;
			continue;
		}

		MeanAndCovarClust_temp.clear();
		MeanAndCovarClust_temp.push_back(vector<double>());
		MeanAndCovarClust_temp.push_back(vector<double>());
		bool XCentered = true;
		bool YCentered = true;
		bool ZCentered = true;
		for(unsigned int i=0;i<GBIons_temp.size();i++){
			GBIons_temp[i].clear();
			vector<unsigned int>().swap(GBIons_temp[i]);
		}
		GBIons_temp.clear();
		vector<vector<unsigned int>>().swap(GBIons_temp);
		current_nbGBAnalyzed = 0;
		double *coord_for_wrap = new double[3*(GBIons[g*3].size()+GBIons[g*3+1].size())];
		for(unsigned int gbid=0;gbid<2;gbid++){
			for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
				coord_for_wrap[(gbid*GBIons[g*3].size()+i)*3] = MySystem.getWrappedPos(GBIons[g*3+gbid][i]).x;
				coord_for_wrap[(gbid*GBIons[g*3].size()+i)*3+1] = MySystem.getWrappedPos(GBIons[g*3+gbid][i]).y;
				coord_for_wrap[(gbid*GBIons[g*3].size()+i)*3+2] = MySystem.getWrappedPos(GBIons[g*3+gbid][i]).z;
			}
		}
		double xtemp,ytemp,ztemp;
		if( ComputeTransVec(coord_for_wrap, GBIons[g*3].size()+GBIons[g*3+1].size(), MySystem.getH1(), MySystem.getH2(), MySystem.getH3(), MySystem.getG1(), MySystem.getG2(), MySystem.getG3(),xtemp, ytemp, ztemp) ){
			TransVec[g*3] = xtemp;
			TransVec[g*3+1] = ytemp;
			TransVec[g*3+2] = ztemp;
		}else{
			cout << "issue" << endl;
		}
		delete[] coord_for_wrap;
		for(unsigned int gbid=0;gbid<2;gbid++){
			vector<Eigen::VectorXd>().swap(coords);
			for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++){
				double xpos = MySystem.getWrappedPos(GBIons[g*3+gbid][i]).x+TransVec[g*3];
				double ypos = MySystem.getWrappedPos(GBIons[g*3+gbid][i]).y+TransVec[g*3+1];
				double zpos = MySystem.getWrappedPos(GBIons[g*3+gbid][i]).z+TransVec[g*3+2];
				double xpos_red = xpos*MySystem.getG1()[0]+ypos*MySystem.getG2()[0]+zpos*MySystem.getG3()[0];
				double ypos_red = xpos*MySystem.getG1()[1]+ypos*MySystem.getG2()[1]+zpos*MySystem.getG3()[1];
				double zpos_red = xpos*MySystem.getG1()[2]+ypos*MySystem.getG2()[2]+zpos*MySystem.getG3()[2];
				if( xpos_red >= 1. || xpos_red < 0. ) xpos_red = xpos_red-floor(xpos_red);
				if( ypos_red >= 1. || xpos_red < 0. ) ypos_red = ypos_red-floor(ypos_red);
				if( zpos_red >= 1. || xpos_red < 0. ) zpos_red = zpos_red-floor(zpos_red);
				xpos = xpos_red*MySystem.getH1()[0]+ypos_red*MySystem.getH2()[0]+zpos_red*MySystem.getH3()[0];
				ypos = xpos_red*MySystem.getH1()[1]+ypos_red*MySystem.getH2()[1]+zpos_red*MySystem.getH3()[1];
				zpos = xpos_red*MySystem.getH1()[2]+ypos_red*MySystem.getH2()[2]+zpos_red*MySystem.getH3()[2];
				coords.push_back(make_3d_vector(xpos,ypos,zpos));
			}
			// find the GMM
			//vector<gauss::gmm::Cluster> clusters = gauss::gmm::ExpectationMaximization(samples, clusters_size);
  			//gauss::gmm::GaussianMixtureModel gmm_model(clusters);
			gauss::TrainSet set(coords);
			//if( fabs(GBId_arr[g]-2.3) < 1e-2 ){ // here is a bug when NbClustToTest is higher than 3
			//	cout << "stop" << endl;
			//}
			gauss::gmm::GaussianMixtureModel new_gmm_model = fitOptimalGMM(set,NbClustToTest);
			unsigned int nbClustFinal = new_gmm_model.getClusters().size();
			Eigen::MatrixXd *Covars = new Eigen::MatrixXd[nbClustFinal];
			Eigen::VectorXd *Means = new Eigen::VectorXd[nbClustFinal];
			Eigen::EigenSolver<Eigen::MatrixXd> Solv;
			vector<double> *Ids = new vector<double>[nbClustFinal];
			for(unsigned int i=0;i<GBIons[g*3+gbid].size();i++) Ids[getMax(new_gmm_model.Classify(coords[i]))].push_back(GBIons[g*3+gbid][i]);
			
			for(unsigned int i=0;i<nbClustFinal;i++){
				if( Ids[i].size() >  nbMinIonsInGB ){
			        	Covars[i] = new_gmm_model.getClusters()[i].distribution->getCovariance();
			        	Means[i] = new_gmm_model.getClusters()[i].distribution->getMean();
			        	Solv.compute(Covars[i]);
					double ev1 = Solv.eigenvalues().col(0)[0].real();
					double ev2 = Solv.eigenvalues().col(0)[1].real();
					double ev3 = Solv.eigenvalues().col(0)[2].real();
			        	//cout << "for cluster " << i << endl;
			        	//cout << "eigenvalue " << Solv.eigenvalues().col(0)[0].real() << ", eigenvector " << Solv.eigenvectors().col(0)[0].real() << "  " << Solv.eigenvectors().col(0)[1].real() << "  " << Solv.eigenvectors().col(0)[2].real() << endl;
			        	//cout << "eigenvalue " << Solv.eigenvalues().col(0)[1].real() << ", eigenvector " << Solv.eigenvectors().col(1)[0].real() << "  " << Solv.eigenvectors().col(1)[1].real() << "  " << Solv.eigenvectors().col(1)[2].real() << endl;
			        	//cout << "eigenvalue " << Solv.eigenvalues().col(0)[2].real() << ", eigenvector " << Solv.eigenvectors().col(2)[0].real() << "  " << Solv.eigenvectors().col(2)[1].real() << "  " << Solv.eigenvectors().col(2)[2].real() << endl;
					bool tostore = false;
					if( ( ev1 < facClust*ev2 ) && ( ev1 < facClust*ev3 ) ){ // todo mean
						MeanAndCovarClust_temp[gbid].push_back(Solv.eigenvectors().col(0)[0].real());
						MeanAndCovarClust_temp[gbid].push_back(Solv.eigenvectors().col(0)[1].real());
						MeanAndCovarClust_temp[gbid].push_back(Solv.eigenvectors().col(0)[2].real());
						tostore = true;
					}else if( ( ev2 < facClust*ev1 ) && ( ev2 < facClust*ev3 ) ){
						MeanAndCovarClust_temp[gbid].push_back(Solv.eigenvectors().col(1)[0].real());
						MeanAndCovarClust_temp[gbid].push_back(Solv.eigenvectors().col(1)[1].real());
						MeanAndCovarClust_temp[gbid].push_back(Solv.eigenvectors().col(1)[2].real());
						tostore = true;
					}else if( ( ev3 < facClust*ev2 ) && ( ev3 < facClust*ev1 ) ){
						MeanAndCovarClust_temp[gbid].push_back(Solv.eigenvectors().col(2)[0].real());
						MeanAndCovarClust_temp[gbid].push_back(Solv.eigenvectors().col(2)[1].real());
						MeanAndCovarClust_temp[gbid].push_back(Solv.eigenvectors().col(2)[2].real());
						tostore = true;
					} 
					if( tostore && gbid == 0 ){
						MeanAndCovarClust_temp[gbid].push_back(Means[i](0));
						MeanAndCovarClust_temp[gbid].push_back(Means[i](1));
						MeanAndCovarClust_temp[gbid].push_back(Means[i](2));
						GBIons_temp.push_back(vector<unsigned int>());
						for(unsigned int n=0;n<Ids[i].size();n++) GBIons_temp[GBIons_temp.size()-1].push_back(Ids[i][n]);
					}else if( tostore && gbid == 1 ){
						MeanAndCovarClust_temp[gbid].push_back(Means[i](0));
						MeanAndCovarClust_temp[gbid].push_back(Means[i](1));
						MeanAndCovarClust_temp[gbid].push_back(Means[i](2));
						// search if there is corresponding GB 1 with almost the same normal
						bool isAssociated = false;
						//double tol_sp = 1-(5e-3);
						double n0, n1, sp, fullsp;
						vector<unsigned int> indexes;
						vector<double> scalarprod;
						unsigned int index_clust_ass;
						for(unsigned int n=0;n<GBIons_temp.size();n++){
							sp = 0.;
							n0 = 0.;
							n1 = 0.;
							for(unsigned int dim=0;dim<3;dim++){
								sp += MeanAndCovarClust_temp[0][6*n+dim]*MeanAndCovarClust_temp[1][((MeanAndCovarClust_temp[1].size()/6)-1)*6+dim];
								n0 += pow(MeanAndCovarClust_temp[0][6*n+dim],2.);
								n1 += pow(MeanAndCovarClust_temp[1][((MeanAndCovarClust_temp[1].size()/6)-1)*6+dim],2.);
							}
							n0 = sqrt(n0);
							n1 = sqrt(n1);
							fullsp = fabs(sp/(n1*n0));
							if( fullsp > tol_sp ){
								isAssociated = true;
								scalarprod.push_back(fullsp);
								indexes.push_back(n);
							}
						}
						if( isAssociated ){
							index_clust_ass = indexes[MT.max(scalarprod)];
							nbGBAnalyzed += 1;
							current_nbGBAnalyzed += 1;
							// Store normals and mean positions for this GB and average of normals
							MeanAndCovarClust_Final.push_back(vector<double>());
							for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back(MeanAndCovarClust_temp[0][6*index_clust_ass+dim]);
							for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back(MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]);
							double x1 = MeanAndCovarClust_temp[0][6*index_clust_ass+0];
							double y1 = MeanAndCovarClust_temp[0][6*index_clust_ass+1];
							double z1 = MeanAndCovarClust_temp[0][6*index_clust_ass+2];
							double x2 = MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+0];
							double y2 = MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+1];
							double z2 = MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+2];
							if( ( fabs(x1) > 1e-1 && fabs(x2) > 1e-1 && x1*x2 < 0. ) ||  ( fabs(y1) > 1e-1 && fabs(y2) > 1e-1 && y1*y2 < 0. ) || ( fabs(z1) > 1e-1 && fabs(z2) > 1e-1 && z1*z2 < 0. ) ){
								double norm = 0.;
								for(unsigned int dim=0;dim<3;dim++) norm += pow((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]-MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/2.,2.);
								norm = sqrt(norm);
								for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]-MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/(2.*norm));
							}else{
								double norm = 0.;
								for(unsigned int dim=0;dim<3;dim++) norm += pow((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]+MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/2.,2.);
								for(unsigned int dim=0;dim<3;dim++) MeanAndCovarClust_Final[nbGBAnalyzed-1].push_back((MeanAndCovarClust_temp[1][6*((MeanAndCovarClust_temp[1].size()/6)-1)+dim]+MeanAndCovarClust_temp[0][6*index_clust_ass+dim])/(2.*norm));
							}
							// Store ions ID for this GB
							GBIons_Final.push_back(vector<unsigned int>());
							GBIons_Final.push_back(vector<unsigned int>());
							GBIons_Final.push_back(vector<unsigned int>());
							for(unsigned int nb=0;nb<GBIons_temp[index_clust_ass].size();nb++) GBIons_Final[(nbGBAnalyzed-1)*3].push_back(GBIons_temp[index_clust_ass][nb]);
							for(unsigned int nb=0;nb<Ids[i].size();nb++) GBIons_Final[(nbGBAnalyzed-1)*3+1].push_back(Ids[i][nb]);

							// Stored ids of grain
							GBId_arr_Final.push_back(GBId_arr[g]);

							// delete this cluster
							for(unsigned int dim=0;dim<6;dim++) MeanAndCovarClust_temp[0].erase(MeanAndCovarClust_temp[0].begin()+6*index_clust_ass);
							GBIons_temp.erase(GBIons_temp.begin()+index_clust_ass);
						}
					}
				}
			}
			delete[] Ids;
			delete[] Covars;
			delete[] Means;
		} // end gbid loop
		cout << current_nbGBAnalyzed << " interfaces analyzed for this GB" << endl;
		// Amorph phase
		if( current_nbGBAnalyzed != 0 && GBIons[g*3+2].size() > nbMinIonsInGB ){
			vector<Eigen::VectorXd>().swap(coords_am);
			for(unsigned int i=0;i<GBIons[g*3+2].size();i++){
				// Compute reduced coordinates
				double xpos = MySystem.getWrappedPos(GBIons[g*3+2][i]).x*MySystem.getG1()[0]+MySystem.getWrappedPos(GBIons[g*3+2][i]).y*MySystem.getG2()[0]+MySystem.getWrappedPos(GBIons[g*3+2][i]).z*MySystem.getG3()[0];
				double ypos = MySystem.getWrappedPos(GBIons[g*3+2][i]).x*MySystem.getG1()[1]+MySystem.getWrappedPos(GBIons[g*3+2][i]).y*MySystem.getG2()[1]+MySystem.getWrappedPos(GBIons[g*3+2][i]).z*MySystem.getG3()[1];
				double zpos = MySystem.getWrappedPos(GBIons[g*3+2][i]).x*MySystem.getG1()[2]+MySystem.getWrappedPos(GBIons[g*3+2][i]).y*MySystem.getG2()[2]+MySystem.getWrappedPos(GBIons[g*3+2][i]).z*MySystem.getG3()[2];
				coords_am.push_back(make_3d_vector(xpos,ypos,zpos));
				if( coords_am[i](0) < tolRedX || coords_am[i](0) > 1.-tolRedX ) XCentered = false;
				if( coords_am[i](1) < tolRedY || coords_am[i](1) > 1.-tolRedY ) YCentered = false;
				if( coords_am[i](2) < tolRedZ || coords_am[i](2) > 1.-tolRedZ ) ZCentered = false;
			}
			unsigned int max_count = 10;
			unsigned int count = 0;
			while( ( !XCentered || !YCentered || !ZCentered ) && ( count < max_count ) ){
				if( !XCentered ){
					XCentered = true;
					YCentered = true;
					ZCentered = true;
					for(unsigned int i=0;i<GBIons[g*3+2].size();i++){
						// Change reduced coordinates
						coords_am[i](0) += 1./( (double) max_count );
						// Wrapped reduced coordinates
						if( coords_am[i](0) >= 1. || coords_am[i](0) < 0. ) coords_am[i](0) = coords_am[i](0)-floor(coords_am[i](0));
						// Wrapped normal coordinates
						for(unsigned int k=0;k<3;k++) TempVec[k] = coords_am[i](0)*MySystem.getH1()[k] + coords_am[i](1)*MySystem.getH2()[k] + coords_am[i](2)*MySystem.getH3()[k];
						// New reduced coordinates
						for(unsigned int k=0;k<3;k++) coords_am[i](k) = TempVec[0]*MySystem.getG1()[k] + TempVec[1]*MySystem.getG2()[k] + TempVec[2]*MySystem.getG3()[k];
						if( coords_am[i](0) < tolRedX || coords_am[i](0) > 1.-tolRedX ) XCentered = false;
						if( coords_am[i](1) < tolRedY || coords_am[i](1) > 1.-tolRedY ) YCentered = false;
						if( coords_am[i](2) < tolRedZ || coords_am[i](2) > 1.-tolRedZ ) ZCentered = false;
					}
				}
				if( !YCentered ){
					XCentered = true;
					YCentered = true;
					ZCentered = true;
					for(unsigned int i=0;i<GBIons[g*3+2].size();i++){
						// Change reduced coordinates
						coords_am[i](1) += 1./( (double) max_count );
						// Wrapped reduced coordinates
						if( coords_am[i](1) >= 1. || coords_am[i](1) < 0. ) coords_am[i](1) = coords_am[i](1)-floor(coords_am[i](1));
						// Wrapped normal coordinates
						for(unsigned int k=0;k<3;k++) TempVec[k] = coords_am[i](0)*MySystem.getH1()[k] + coords_am[i](1)*MySystem.getH2()[k] + coords_am[i](2)*MySystem.getH3()[k];
						// New reduced coordinates
						for(unsigned int k=0;k<3;k++) coords_am[i](k) = TempVec[0]*MySystem.getG1()[k] + TempVec[1]*MySystem.getG2()[k] + TempVec[2]*MySystem.getG3()[k];
						if( coords_am[i](0) < tolRedX || coords_am[i](0) > 1.-tolRedX ) XCentered = false;
						if( coords_am[i](1) < tolRedY || coords_am[i](1) > 1.-tolRedY ) YCentered = false;
						if( coords_am[i](2) < tolRedZ || coords_am[i](2) > 1.-tolRedZ ) ZCentered = false;
					}
				}
				if( !ZCentered ){
					XCentered = true;
					YCentered = true;
					ZCentered = true;
					for(unsigned int i=0;i<GBIons[g*3+2].size();i++){
						// Change reduced coordinates
						coords_am[i](2) += 1./( (double) max_count );
						// Wrapped reduced coordinates
						if( coords_am[i](2) >= 1. || coords_am[i](2) < 0. ) coords_am[i](2) = coords_am[i](2)-floor(coords_am[i](2));
						// Wrapped normal coordinates
						for(unsigned int k=0;k<3;k++) TempVec[k] = coords_am[i](0)*MySystem.getH1()[k] + coords_am[i](1)*MySystem.getH2()[k] + coords_am[i](2)*MySystem.getH3()[k];
						// New reduced coordinates
						for(unsigned int k=0;k<3;k++) coords_am[i](k) = TempVec[0]*MySystem.getG1()[k] + TempVec[1]*MySystem.getG2()[k] + TempVec[2]*MySystem.getG3()[k];
						if( coords_am[i](0) < tolRedX || coords_am[i](0) > 1.-tolRedX ) XCentered = false;
						if( coords_am[i](1) < tolRedY || coords_am[i](1) > 1.-tolRedY ) YCentered = false;
						if( coords_am[i](2) < tolRedZ || coords_am[i](2) > 1.-tolRedZ ) ZCentered = false;
					}
				}
			count += 1;	
			} // end while
			if( count == max_count ) cout << "We didn't find translation" << endl;
			//else{
			//	for(unsigned int k=0;k<3;k++) TempVec[k] = TransVec[g*6+(2*3)]*MySystem.getH1()[k] + TransVec[g*6+(2*3)+1]*MySystem.getH2()[k] + TransVec[g*6+(2*3)+2]*MySystem.getH3()[k];
			//	for(unsigned int k=0;k<3;k++) TransVec[g*6+(2*3)+k] = TempVec[k]; 
			//	//cout << "Translation find with vector : " << TransVec[g*6+(gbid*3)] << " " << TransVec[g*6+(gbid*3)+1] << " " << TransVec[g*6+(gbid*3)+2] << endl; 
			//}

			// find the GMM
			for(unsigned int i=0;i<GBIons[g*3+2].size();i++){
				for(unsigned int k=0;k<3;k++) TempVec[k] = coords_am[i](0)*MySystem.getH1()[k] + coords_am[i](1)*MySystem.getH2()[k] + coords_am[i](2)*MySystem.getH3()[k];
				for(unsigned int k=0;k<3;k++) coords_am[i](k) = TempVec[k]; 
			}
			gauss::TrainSet set_am(coords_am);
			gauss::gmm::GaussianMixtureModel new_gmm_model_am = fitOptimalGMM_am(set_am,NbClustToTest);
			unsigned int nbClustFinal = new_gmm_model_am.getClusters().size();
			Eigen::MatrixXd *Covars = new Eigen::MatrixXd[nbClustFinal];
			Eigen::VectorXd *Means = new Eigen::VectorXd[nbClustFinal];
			Eigen::EigenSolver<Eigen::MatrixXd> Solv;
			vector<double> *Ids = new vector<double>[nbClustFinal];
			double *normal_vec = new double[3*nbClustFinal];
			for(unsigned int i=0;i<GBIons[g*3+2].size();i++) Ids[getMax(new_gmm_model_am.Classify(coords_am[i]))].push_back(GBIons[g*3+2][i]);
			for(unsigned int i=0;i<nbClustFinal;i++){
			        Covars[i] = new_gmm_model_am.getClusters()[i].distribution->getCovariance();
			        Means[i] = new_gmm_model_am.getClusters()[i].distribution->getMean();
			        Solv.compute(Covars[i]);
				double ev1 = Solv.eigenvalues().col(0)[0].real();
				double ev2 = Solv.eigenvalues().col(0)[1].real();
				double ev3 = Solv.eigenvalues().col(0)[2].real();
				if( ( ev1 < ev2 ) && ( ev1 < ev3 ) ){ 
					for(unsigned int dim=0;dim<3;dim++) normal_vec[i*3+dim] = Solv.eigenvectors().col(0)[dim].real();
				}else if( ( ev2 < ev1 ) && ( ev2 < ev3 ) ){
					for(unsigned int dim=0;dim<3;dim++) normal_vec[i*3+dim] = Solv.eigenvectors().col(1)[dim].real();
				}else if( ( ev3 < ev2 ) && ( ev3 < ev1 ) ){
					for(unsigned int dim=0;dim<3;dim++) normal_vec[i*3+dim] = Solv.eigenvectors().col(2)[dim].real();
				}
			}
			double n0, n1, sp, fullsp;
			vector<unsigned int> indexes;
			vector<unsigned int> already_used;
			vector<double> scalarprod;
			unsigned int index_clust_ass;  
			for(unsigned int n=0;n<current_nbGBAnalyzed;n++){
				scalarprod.clear();
				indexes.clear();
				for(unsigned int i=0;i<nbClustFinal;i++){
					bool already = false;
					for(unsigned aa=0;aa<already_used.size();aa++){
						if( i == already_used[aa] ){
							already = true;
							break;
						}
					}
					if( !already ){
						sp = 0.;
						n0 = 0.;
						n1 = 0.;
						for(unsigned int dim=0;dim<3;dim++){
							sp += normal_vec[i*3+dim]*MeanAndCovarClust_Final[nbGBAnalyzed-1-n][12+dim];
							n0 += pow(MeanAndCovarClust_Final[nbGBAnalyzed-1-n][12+dim],2.);
							n1 += pow(normal_vec[i*3+dim],2.);
						}
						n0 = sqrt(n0);
						n1 = sqrt(n1);
						fullsp = fabs(sp/(n1*n0));
						scalarprod.push_back(fullsp);
						indexes.push_back(i);
					}
				}
				index_clust_ass = indexes[MT.max(scalarprod)];
				// Store ions
				for(unsigned int nb=0;nb<Ids[index_clust_ass].size();nb++) GBIons_Final[(nbGBAnalyzed-1-n)*3+2].push_back(Ids[index_clust_ass][nb]);
			}
			delete[] normal_vec;
			delete[] Covars;
			delete[] Means;
			delete[] Ids;
		}// end amorph
	}

	double *MeanStrain_1 = new double[8];
	double *MeanStrain_2 = new double[8];
	double *MeanStrain_am = new double[8];
	long double *MeanStress_1 = new long double[6];
	long double *MeanStress_2 = new long double[6];
	long double *MeanStress_am = new long double[6];

	std::ofstream f_out_res(OutputFilename);
	f_out_res << "nx(1) ny(2) nz(3) G1(4) G2(5) AmThick(6) Eps_xx_1(7) Eps_yy_1(8) Eps_zz_1(9) Eps_xy_1(10) Eps_xz_1(11) Eps_yz_1(12) ShearInv_1(13) SphInv_1(14) Eps_xx_2(15) Eps_yy_2(16) Eps_zz_2(17) Eps_xy_2(18) Eps_xz_2(19) Eps_yz_2(20) ShearInv_2(21) SphInv_2(22) Eps_xx_am(23) Eps_yy_am(24) Eps_zz_am(25) Eps_xy_am(26) Eps_xz_am(27) Eps_yz_am(28) ShearInv_am(29) SphInv_am(30) sigma_xx_1(31) sigma_yy_1(32) sigma_zz_1(33) sigma_xy_1(34) sigma_xz_1(35) sigma_yz_1(36) sigma_xx_2(37) sigma_yy_2(38) sigma_zz_2(39) sigma_xy_2(40) sigma_xz_2(41) sigma_yz_2(42) sigma_xx_am(43) sigma_yy_am(44) sigma_zz_am(45) sigma_xy_am(46) sigma_xz_am(47) sigma_yz_am(48)" << endl;
	for(unsigned int i=0;i<nbGBAnalyzed;i++){
		//cout << "GB " << i << " " << GBId_arr_Final[i] << endl;
		//cout << "N1 = " << MeanAndCovarClust_Final[i][0] << " " << MeanAndCovarClust_Final[i][1] << " " << MeanAndCovarClust_Final[i][2] << endl; 
		//cout << "N2 = " << MeanAndCovarClust_Final[i][6] << " " << MeanAndCovarClust_Final[i][7] << " " << MeanAndCovarClust_Final[i][8] << endl; 
		//cout << "Average normal = "  << MeanAndCovarClust_Final[i][12] << " " << MeanAndCovarClust_Final[i][13] << " " << MeanAndCovarClust_Final[i][14] << endl; 
		//cout << "M1 = " << MeanAndCovarClust_Final[i][3] << " " << MeanAndCovarClust_Final[i][4] << " " << MeanAndCovarClust_Final[i][5] << endl; 
		//cout << "M2 = " << MeanAndCovarClust_Final[i][9] << " " << MeanAndCovarClust_Final[i][10] << " " << MeanAndCovarClust_Final[i][11] << endl; 
		double d1 = 0.;
		double d2 = 0.;
		for(unsigned int dim=0;dim<3;dim++){
			//MeanAndCovarClust_Final[i][3+dim] -= TransVec[i*6+dim];
			//MeanAndCovarClust_Final[i][9+dim] -= TransVec[i*6+3+dim];
			d1 += MeanAndCovarClust_Final[i][3+dim]*MeanAndCovarClust_Final[i][12+dim];
			d2 += MeanAndCovarClust_Final[i][9+dim]*MeanAndCovarClust_Final[i][12+dim];
		}
		double am_thick;
		if( d1 > d2 ) am_thick = d1-d2;
		else am_thick = d2-d1;
		//cout << "Amorphous thickness : " << am_thick << endl;
		for(unsigned int ii=0;ii<size_AtStrain;ii++){
			MeanStrain_1[ii] = 0.;
			MeanStrain_2[ii] = 0.;
			MeanStrain_am[ii] = 0.;
		}
		for(unsigned int ii=0;ii<size_Stress;ii++){
			MeanStress_1[ii] = 0.;
			MeanStress_2[ii] = 0.;
			MeanStress_am[ii] = 0.;
		}
		// compute average strain and stress for interface 1
		double temp_vol = 0.;
		for(unsigned int n=0;n<GBIons_Final[i*3].size();n++){
			for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_1[ii] += MySystem.getAux(AtStrain_ind)[GBIons_Final[i*3][n]*size_AtStrain+ii];
			for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_1[ii] += (long double) MySystem.getAux(Stress_ind)[GBIons_Final[i*3][n]*size_Stress+ii];
		       	temp_vol += MySystem.getAux(AtVol_ind)[GBIons_Final[i*3][n]*size_AtVol];
		}
		for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_1[ii] /= temp_vol;
		temp_vol = 0.;
		for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_1[ii] /= GBIons_Final[i*3].size();
		// compute average strain and stress for interface 2
		for(unsigned int n=0;n<GBIons_Final[i*3+1].size();n++){
			for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_2[ii] += MySystem.getAux(AtStrain_ind)[GBIons_Final[i*3+1][n]*size_AtStrain+ii];
			for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_2[ii] += (long double) MySystem.getAux(Stress_ind)[GBIons_Final[i*3+1][n]*size_Stress+ii];
		        temp_vol += MySystem.getAux(AtVol_ind)[GBIons_Final[i*3+1][n]*size_AtVol];
		}
		for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_2[ii] /= temp_vol;
		temp_vol = 0.;
		for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_2[ii] /= GBIons_Final[i*3+1].size();
		// compute average strain and stress for amorph 
		for(unsigned int n=0;n<GBIons_Final[i*3+2].size();n++){
			for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_am[ii] += MySystem.getAux(AtStrain_ind)[GBIons_Final[i*3+2][n]*size_AtStrain+ii];
			for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_am[ii] += (long double) MySystem.getAux(Stress_ind)[GBIons_Final[i*3+2][n]*size_Stress+ii];
		        temp_vol += MySystem.getAux(AtVol_ind)[GBIons_Final[i*3+2][n]*size_AtVol];
		}
		for(unsigned int ii=0;ii<size_Stress;ii++) MeanStress_am[ii] /= temp_vol;
		for(unsigned int ii=0;ii<size_AtStrain;ii++) MeanStrain_am[ii] /= GBIons_Final[i*3+2].size();

		f_out_res << MeanAndCovarClust_Final[i][12] << " " << MeanAndCovarClust_Final[i][13] << " " << MeanAndCovarClust_Final[i][14] << " " << GBId_arr_Final[i] - (GBId_arr_Final[i]-floor(GBId_arr_Final[i])) << " " << (int) 10.*(GBId_arr_Final[i]-floor(GBId_arr_Final[i])) << " " << am_thick ;
		for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_1[ii]; 
		for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_2[ii]; 
		for(unsigned int ii=0;ii<8;ii++) f_out_res << " " << MeanStrain_am[ii]; 
		for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_1[ii]; 
		for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_2[ii]; 
		for(unsigned int ii=0;ii<6;ii++) f_out_res << " " << MeanStress_am[ii]; 
		f_out_res << endl;


		string path_1 = "./AtomicSystems_";
		string path_2 = "/";
		string fullpath = path_1+timestep+path_2;
		std::filesystem::create_directory(fullpath);
		//_mkdir(fullpath.c_str());
		string name = "TreatedGB_";
		string ext = ".cfg";
		string und = "_";
		auto ii = to_string(i);
		string fullname = fullpath+name+ii+ext;
		std::ofstream f_out(fullname);
		f_out << "ITEM: TIMESTEP" << endl;
		f_out << timestep << endl;
		f_out << "ITEM: NUMBER OF ATOMS" << endl;
		f_out << GBIons_Final[i*3].size()+GBIons_Final[i*3+1].size()+GBIons_Final[i*3+2].size() << endl;
		double arr[4] = {0.,MySystem.getH2()[0],MySystem.getH3()[0],MySystem.getH2()[0]+MySystem.getH3()[0]};
                double arr_2[2] = {0.,MySystem.getH3()[1]};
		f_out << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n" << MT.min(arr,4) << "\t" << MySystem.getH1()[0]+MT.max(arr,4) << "\t" << MySystem.getH2()[0] << "\n" << MT.min(arr_2,2) << "\t" << MySystem.getH2()[1]+MT.max(arr_2,2) << "\t" << MySystem.getH3()[0] << "\n0\t" << MySystem.getH3()[2] << "\t" << MySystem.getH3()[1] << "\n";
		f_out << "ITEM: ATOMS id x y z clustId" << endl;
		for(unsigned int n=0;n<GBIons_Final[i*3].size();n++){
			f_out << n+1 << " " << MySystem.getWrappedPos(GBIons_Final[i*3][n]).x << " " << MySystem.getWrappedPos(GBIons_Final[i*3][n]).y << " " << MySystem.getWrappedPos(GBIons_Final[i*3][n]).z << " " << 0 << endl;
		}
		for(unsigned int n=0;n<GBIons_Final[i*3+1].size();n++){
			f_out << n+1+GBIons_Final[i*3].size() << " " << MySystem.getWrappedPos(GBIons_Final[i*3+1][n]).x << " " << MySystem.getWrappedPos(GBIons_Final[i*3+1][n]).y << " " << MySystem.getWrappedPos(GBIons_Final[i*3+1][n]).z << " " << 1 << endl;
		}
		for(unsigned int n=0;n<GBIons_Final[i*3+2].size();n++){
			f_out << n+1+GBIons_Final[i*3].size()+GBIons_Final[i*3+1].size() << " " << MySystem.getWrappedPos(GBIons_Final[i*3+2][n]).x << " " << MySystem.getWrappedPos(GBIons_Final[i*3+2][n]).y << " " << MySystem.getWrappedPos(GBIons_Final[i*3+2][n]).z << " " << 2 << endl;
		}
		f_out.close();

	}
	f_out_res.close();
	
	delete[] TransVec;
	delete[] TempVec;
	delete[] MeanStrain_1;
	delete[] MeanStrain_2;
	delete[] MeanStrain_am;
	delete[] MeanStress_1;
	delete[] MeanStress_2;
	delete[] MeanStress_am;
	//cout << "Compute GB plane normals" << endl;
	//// Compute the normal plane to each GB using only interface ions (maybe compute here stress and strain for interfaces)
	//// Pb: multiple planes for a given GBId (i.e. an other one formed by the periodic replica of the grain)
	//double *Normals = new double[GBId_arr.size()*4];
	//MathTools MT;
	//vector<vector<double>> coords;

	//double afit,bfit,cfit;
	//for(unsigned int g=0;g<GBId_arr.size();g++){
	//	coords.clear();
	//	for(unsigned int i=0;i<GBIons[g*2].size();i++){
	//		coords.push_back(vector<double>());
	//		coords[i].push_back(MySystem.getWrappedPos(GBIons[g*2][i]).x);
	//		coords[i].push_back(MySystem.getWrappedPos(GBIons[g*2][i]).y);
	//		coords[i].push_back(MySystem.getWrappedPos(GBIons[g*2][i]).z);
	//	}
	//	MT.plane_fit(coords,afit,bfit,cfit);
	//	// Normalize
	//	double normFac = sqrt(pow(afit,2.)+pow(bfit,2.)+1);
	//	Normals[g*4] = afit/normFac;
	//	Normals[g*4+1] = bfit/normFac;
	//	Normals[g*4+2] = -1./normFac;
	//	Normals[g*4+3] = cfit/normFac;
	//	cout << "GBId : " << GBId_arr[g] << ", plane coeff are : " << Normals[g*4] << ", " << Normals[g*4+1] << ", " << Normals[g*4+2] << ", " << Normals[g*4+3] << ", nb of ion considered : " << GBIons[g*2].size() << endl; 
	//}

	//cout << "Done" << endl;

	//// z = ax+by+c <=> apx+bpy+cpz+dp=0 and a,b,c normed
	//// apx+bpy+cpax+cpbby+ccp+dp=0
	////Points[count].push_back(((double) x));
	////Points[count].push_back(((double) y));
	////Points[count].push_back((a*(double) x)+(b*(double) y)+c+(((double) (rand()%100) - 50.)/10));


	//delete[] Normals;
	return 0;
}
