
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "Descriptors.h"
#include <Displays.h>

using namespace std;
using Eigen::Index;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 6 ){
		cerr << "Usage: ./FitAndSaveGMM InputDirectory (NameOfDescriptor) MaxDimension MaxCorrelation minScore yaml_filename" << endl;
		cerr << "InputDirectory should contains subdirectories having the names of the labels and containing either: " << endl;
	        cerr << "\t- atomic files (in .cfg format such as printed by AtomHIC) where the name of the auxiliary property to use as descriptor for GMM has to be provided in the NameOfDescriptor argument, in this case it is recommanded to provide the properties of the descriptor in a \"InputDirectory/DescriptorProperties.ath\" file (an example file could be found in /data/ExampleFiles/), the descriptors can be filtered (by type or element for instance) by setting the value of DESCRIPTORS_FILTERING_TYPE in FixedParameters.ath file in the working directory (if not default value will be used)" << endl;
		cerr << "\t- or files listing only the descriptors" << endl;

		return EXIT_FAILURE;
	}
	string InputDir = argv[1];
	string DescriptorName, yaml_filename;
	unsigned int current_ind_read = 2;
	if( argc == 7 ){
		DescriptorName = argv[2];
		current_ind_read++;
	}
	unsigned int MaxDim;
        istringstream iss_md(argv[current_ind_read]);
	iss_md >> MaxDim;
	current_ind_read++;
	double MaxCor;
        istringstream iss_mc(argv[current_ind_read]);
	iss_mc >> MaxCor;
	current_ind_read++;
	double MinScore;
        istringstream iss_ms(argv[current_ind_read]);
	iss_ms >> MinScore;
	current_ind_read++;
	yaml_filename = argv[current_ind_read];

	Descriptors *MyDes;
	if( argc == 7 ) MyDes = new Descriptors(InputDir,DescriptorName);
	else MyDes = new Descriptors(InputDir);

	// check inputs
	if( MaxCor < 0. || MaxCor > 1. ){
        	cerr << "MaxCorrelation must be between 0 and 1." << endl;
		exit(EXIT_FAILURE);
	}
	if( MinScore < 0. || MinScore > 1. ){
        	cerr << "MinScore must be between 0 and 1." << endl;
		exit(EXIT_FAILURE);
	}
	
	// get descriptor properties
	Index nbDim = MyDes->getDim();
    	Index nbLab = MyDes->getNbLabels();
	unsigned int nbFilter = MyDes->getNbFilter();
	MyDes->constructSubarrays(true,true);

	vector<vector<unsigned int>> selected_dim;	
	// principal loop on the different filters (really need to be seriously tested after removing the corresponding nbody section in the yaml file!!)
	for(unsigned int f=0;f<nbFilter;f++){
		cout << "Treating " << MyDes->getFilterValue(f) << endl;
		// compute global average and class average
		Eigen::RowVectorXd globalMean = Eigen::RowVectorXd::Zero(nbDim);
    		vector<Eigen::RowVectorXd> classMeans;
    		classMeans.reserve(nbLab);
		MatrixXd *current_Des;
		unsigned int ntot = 0;
		for(unsigned int l=0;l<nbLab;l++){
			unsigned int n;
			current_Des = MyDes->getSubarray(n,MyDes->getFilterValue(f),MyDes->getLabels(l));
        		Eigen::RowVectorXd mean = Eigen::RowVectorXd::Zero(nbDim);
        		if( n > 0 ) mean = current_Des->colwise().mean();
			classMeans.push_back(mean);
			globalMean += n * mean;
			ntot += n;
		}
		globalMean /= static_cast<double>(ntot);
		// compute sum of squares between classes and total
		Eigen::RowVectorXd between = Eigen::RowVectorXd::Zero(nbDim);
		Eigen::RowVectorXd total = Eigen::RowVectorXd::Zero(nbDim);
		for(unsigned int l=0;l<nbLab;l++){
			unsigned int n;
			current_Des = MyDes->getSubarray(n,MyDes->getFilterValue(f),MyDes->getLabels(l));
			if( n == 0 ) continue;
        		const Eigen::RowVectorXd delta = classMeans[l] - globalMean;
        		between.array() += static_cast<double>(n) * delta.array().square();
        		const Eigen::MatrixXd centered = current_Des->rowwise() - globalMean;
			total.array() += centered.array().square().colwise().sum();
		}
		// compute discrimination score for each dimension
		Eigen::RowVectorXd scores = Eigen::RowVectorXd::Zero(nbDim);
		Eigen::RowVectorXd totalStd = Eigen::RowVectorXd::Zero(nbDim);
		double sum_var = 0.;
		for(unsigned int d=0;d<nbDim;d++){
			if (total(d) > 0.0){
				scores(d) = between(d) / total(d);
				totalStd(d) = sqrt(total(d));
				sum_var += totalStd(d);
			}
		}
		// sort candidates dimensions
    		vector<Index> candidates(nbDim);
		iota(candidates.begin(),candidates.end(),Index{0});
		sort(candidates.begin(),candidates.end(),[&scores](Index a, Index b){return scores(a) > scores(b);});
    		// Greedy selection with correlation filtering
    		vector<unsigned int> selected;
    		selected.reserve(min(MaxDim, static_cast<unsigned int>(nbDim)));
		unsigned int cur_nb_sel = 0;
		double var_expl = 0.;
		unsigned int nbRed = 0;
		for(unsigned int d=0;d<nbDim;d++){
			Index cur_ind = candidates[d];
        		// Ignore constant or insufficiently discriminative features.
        		if( total(cur_ind) <= 0.0 || scores(cur_ind) < MinScore ) continue;
        		bool redundant = false;
			for(unsigned int s=0;s<cur_nb_sel;s++){
				Index sel_ind = selected[s];
            			double covariance = 0.0;
            			// Compute covariance
				for(unsigned int l=0;l<nbLab;l++){
					unsigned int n;
					current_Des = MyDes->getSubarray(n,MyDes->getFilterValue(f),MyDes->getLabels(l));
					if( n == 0 ) continue;

            			    covariance += ( (current_Des->col(cur_ind).array() - globalMean(cur_ind)) * (current_Des->col(sel_ind).array() - globalMean(sel_ind)) ).sum();
            			}

            			double correlation = covariance / (totalStd(cur_ind) * totalStd(sel_ind));
            			if( abs(correlation) > MaxCor ){
            				redundant = true;
					nbRed++;
            				break;
            			}
			}
		        if( !redundant ){
				selected.push_back(static_cast<unsigned int>(cur_ind));
				cur_nb_sel++;
				var_expl += totalStd(cur_ind);
			}
			if( selected.size() >= MaxDim ) break;
		}
		cout << "The selected dimensions allows for explaining " << var_expl*100./sum_var << "\% of the variance" << endl;
		cout << nbRed << " dimensions have been excluded because of too high correlation (>" << MaxCor << ") with already selected dimensions" << endl;
		if( selected.size() != MaxDim ) cout << "WARNING the final number of dimension is lower than the maximum number of dimension provided, you may decrease the maximum correlation or the minimum score" << endl;
		string fname = "SelectedDim_"+MyDes->getFilterValue(f)+".dat";
	       	ofstream writefile(fname);
		for(unsigned int d=0;d<cur_nb_sel;d++) writefile << selected[d] << endl;
		writefile.close();
		selected_dim.push_back(selected);
	} // end loop on nbFilter
	
	cout << "Modifying yaml file accordingly" << endl;
	ifstream file(yaml_filename, ios::in);
	vector<string> to_read;
	string line;
	// first read and store all the yaml file
	if(file){
		while(getline(file,line)){
			to_read.push_back(line);
			if( !file ) break;
		}
		file.close();
	}else{
		cout << "The yaml file " << yaml_filename << " cannot be openned" << endl;
		exit(EXIT_FAILURE);
		return false;
	}
	vector<string> to_print = to_read;
	unsigned int nblines = to_read.size();
	vector<size_t> pos(nbFilter);
	vector<unsigned int> cur_dim(nbFilter,0);
	string buffer_s;
	vector<unsigned int> lines2rm;
	for(unsigned int l=0;l<nblines;l++){
		for(int f=0;f<nbFilter;f++){
			bool continue_read = true;
			pos[f] = to_read[l].find("- speciesblock: "+MyDes->getFilterValue(f));
			if(pos[f]!=string::npos){
				l += 13; // go to the nobdy section
				continue_read = true;
				while( continue_read ){
					istringstream iss(to_read[l]);
					iss >> buffer_s >> buffer_s;
					if( buffer_s == "{type:" ){
						bool tokeep = false;
						for(unsigned int d=0;d<selected_dim[f].size();d++){
							if( selected_dim[f][d] == cur_dim[f] ){
								tokeep = true;
								selected_dim[f].erase(selected_dim[f].begin()+d);
								break;
							}
						}
						if( !tokeep ) lines2rm.push_back(l);
						l++;
						cur_dim[f]++;
					}else{
						continue_read = false;
						break;
					}
					if( l == nblines ) break;
				}
			}
			if( l == nblines ) break;
			if( !continue_read ) f--;
		} // end loop on filters
		if( l == nblines ) break;
	}

	// remove the corresponding lines and print the new yaml file
	sort(lines2rm.begin(),lines2rm.end(),[](unsigned int a, unsigned int b){return (a>b);});
	for(unsigned int i=0;i<lines2rm.size();i++) to_read.erase(to_read.begin()+lines2rm[i]);
	ofstream writefile("Modified_"+yaml_filename);
	for(unsigned int i=0;i<to_read.size();i++) writefile << to_read[i] << endl;
	writefile.close();
	cout << "Modified_"+yaml_filename+" successfully writed" << endl;

	delete MyDes;
}
