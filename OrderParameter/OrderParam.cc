#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <stdio.h>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string.h>
#include <vector>
#include <array>
#include <random>
#include <complex>
#include <iomanip>

std::complex<double> spherical_harmonics(unsigned int l, int m, double theta, double phi){
	int mabs;
	if( m < 0 ) mabs = -m;
	else mabs = m;
	double leg = std::sph_legendre(l, mabs, theta);
	if( m < 0 ) leg *= pow(-1., .5*(m-mabs));
	std::complex<double> sph_harm(leg*cos(phi*m),leg*sin(phi*m));
	return sph_harm;
}

double gaussian(double x, double mu, double sigma);

using namespace std;

int main(){

        string fichier ="posfile";
	double rc = 5.; // cutting radius for neighbor search
	int lsph = 20; // spherical harmonic degree
	int nbPer = 4; // number of periodic images
	double tolSites = 5e-0;
	double zeronum = 1e-8;
	double PI = 3.141592654;
	const int bar_length = 30;
	double prog=0;
	unsigned int nbCellX(0), nbCellY(0), nbCellZ(0);
	double CellSizeX, CellSizeY, CellSizeZ;
	vector<vector<unsigned int>> Cells; // contains ions id belonging to the different cells

	char filename[15];
	char timestep_c[20];
	strcpy(filename, "Order_");
	char filename_distrib[30];
	strcpy(filename_distrib, "Disorder_");
	vector<double> Mg1NormFac, Mg2NormFac, SiNormFac, O1NormFac, O2NormFac, O3NormFac, MgNormFac, ONormFac, buffer_vector; // normalization factors for the different species and sites
	vector<unsigned int> buffer_vector_uint;
	vector<vector<double>> AtData, Nalpha; // contain the list of neighboring atom of all species (colatitude, longitude)
	vector<vector<unsigned int>> Malpha;  // contain the id number of neighboring atoms of the same species
	vector<vector<complex<double>>> Qalpha; // contain the qalpha_l_m values for each ion 
	int h_lo, h_hi, k_lo, k_hi, l_lo, l_hi;
	unsigned int nbAt, jmin, nbN, count, countSi, countMg, countO, countMg1, countMg2, countO1, countO2, countO3, timestep;
	double x, y, z, type, lx, ly, lz, d_squared, xp, yp, colat, longit, buffer_d, buffer_d_1, buffer_d_2, NormFac1_Mg1, NormFac1_Mg2, NormFac1_O1, NormFac1_O2, NormFac1_O3, NormFac2_Si, NormFac2_Mg1, NormFac2_Mg2, NormFac2_O1, NormFac2_O2, NormFac2_O3;
	bool buffer_bool;
	double rc_squared = pow(rc,2.);
	// get atomic data
	count = 0;
	lx = 0.;
	ly = 0.;
	lz = 0.;
        ifstream id1(fichier, ios::in);
        if(id1) {
		buffer_vector.push_back(0.);
		buffer_vector.push_back(0.);
		buffer_vector.push_back(0.);
		buffer_vector.push_back(0.);
                while(id1) {
			if( count == 0 ){
			       id1 >> timestep >> lx >> ly >> lz;
			       count = 1;
			}else{
				id1 >> type >> x >> y >> z;
                        	if(id1) {
					buffer_vector[0] = x;
					buffer_vector[1] = y;
					buffer_vector[2] = z;
					buffer_vector[3] = type;
                                	AtData.push_back(buffer_vector);
				}
                        }
                }
                id1.close();
        }
        else { cout << "Missing atomic position file (should be named posfile)" << endl; }
	if( lx == 0 || ly == 0 || lz == 0 ) cout << "Missing box lenghts in the position file (should be given at the first line of posfile" << endl;
	sprintf(timestep_c, "%d", timestep);
	strcat(filename, timestep_c);
	strcat(filename, ".xsf");
	strcat(filename_distrib, timestep_c);
	strcat(filename_distrib, ".txt");
	nbAt = AtData.size();

	// declare array which depend on number of ions
	double *Calpha = new double[nbAt];
	unsigned int *AtSites = new unsigned int[nbAt];
	double *Order = new double[nbAt];
	double *Disorder = new double[nbAt];
	int AtOverTen=(nbAt/200);
	// wrap ions
	for(unsigned int i=0;i<nbAt;i++){
		if( ( AtData[i][0] < 0 ) || ( AtData[i][0] > lx ) ){
			buffer_d = lx*((AtData[i][0]/lx)-floor(AtData[i][0]/lx));
			AtData[i][0] = buffer_d;
		}
		if( ( AtData[i][1] < 0 ) || ( AtData[i][1] > ly ) ){
			buffer_d = ly*((AtData[i][1]/ly)-floor(AtData[i][1]/ly));
			AtData[i][1] = buffer_d;
		}
		// Does not wrap in z direction as it can be free surface //TODO : find solution for computing gaussian at the end, as we need to wrap for neighbor search
		if( ( AtData[i][2] < 0 ) || ( AtData[i][2] >= lz ) ){
			buffer_d = lz*((AtData[i][2]/lz)-floor(AtData[i][2]/lz));
			AtData[i][2] = buffer_d;
		}
	}

	// search neighbors and store their colatitude and azimuth
	// decide using brute force or cell list algorithm for neighboring research
	double AtDensity_Fo;
	AtDensity_Fo = 28./(4.8435*10.1864*6.001);
	unsigned int NbAt_CellList; // number of atom above which we use cell list algo for neigh search (correspond to 4 times 27 times a cell size)
	NbAt_CellList = rc*rc*rc*27*2*AtDensity_Fo;
	unsigned int BrutForce = 0;
	unsigned int MulForCell = 1;
	if( lx < rc ) MulForCell *= ceil(rc/lx);
	if( ly < rc ) MulForCell *= ceil(rc/ly);
	if( lz < rc ) MulForCell *= ceil(rc/lz);
	cout << MulForCell << endl;
	//if( nbAt > NbAt_CellList ) BrutForce = 0;
	//else BrutForce = 1;
	cout << NbAt_CellList << endl;
	cout << "searching neigbors" << endl;
	// Brut force algorithm for small systems
	double start_s=clock();
	if( BrutForce == 1 ){
	cout << "we use brute force" << endl;
	cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
	for(unsigned int i=0;i<nbAt;i++){
		if(i%AtOverTen == 0){
			prog = double(i)/double(nbAt);
			cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
		}
		h_lo = 0;
		h_hi = 0;
		k_lo = 0;
		k_hi = 0;
		l_lo = 0;
		l_hi = 0;
		Nalpha.push_back(vector<double>());
		Malpha.push_back(vector<unsigned int>());
		for(unsigned int k=0;k<nbPer+1;k++){
                	if( (AtData[i][0]-rc) < -(k*lx) ) h_lo = -k-1;
                	if( (AtData[i][0]+rc) > lx*(1+k) ) h_hi = k+1;
                	if( (AtData[i][1]-rc) < -(k*ly) ) k_lo = -k-1;
                	if( (AtData[i][1]+rc) > ly*(1+k) ) k_hi = k+1;
                	if( (AtData[i][2]-rc) < -(k*lz) ) l_lo = -k-1;
                	if( (AtData[i][2]+rc) > lz*(1+k) ) l_hi = k+1;
		}
		for(unsigned int j=0;j<nbAt;j++){
			for(int xcl=h_lo;xcl<h_hi+1;xcl++){
				for(int ycl=k_lo;ycl<k_hi+1;ycl++){
					for(int zcl=l_lo;zcl<l_hi+1;zcl++){
						d_squared = pow(AtData[j][0]-AtData[i][0]+(xcl*lx),2) + pow(AtData[j][1]-AtData[i][1]+(ycl*ly),2) + pow(AtData[j][2]-AtData[i][2]+(zcl*lz),2);
						if( d_squared > zeronum && d_squared < rc_squared ){
							colat = acos((AtData[j][2]-AtData[i][2]+(zcl*lz))/pow(d_squared,.5));
							xp = AtData[j][0]-AtData[i][0]+(xcl*lx);
							yp = AtData[j][1]-AtData[i][1]+(ycl*ly);
				                        if( xp > 0 ) longit = atan(yp/xp);
	                    				else if( ( xp < 0 ) && ( yp >= 0 ) ) longit = atan(yp/xp) + PI;
	                    				else if( ( xp < 0 ) and ( yp < 0 ) ) longit = atan(yp/xp) - PI;
	                    				else if( ( fabs(xp) < zeronum ) and ( yp > 0 ) ) longit = PI/2.;
	                    				else if( ( fabs(xp) < zeronum ) and ( yp < 0 ) ) longit = -PI/2.;
	                    				else if( ( fabs(xp) < zeronum ) and ( fabs(yp) < zeronum ) ) longit = 0.;
							Nalpha[i].push_back(colat);
							Nalpha[i].push_back(longit);
							if( AtData[i][3] == AtData[j][3] ) Malpha[i].push_back(j);
						}
					}
				}
			}
		}
	}
	}else{
	cout << "we use cell list" << endl;
	// Cell list algorithm for wider systems
	// construct cells
	int NeighCellX, NeighCellY, NeighCellZ;
	nbCellX = floor(lx/rc);
	if( nbCellX == 0 ){
		nbCellX = 1;
		CellSizeX = lx;
		NeighCellX = ceil(rc/lx);
	}else{
		CellSizeX = lx/nbCellX;
		NeighCellX = 1;
	}
	nbCellY = floor(ly/rc);
	if( nbCellY == 0 ){
		nbCellY = 1;
		CellSizeY = ly;
		NeighCellY = ceil(rc/ly);
	}else{
		CellSizeY = ly/nbCellY;
		NeighCellY = 1;
	}
	nbCellZ = floor(lz/rc);
	if( nbCellZ == 0 ){
		nbCellZ = 1;
		CellSizeZ = lz;
		NeighCellZ = ceil(rc/lz);
	}else{
		CellSizeZ = lz/nbCellZ;
		NeighCellZ = 1;
	}
	for(unsigned int i=0;i<nbCellX;i++){
		for(unsigned int j=0;j<nbCellY;j++){
			for(unsigned int k=0;k<nbCellZ;k++){
				Cells.push_back(vector<unsigned int>());
				if( (i == nbCellX-1) && (j == nbCellY-1) && (k == nbCellZ-1) ){
					for(unsigned at=0;at<nbAt;at++){
						if( (AtData[at][0]>=i*CellSizeX) && (AtData[at][0]<=(i+1)*CellSizeX) && (AtData[at][1]>=j*CellSizeY) && (AtData[at][1]<=(j+1)*CellSizeY) && (AtData[at][2]>=k*CellSizeZ) && (AtData[at][2]<=(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
					}
				}else if( (i == nbCellX-1) && (j == nbCellY-1) ){  
					for(unsigned at=0;at<nbAt;at++){
						if( (AtData[at][0]>=i*CellSizeX) && (AtData[at][0]<=(i+1)*CellSizeX) && (AtData[at][1]>=j*CellSizeY) && (AtData[at][1]<=(j+1)*CellSizeY) && (AtData[at][2]>=k*CellSizeZ) && (AtData[at][2]<(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
					}
				}else if( (j == nbCellY-1) && (k == nbCellZ-1) ){  
					for(unsigned at=0;at<nbAt;at++){
						if( (AtData[at][0]>=i*CellSizeX) && (AtData[at][0]<(i+1)*CellSizeX) && (AtData[at][1]>=j*CellSizeY) && (AtData[at][1]<=(j+1)*CellSizeY) && (AtData[at][2]>=k*CellSizeZ) && (AtData[at][2]<=(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
					}
				}else if( (i == nbCellX-1) && (k == nbCellZ-1) ){  
					for(unsigned at=0;at<nbAt;at++){
						if( (AtData[at][0]>=i*CellSizeX) && (AtData[at][0]<=(i+1)*CellSizeX) && (AtData[at][1]>=j*CellSizeY) && (AtData[at][1]<(j+1)*CellSizeY) && (AtData[at][2]>=k*CellSizeZ) && (AtData[at][2]<=(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
					}
				}else if(i == nbCellX-1){  
					for(unsigned at=0;at<nbAt;at++){
						if( (AtData[at][0]>=i*CellSizeX) && (AtData[at][0]<=(i+1)*CellSizeX) && (AtData[at][1]>=j*CellSizeY) && (AtData[at][1]<(j+1)*CellSizeY) && (AtData[at][2]>=k*CellSizeZ) && (AtData[at][2]<(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
					}
				}else if(j == nbCellY-1){  
					for(unsigned at=0;at<nbAt;at++){
						if( (AtData[at][0]>=i*CellSizeX) && (AtData[at][0]<(i+1)*CellSizeX) && (AtData[at][1]>=j*CellSizeY) && (AtData[at][1]<=(j+1)*CellSizeY) && (AtData[at][2]>=k*CellSizeZ) && (AtData[at][2]<(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
					}
				}else if(k == nbCellZ-1){  
					for(unsigned at=0;at<nbAt;at++){
						if( (AtData[at][0]>=i*CellSizeX) && (AtData[at][0]<(i+1)*CellSizeX) && (AtData[at][1]>=j*CellSizeY) && (AtData[at][1]<(j+1)*CellSizeY) && (AtData[at][2]>=k*CellSizeZ) && (AtData[at][2]<=(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
					}
				}else{
					for(unsigned at=0;at<nbAt;at++){
						if( (AtData[at][0]>=i*CellSizeX) && (AtData[at][0]<(i+1)*CellSizeX) && (AtData[at][1]>=j*CellSizeY) && (AtData[at][1]<(j+1)*CellSizeY) && (AtData[at][2]>=k*CellSizeZ) && (AtData[at][2]<(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
					}
				}
			}
		}
	}
	cout << endl;
	// Initialize vectors containing neighboring infos
	for(unsigned int i=0;i<nbAt;i++){
		Nalpha.push_back(vector<double>());
		Malpha.push_back(vector<unsigned int>());
	}
	// Perform neighboring research
	double xpos,ypos,zpos;
	int ibx, jby, kbz;
	int Nclx, Ncly, Nclz;
	cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
	for(unsigned int i=0;i<nbCellX;i++){
		for(unsigned int j=0;j<nbCellY;j++){
			for(unsigned int k=0;k<nbCellZ;k++){
				prog = double(i*nbCellZ*nbCellY+j*nbCellZ+k)/double(nbCellX*nbCellY*nbCellZ);
				cout << "\r[" << string(floor(bar_length*prog),'X') << string(ceil(bar_length*(1-prog)),'-') << "] " << setprecision(3) << 100*prog << "%";
				for(unsigned int at1 = 0; at1<Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].size(); at1++){
					xpos = AtData[Cells[i*nbCellZ*nbCellY+j*nbCellZ+k][at1]][0];
					ypos = AtData[Cells[i*nbCellZ*nbCellY+j*nbCellZ+k][at1]][1];
					zpos = AtData[Cells[i*nbCellZ*nbCellY+j*nbCellZ+k][at1]][2];
					for(int bx=-NeighCellX;bx<NeighCellX+1;bx++){
						for(int by=-NeighCellY;by<NeighCellY+1;by++){
							for(int bz=-NeighCellZ;bz<NeighCellZ+1;bz++){
								// Search using periodic BC which cell use in case of border cell and what CL applied to ion pos
								// for a non border cell =>
								Nclx = 0;
								Ncly = 0;
								Nclz = 0;
								ibx = bx+i;
								jby = by+j;
								kbz = bz+k;
								// border cells
								if( i == nbCellX-1 && bx == 1 ){
									ibx = 0;
									Nclx = 1;
								}else if( i == 0 && bx == -1 ){
									ibx = nbCellX-1;
									Nclx = -1;
								}
								if( j == nbCellY-1 && by == 1 ){
									jby = 0;
									Ncly = 1;
								}else if( j == 0 && by == -1 ){
									jby = nbCellY-1;
									Ncly = -1;
								}
								if( k == nbCellZ-1 && bz == 1 ){
									kbz = 0;
									Nclz = 1;
								}else if( k == 0 && bz == -1 ){
									kbz = nbCellZ-1;
									Nclz = -1;
								}
								for(unsigned int at2=0;at2<Cells[ibx*nbCellZ*nbCellY+jby*nbCellZ+kbz].size();at2++){
									d_squared = pow(AtData[Cells[ibx*nbCellY*nbCellZ+jby*nbCellZ+kbz][at2]][0]-xpos+Nclx*lx,2.)+pow(AtData[Cells[ibx*nbCellY*nbCellZ+jby*nbCellZ+kbz][at2]][1]-ypos+Ncly*ly,2.)+pow(AtData[Cells[ibx*nbCellZ*nbCellY+jby*nbCellZ+kbz][at2]][2]-zpos+Nclz*lz,2.);
									if( d_squared > zeronum && d_squared < rc_squared ){
										colat = acos((AtData[Cells[ibx*nbCellZ*nbCellY+jby*nbCellZ+kbz][at2]][2]-zpos+Nclz*lz)/pow(d_squared,.5));
										xp = AtData[Cells[ibx*nbCellZ*nbCellY+jby*nbCellZ+kbz][at2]][0]-xpos+Nclx*lx;
										yp = AtData[Cells[ibx*nbCellZ*nbCellY+jby*nbCellZ+kbz][at2]][1]-ypos+Ncly*ly;
				                			        if( xp > 0 ) longit = atan(yp/xp);
	                    							else if( ( xp < 0 ) && ( yp >= 0 ) ) longit = atan(yp/xp) + M_PI;
	                    							else if( ( xp < 0 ) and ( yp < 0 ) ) longit = atan(yp/xp) - M_PI;
	                    							else if( ( fabs(xp) < zeronum ) and ( yp > 0 ) ) longit = M_PI/2.;
	                    							else if( ( fabs(xp) < zeronum ) and ( yp < 0 ) ) longit = -M_PI/2.;
	                    							else if( ( fabs(xp) < zeronum ) and ( fabs(yp) < zeronum ) ) longit = 0.;
										Nalpha[Cells[i*nbCellZ*nbCellY+j*nbCellZ+k][at1]].push_back(colat);
										Nalpha[Cells[i*nbCellZ*nbCellY+j*nbCellZ+k][at1]].push_back(longit);
										if( AtData[Cells[i*nbCellZ*nbCellY+j*nbCellZ+k][at1]][3] == AtData[Cells[ibx*nbCellY*nbCellZ+jby*nbCellZ+kbz][at2]][3] ) Malpha[Cells[i*nbCellZ*nbCellY+j*nbCellZ+k][at1]].push_back(Cells[ibx*nbCellZ*nbCellY+jby*nbCellZ+kbz][at2]);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	cout << endl;
	} // end neighbor search
	cout << endl;
	double stop_s=clock();
	cout << "time = " << stop_s-start_s << endl;
	cout << "compute spherical harmonics" << endl;
	for(unsigned int i=0;i<nbAt;i++){
		if(i%AtOverTen == 0){
			prog = double(i)/double(nbAt);
			cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
		}
		Qalpha.push_back(vector<complex<double>>());
		Calpha[i] = 0;
		for(int l=-lsph;l<lsph+1;l++){
			Qalpha[i].push_back((0.,0.));
			for(unsigned int N=0;N<Nalpha[i].size()/2;N++){
				Qalpha[i][l+lsph] += spherical_harmonics(lsph,l,Nalpha[i][N*2],Nalpha[i][N*2+1]);
			}
			Calpha[i] += (pow(Qalpha[i][l+lsph].real(), 2.) + pow(Qalpha[i][l+lsph].imag(), 2.));
		}
		buffer_d = Calpha[i];
	}
	cout << endl;
	
	countMg = 0;
	countO = 0;
	for(unsigned int i=0;i<nbAt;i++){
		if(AtData[i][3]==1){
			if(countMg==0){
				MgNormFac.push_back(Calpha[i]);
				MgNormFac.push_back(1);
				countMg += 1;
			}else{
				buffer_bool = false;
				for(unsigned int j=0;j<MgNormFac.size()/2;j++){
					if( fabs(Calpha[i]-MgNormFac[j*2]) < tolSites ){
						MgNormFac[j*2] *= MgNormFac[j*2+1];
						MgNormFac[j*2] += Calpha[i];
						MgNormFac[j*2+1] += 1;
						MgNormFac[j*2] /= MgNormFac[j*2+1];
						buffer_bool = true;
						break;
					}
				}
				if(!buffer_bool){
					MgNormFac.push_back(Calpha[i]);
					MgNormFac.push_back(1);
				}
			}
		}else if(AtData[i][3] == 3){
			if(countO==0){
				ONormFac.push_back(Calpha[i]);
				ONormFac.push_back(1);
				countO += 1;
			}else{
				buffer_bool = false;
				for(unsigned int j=0;j<ONormFac.size()/2;j++){
					if( fabs(Calpha[i]-ONormFac[j*2]) < tolSites ){
						ONormFac[j*2] *= ONormFac[j*2+1];
						ONormFac[j*2] += Calpha[i];
						ONormFac[j*2+1] += 1;
						ONormFac[j*2] /= ONormFac[j*2+1];
						buffer_bool = true;
						break;
					}
				}
				if(!buffer_bool){
					ONormFac.push_back(Calpha[i]);
					ONormFac.push_back(1);
				}
			}
		}
	}

	count = 0;
	for(unsigned int n=0;n<2;n++){
		buffer_d = MgNormFac[1];
		buffer_d_1 = MgNormFac[0];
		buffer_d_2 = 0;
		for(unsigned int i=1;i<MgNormFac.size()/2;i++){
			if(MgNormFac[i*2+1] > buffer_d){
				buffer_d = MgNormFac[i*2+1];
				buffer_d_1 = MgNormFac[i*2];
				buffer_d_2 = i*2;
			}
		}
		if(count==0) NormFac1_Mg1 = buffer_d_1;
		else if(count==1) NormFac1_Mg2 = buffer_d_1;
		MgNormFac.erase(MgNormFac.begin()+buffer_d_2+1);
		MgNormFac.erase(MgNormFac.begin()+buffer_d_2);
		count += 1;
	}

	count = 0;
	for(unsigned int n=0;n<3;n++){
		buffer_d = ONormFac[1];
		buffer_d_1 = ONormFac[0];
		buffer_d_2 = 0;
		for(unsigned int i=1;i<ONormFac.size()/2;i++){
			if(ONormFac[i*2+1] > buffer_d){
				buffer_d = ONormFac[i*2+1];
				buffer_d_1 = ONormFac[i*2];
				buffer_d_2 = i*2;
			}
		}
		if(count==0) NormFac1_O1 = buffer_d_1;
		else if(count==1) NormFac1_O2 = buffer_d_1;
		else if(count==2) NormFac1_O3 = buffer_d_1;
		ONormFac.erase(ONormFac.begin()+buffer_d_2+1);
		ONormFac.erase(ONormFac.begin()+buffer_d_2);
		count += 1;
	}

	// assign new id to all corresponding to the site => 1 Si, 2 Mg1, 3 Mg2, 4 O1, 5 O2, 6 O3
	for(unsigned int i=0;i<nbAt;i++){
		if( AtData[i][3] == 2 ) AtSites[i] = 1;
		else if( AtData[i][3] == 1 ){
			buffer_d = fabs(Calpha[i]-NormFac1_Mg1);
			if( fabs(Calpha[i]-NormFac1_Mg2) < buffer_d ) AtSites[i] = 3;
			else AtSites[i] = 2;
		}else if( AtData[i][3] == 3 ){
			buffer_d = fabs(Calpha[i]-NormFac1_O1);
			AtSites[i] = 4;
			if( fabs(Calpha[i]-NormFac1_O2) < buffer_d ){
				buffer_d = fabs(Calpha[i]-NormFac1_O2);
				AtSites[i] = 5;
			}
			if( fabs(Calpha[i]-NormFac1_O3) < buffer_d ){
				AtSites[i] = 6;
			}
		}
	}

	//update neighbor list for corresponding to ions of the same site
	for(unsigned int i=0;i<nbAt;i++){
		buffer_vector_uint.clear();
		for(unsigned int j=0;j<Malpha[i].size();j++){
			if( AtSites[i] == AtSites[Malpha[i][j]] ) buffer_vector_uint.push_back(Malpha[i][j]);
		}
		Malpha[i].clear();
		for(unsigned int j=0;j<buffer_vector_uint.size();j++) Malpha[i].push_back(buffer_vector_uint[j]);
	}

	// Compute order parameter using the classic formulation
	for(unsigned int i=0;i<nbAt;i++){
		Order[i] = 0;
		nbN = Malpha[i].size();
		for(unsigned int j=0;j<nbN;j++){
			for(unsigned int l=0;l<(lsph*2+1);l++){
				Order[i] += ((Qalpha[i][l].real()*Qalpha[Malpha[i][j]][l].real())+Qalpha[i][l].imag()*Qalpha[Malpha[i][j]][l].imag())/(pow(Calpha[i],.5)*pow(Calpha[Malpha[i][j]],.5)); 
			}
		}
		if( nbN == 0 ) Order[i] = 0;
		else Order[i] /= nbN;
	}

	// renormalize according to the different sites
	// FIRST TEST : consider that sites established before are the right ones
	countSi = 0;
	countMg1 = 0;
	countMg2 = 0;
	countO1 = 0;
	countO2 = 0;
	countO3 = 0;
	tolSites = 1e-2;
	for(unsigned int i=0;i<nbAt;i++){
		if(AtSites[i]==1){
			if(countSi==0){
				SiNormFac.push_back(Order[i]);
				SiNormFac.push_back(1);
				countSi += 1;
			}else{
				buffer_bool = false;
				for(unsigned int j=0;j<SiNormFac.size()/2;j++){
					if( fabs(Order[i]-SiNormFac[j*2]) < tolSites ){
						SiNormFac[j*2] *= SiNormFac[j*2+1];
						SiNormFac[j*2] += Order[i];
						SiNormFac[j*2+1] += 1;
						SiNormFac[j*2] /= SiNormFac[j*2+1];
						buffer_bool = true;
						break;
					}
				}
				if(!buffer_bool){
					SiNormFac.push_back(Order[i]);
					SiNormFac.push_back(1);
				}
			}
		}else if(AtSites[i]==2){
			if(countMg1 ==0){
				Mg1NormFac.push_back(Order[i]);
				Mg1NormFac.push_back(1);
				countMg1 += 1;
			}else{
				buffer_bool = false;
				for(unsigned int j=0;j<Mg1NormFac.size()/2;j++){
					if( fabs(Order[i]-Mg1NormFac[j*2]) < tolSites ){
						Mg1NormFac[j*2] *= Mg1NormFac[j*2+1];
						Mg1NormFac[j*2] += Order[i];
						Mg1NormFac[j*2+1] += 1;
						Mg1NormFac[j*2] /= Mg1NormFac[j*2+1];
						buffer_bool = true;
						break;
					}
				}
				if(!buffer_bool){
					Mg1NormFac.push_back(Order[i]);
					Mg1NormFac.push_back(1);
				}
			}
		}else if(AtSites[i]==3){
			if(countMg2 ==0){
				Mg2NormFac.push_back(Order[i]);
				Mg2NormFac.push_back(1);
				countMg2 += 1;
			}else{
				buffer_bool = false;
				for(unsigned int j=0;j<Mg2NormFac.size()/2;j++){
					if( fabs(Order[i]-Mg2NormFac[j*2]) < tolSites ){
						Mg2NormFac[j*2] *= Mg2NormFac[j*2+1];
						Mg2NormFac[j*2] += Order[i];
						Mg2NormFac[j*2+1] += 1;
						Mg2NormFac[j*2] /= Mg2NormFac[j*2+1];
						buffer_bool = true;
						break;
					}
				}
				if(!buffer_bool){
					Mg2NormFac.push_back(Order[i]);
					Mg2NormFac.push_back(1);
				}
			}
		}else if(AtSites[i]==4){
			if(countO1 ==0){
				O1NormFac.push_back(Order[i]);
				O1NormFac.push_back(1);
				countO1 += 1;
			}else{
				buffer_bool = false;
				for(unsigned int j=0;j<O1NormFac.size()/2;j++){
					if( fabs(Order[i]-O1NormFac[j*2]) < tolSites ){
						O1NormFac[j*2] *= O1NormFac[j*2+1];
						O1NormFac[j*2] += Order[i];
						O1NormFac[j*2+1] += 1;
						O1NormFac[j*2] /= O1NormFac[j*2+1];
						buffer_bool = true;
						break;
					}
				}
				if(!buffer_bool){
					O1NormFac.push_back(Order[i]);
					O1NormFac.push_back(1);
				}
			}
		}else if(AtSites[i]==5){
			if(countO2 ==0){
				O2NormFac.push_back(Order[i]);
				O2NormFac.push_back(1);
				countO2 += 1;
			}else{
				buffer_bool = false;
				for(unsigned int j=0;j<O2NormFac.size()/2;j++){
					if( fabs(Order[i]-O2NormFac[j*2]) < tolSites ){
						O2NormFac[j*2] *= O2NormFac[j*2+1];
						O2NormFac[j*2] += Order[i];
						O2NormFac[j*2+1] += 1;
						O2NormFac[j*2] /= O2NormFac[j*2+1];
						buffer_bool = true;
						break;
					}
				}
				if(!buffer_bool){
					O2NormFac.push_back(Order[i]);
					O2NormFac.push_back(1);
				}
			}
		}else if(AtSites[i]==6){
			if(countO3 ==0){
				O3NormFac.push_back(Order[i]);
				O3NormFac.push_back(1);
				countO3 += 1;
			}else{
				buffer_bool = false;
				for(unsigned int j=0;j<O3NormFac.size()/2;j++){
					if( fabs(Order[i]-O3NormFac[j*2]) < tolSites ){
						O3NormFac[j*2] *= O3NormFac[j*2+1];
						O3NormFac[j*2] += Order[i];
						O3NormFac[j*2+1] += 1;
						O3NormFac[j*2] /= O3NormFac[j*2+1];
						buffer_bool = true;
						break;
					}
				}
				if(!buffer_bool){
					O3NormFac.push_back(Order[i]);
					O3NormFac.push_back(1);
				}
			}
		}
	}

	// find the most represented norm factors
	NormFac2_Si = SiNormFac[0];
	buffer_d = SiNormFac[1];
	for(unsigned int i=1;i<SiNormFac.size()/2;i++){
		if( SiNormFac[i*2+1] >  buffer_d ){
			buffer_d = SiNormFac[i*2+1];
			NormFac2_Si = SiNormFac[i*2];
		}
	}
	NormFac2_Mg1 = Mg1NormFac[0];
	buffer_d = Mg1NormFac[1];
	for(unsigned int i=1;i<Mg1NormFac.size()/2;i++){
		if( Mg1NormFac[i*2+1] >  buffer_d ){
			buffer_d = Mg1NormFac[i*2+1];
			NormFac2_Mg1 = Mg1NormFac[i*2];
		}
	}
	NormFac2_Mg2 = Mg2NormFac[0];
	buffer_d = Mg2NormFac[1];
	for(unsigned int i=1;i<Mg2NormFac.size()/2;i++){
		if( Mg2NormFac[i*2+1] >  buffer_d ){
			buffer_d = Mg2NormFac[i*2+1];
			NormFac2_Mg2 = Mg2NormFac[i*2];
		}
	}
	NormFac2_O1 = O1NormFac[0];
	buffer_d = O1NormFac[1];
	for(unsigned int i=1;i<O1NormFac.size()/2;i++){
		if( O1NormFac[i*2+1] >  buffer_d ){
			buffer_d = O1NormFac[i*2+1];
			NormFac2_O1 = O1NormFac[i*2];
		}
	}
	NormFac2_O2 = O2NormFac[0];
	buffer_d = O2NormFac[1];
	for(unsigned int i=1;i<O2NormFac.size()/2;i++){
		if( O2NormFac[i*2+1] >  buffer_d ){
			buffer_d = O2NormFac[i*2+1];
			NormFac2_O2 = O2NormFac[i*2];
		}
	}
	NormFac2_O3 = O3NormFac[0];
	buffer_d = O3NormFac[1];
	for(unsigned int i=1;i<O3NormFac.size()/2;i++){
		if( O3NormFac[i*2+1] >  buffer_d ){
			buffer_d = O3NormFac[i*2+1];
			NormFac2_O3 = O3NormFac[i*2];
		}
	}
	// and renormalize order parameters
	for(unsigned int i=0;i<nbAt;i++){
		if(AtSites[i]==1) Order[i] /= NormFac2_Si;
		else if(AtSites[i]==2) Order[i] /= NormFac2_Mg1;
		else if(AtSites[i]==3) Order[i] /= NormFac2_Mg2;
		else if(AtSites[i]==4) Order[i] /= NormFac2_O1;
		else if(AtSites[i]==5) Order[i] /= NormFac2_O2;
		else if(AtSites[i]==6) Order[i] /= NormFac2_O3;
	}


	// END FIRST TEST
	// SECOND TEST : reassign site by identifying the most represented order parameters
	//countSi = 0;
	//countMg = 0;
	//countO = 0;
	//tolSites = 1e-1;
	//MgNormFac.clear();
	//ONormFac.clear();
	//for(unsigned int i=0;i<nbAt;i++){
	//	if(AtData[i][3]==2){
	//		if(countSi==0){
	//			SiNormFac.push_back(Order[i]);
	//			SiNormFac.push_back(1);
	//			countSi += 1;
	//		}else{
	//			buffer_bool = false;
	//			for(unsigned int j=0;j<SiNormFac.size()/2;j++){
	//				if( fabs(Order[i]-SiNormFac[j*2]) < tolSites ){
	//					SiNormFac[j*2] *= SiNormFac[j*2+1];
	//					SiNormFac[j*2] += Order[i];
	//					SiNormFac[j*2+1] += 1;
	//					SiNormFac[j*2] /= SiNormFac[j*2+1];
	//					buffer_bool = true;
	//					break;
	//				}
	//			}
	//			if(!buffer_bool){
	//				SiNormFac.push_back(Order[i]);
	//				SiNormFac.push_back(1);
	//			}
	//		}
	//	}else if(AtData[i][3]==1){
	//		if(countMg ==0){
	//			MgNormFac.push_back(Order[i]);
	//			MgNormFac.push_back(1);
	//			countMg += 1;
	//		}else{
	//			buffer_bool = false;
	//			for(unsigned int j=0;j<MgNormFac.size()/2;j++){
	//				if( fabs(Order[i]-MgNormFac[j*2]) < tolSites ){
	//					MgNormFac[j*2] *= MgNormFac[j*2+1];
	//					MgNormFac[j*2] += Order[i];
	//					MgNormFac[j*2+1] += 1;
	//					MgNormFac[j*2] /= MgNormFac[j*2+1];
	//					buffer_bool = true;
	//					break;
	//				}
	//			}
	//			if(!buffer_bool){
	//				MgNormFac.push_back(Order[i]);
	//				MgNormFac.push_back(1);
	//			}
	//		}
	//	}else if(AtData[i][3]==3){
	//		if(countO ==0){
	//			ONormFac.push_back(Order[i]);
	//			ONormFac.push_back(1);
	//			countO += 1;
	//		}else{
	//			buffer_bool = false;
	//			for(unsigned int j=0;j<ONormFac.size()/2;j++){
	//				if( fabs(Order[i]-ONormFac[j*2]) < tolSites ){
	//					ONormFac[j*2] *= ONormFac[j*2+1];
	//					ONormFac[j*2] += Order[i];
	//					ONormFac[j*2+1] += 1;
	//					ONormFac[j*2] /= ONormFac[j*2+1];
	//					buffer_bool = true;
	//					break;
	//				}
	//			}
	//			if(!buffer_bool){
	//				ONormFac.push_back(Order[i]);
	//				ONormFac.push_back(1);
	//			}
	//		}
	//	}
	//}

	//NormFac2_Si = SiNormFac[0];
	//buffer_d = SiNormFac[1];
	//for(unsigned int i=1;i<SiNormFac.size()/2;i++){
	//	if( SiNormFac[i*2+1] >  buffer_d ){
	//		buffer_d = SiNormFac[i*2+1];
	//		NormFac2_Si = SiNormFac[i*2];
	//	}
	//}

	//count = 0;
	//for(unsigned int n=0;n<2;n++){
	//	buffer_d = MgNormFac[1];
	//	buffer_d_1 = MgNormFac[0];
	//	buffer_d_2 = 0;
	//	for(unsigned int i=1;i<MgNormFac.size()/2;i++){
	//		if(MgNormFac[i*2+1] > buffer_d){
	//			buffer_d = MgNormFac[i*2+1];
	//			buffer_d_1 = MgNormFac[i*2];
	//			buffer_d_2 = i*2;
	//		}
	//	}
	//	if(count==0) NormFac2_Mg1 = buffer_d_1;
	//	else if(count==1) NormFac2_Mg2 = buffer_d_1;
	//	MgNormFac.erase(MgNormFac.begin()+buffer_d_2+1);
	//	MgNormFac.erase(MgNormFac.begin()+buffer_d_2);
	//	count += 1;
	//}

	//count = 0;
	//for(unsigned int n=0;n<3;n++){
	//	buffer_d = ONormFac[1];
	//	buffer_d_1 = ONormFac[0];
	//	buffer_d_2 = 0;
	//	for(unsigned int i=1;i<ONormFac.size()/2;i++){
	//		if(ONormFac[i*2+1] > buffer_d){
	//			buffer_d = ONormFac[i*2+1];
	//			buffer_d_1 = ONormFac[i*2];
	//			buffer_d_2 = i*2;
	//		}
	//	}
	//	if(count==0) NormFac2_O1 = buffer_d_1;
	//	else if(count==1) NormFac2_O2 = buffer_d_1;
	//	else if(count==2) NormFac2_O3 = buffer_d_1;
	//	ONormFac.erase(ONormFac.begin()+buffer_d_2+1);
	//	ONormFac.erase(ONormFac.begin()+buffer_d_2);
	//	count += 1;
	//}
	//// renormalize by searching the closest normalization factor for the different species
	//double tolRenorm = 1.05;
	//for(unsigned int i=0;i<nbAt;i++){
	//	if(AtData[i][3] == 1){
	//		buffer_d = fabs(Order[i]-NormFac2_Mg1);
	//		buffer_d_1 = NormFac2_Mg1;
	//		if( ( (Order[i] < NormFac2_Mg2*tolRenorm) && ( fabs(Order[i]-NormFac2_Mg2) < buffer_d ) ) || (Order[i] > NormFac2_Mg1*tolRenorm) ) buffer_d_1 = NormFac2_Mg2;
	//		Order[i] /= buffer_d_1;
	//	}else if(AtData[i][3] == 2){
	//		Order[i] /= NormFac2_Si;
	//	}else if(AtData[i][3] == 3){
	//		buffer_d = fabs(Order[i]-NormFac2_O1);
	//		buffer_d_1 = NormFac2_O1;
	//		buffer_d_2 = (NormFac2_O1*tolRenorm)-Order[i];
	//		if( ( (Order[i] < NormFac2_O2*tolRenorm) && ( fabs(Order[i]-NormFac2_O2) < buffer_d ) ) || (0. > buffer_d_2) ){
	//			buffer_d = fabs(Order[i]-NormFac2_O2);
	//		        buffer_d_1 = NormFac2_Mg2;
	//			buffer_d_2 = (NormFac2_O2*tolRenorm)-Order[i];
	//		}
	//		if( ( (Order[i] < NormFac2_O3*tolRenorm) && ( fabs(Order[i]-NormFac2_O3) < buffer_d ) ) || (0. > buffer_d_2) ){
	//			buffer_d = fabs(Order[i]-NormFac2_O3);
	//		        buffer_d_1 = NormFac2_O3;
	//		}
	//		Order[i] /= buffer_d_1;
	//	}
	//}
	// END SECOND TEST

	for(unsigned int i=0;i<nbAt;i++){
	       if( Order[i] > 1 ) Order[i] = 1;
       	       Disorder[i] = 1-fabs(Order[i]);
	}
	
	// Compute z distribution of disorder assimiling each ion to gaussian
	// get min and max z values 
	double MinZ = AtData[0][2];
	double MaxZ = AtData[0][2];
	for(unsigned int i=1;i<nbAt;i++){
		if(AtData[i][2] < MinZ) MinZ = AtData[i][2];
		if(AtData[i][2] > MaxZ) MaxZ = AtData[i][2];
	}
	double sigma = 1.; // A	
	unsigned int nbPts = 100;
	double *z_pos = new double[nbPts];
	double *diso_pos = new double[nbPts];
	double MinDiso;
	double ZMaxDiso;
	double MaxDiso;
	double NormFacGauss;
	for(unsigned int i=0;i<nbPts;i++){
		z_pos[i] = MinZ+ (0.2*(MaxZ-MinZ)) + ( ( 0.6*(MaxZ-MinZ) )/(nbPts-1.))*i;
		diso_pos[i] = 0.;
		for(unsigned int j=0;j<nbAt;j++){
			diso_pos[i] += gaussian(z_pos[i], AtData[j][2], sigma)*Disorder[j];
		}
		if(i == 0){
			MinDiso = diso_pos[i];
			MaxDiso = diso_pos[i];
			ZMaxDiso = z_pos[i];
		}else if( diso_pos[i] < MinDiso ){
			MinDiso = diso_pos[i];
		}else if( diso_pos[i] > MaxDiso ){
			MaxDiso = diso_pos[i];
			ZMaxDiso = z_pos[i];
		}
	}
	// normalize and shift the distrib to the lowest values
	NormFacGauss = 0.;
	for(unsigned int i=0;i<nbPts;i++){
		diso_pos[i] -= MinDiso;
		NormFacGauss += diso_pos[i];
	}
	double Mu_GB; // mean of Gaussian distrib
	double Sigma_GB; // stdev of Gaussian distrib
	Mu_GB = 0.;
	Sigma_GB = 0.;
	for(unsigned int i=0;i<nbPts;i++) Mu_GB += z_pos[i]*diso_pos[i]/NormFacGauss;
	for(unsigned int i=0;i<nbPts;i++) Sigma_GB += pow((z_pos[i]*diso_pos[i]/NormFacGauss)-Mu_GB,2.);
	cout << Mu_GB << " " << Sigma_GB << endl;
	Sigma_GB = pow(Sigma_GB,.5)/nbPts;
	cout << Mu_GB << " " << Sigma_GB << endl;
	double *EstDisoPos = new double[nbPts];
	for(unsigned int i=0;i<nbPts;i++) EstDisoPos[i] = gaussian(z_pos[i], Mu_GB, Sigma_GB);
	
	// write output file
	ofstream writefile(filename);
	writefile << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n" << nbAt << "\nITEM: BOX BOUNDS xy xz yz pp pp pp\n0 " << lx << " 0\n0 "<< ly << " 0\n0 " << lz << " 0\nITEM: ATOMS id type xu yu zu order disorder\n";
	for(unsigned int i=0;i<nbAt;i++) writefile << i+1 << " " << AtData[i][3] << " " << AtData[i][0] << " " << AtData[i][1] << " " << AtData[i][2] << " " << Order[i] << " " << Disorder[i] << "\n";
	writefile.close();

	ofstream writefile_distrib(filename_distrib);
	for(unsigned int i=0;i<nbPts;i++) writefile_distrib << z_pos[i] << " " << diso_pos[i]/NormFacGauss << " " << EstDisoPos[i] << "\n";
	writefile_distrib.close();
	
	// Free allocated memory
	delete[] z_pos;
	delete[] diso_pos;
	delete[] Calpha;
	delete[] AtSites;
	delete[] Order;
	delete[] Disorder;

	return 0;
}

double gaussian(double x, double mu, double sigma){
	return (1./(sigma*sqrt(M_PI*2.)))*exp(-pow(x-mu, 2.)/(2.*pow(sigma,2.)));
}

