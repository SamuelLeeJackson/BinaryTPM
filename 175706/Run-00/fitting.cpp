#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <stdlib.h>
using namespace std;

// Input parameters

const int num_lightcurves = 1;										// Number of data sets
const int num_wavelengths[num_lightcurves] = {3};				// Number of wavelengths in each data set
const int num_data_points = 48;									// Total number of data points
const int num_rand_data_points = 0;									// Number of randomly generated data points if performing Monte Carlo study (set to 0 if not)
const int minimum_thermal_inertia = 0;								// Minimum thermal inertia value to run model for
const int maximum_thermal_inertia = 400;							// Maximum thermal inertia value to run model for (set to minimum thermal inertia value if you want to run for a single and fixed thermal inertia value)
const int thermal_inertia_step = 10;								// Thermal inertia step size
const int num_roughness_steps = 21;									// Number of roughness steps (varies from 0 to 1 in this number of increments)
const int num_diameter_pv_bins = 20;								// Number of geometric albedo bins to bin the accepted clones into
const int NTAU = 2925;												// Number of modelling time steps
string obs_data_file = "175706-Run-00";								// Data input file name
string model_data_file = "1996FG3_binary_lowres";									// Shape model file name
int phase_offset0 = 1863;												// Phase offset constant in time steps (to adjust observation and modelling phasing)
double min_diameter = 1000.0;										// Minimum possible diameter of asteroid in metres
double max_diameter = 3000.0;										// Maximum possible diameter of asteroid in meteres
double diameter_step = 10.0;											// Diameter step in metres
double chi_square_threshold = 14.2;									// Chi-square fitting threshold (for 3 free parameters it is 3.53 for 1-sigma, 8.02 for 2-sigma, and 14.2 for 3-sigma)
double H = 17.833;													// Absolute visual magnitude of asteroid
double H_uncertainty = 0.024;											// Absolute visual magnitude uncertainty
double G = -0.041;													// Phase slope of asteroid
double G_uncertainty = 0.005;											// Phase slope uncertainty
int num_rand_magnitudes = 1;										// Number of randomly generated absolute magnitudes and phase slopes
double model_diameter =	1.70055;								// Unscaled effective diameter of asteroid shape model
double model_albedo = 0.01;											// Bond albedo of a smooth flat surface used in thermophysical modelling
int seed = 18;														// Seed value for random number generation
double emissivity_multiplier = 0.95/0.95;

// Global variables

const int num_model_TIs = (maximum_thermal_inertia-minimum_thermal_inertia)/thermal_inertia_step;


// Uniform random number generator


struct Ran {
	unsigned long long int u,v,w;
	Ran (unsigned long long int j) : v(4101842887655102017LL), w(1) {
		u = j^v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline unsigned long long int int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		unsigned long long int x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	inline unsigned int int32() { return (unsigned int)int64(); }
};


// Normal random number generator


struct Normaldev : Ran {
	double mu,sig;
	Normaldev(double mmu, double ssig, unsigned long long int i) : Ran (i), mu(mmu), sig(ssig) {}
	double dev() {
		double u,v,x,y,q;
		do {
			u = doub();
			v = 1.7156*(doub()-0.5);
			x = u - 0.449871;
			y = abs(v) + 0.386595;
			q = (x*x) + y*(0.19600*y-0.25472*x);
		} while (q > 0.27597
			&& (q > 0.27846 || (v*v) > -4.*log(u)*(u*u)));
		return mu + sig*v/u;
	}
};


// Class Obs_Data to store and manipulate observational thermal data


class Obs_Data {

private:

	int data_point_lightcurves[num_data_points];
	int data_point_wavelengths[num_data_points];
	int data_point_phases[num_data_points];
	int data_point_phase_uncertainties[num_data_points];
	double data_points[num_rand_data_points+1][num_data_points];
	double data_point_errors[num_data_points];

public:

	int &rtn_data_point_lightcurves(int n) { return data_point_lightcurves[n]; }
	int &rtn_data_point_wavelengths(int n) { return data_point_wavelengths[n]; }
	int &rtn_data_point_phases(int n) { return data_point_phases[n]; }
	int &rtn_data_point_phase_uncertainties(int n) { return data_point_phase_uncertainties[n]; }
	double &rtn_data_points(int n, int m) { return data_points[n][m]; }
	double &rtn_data_point_errors(int n) { return data_point_errors[n]; }

	void read_data() {
		string filename;
		filename=obs_data_file+".txt";
		ifstream read(filename.c_str());
		for (int i=0;i!=num_data_points;++i) {
			read >> data_point_lightcurves[i];
			read >> data_point_wavelengths[i];
			read >> data_point_phases[i];
			read >> data_point_phase_uncertainties[i];
			read >> data_points[0][i];
			read >> data_point_errors[i];
		}
	}

	void write_data() {
		string filename;
		filename=obs_data_file+"_test.txt";
		ofstream write(filename.c_str());
		for (int i=0;i!=num_data_points;++i) {
			write << data_point_lightcurves[i] << "\t";
			write << data_point_wavelengths[i] << "\t";
			write << data_point_phases[i] << "\t";
			write << data_point_phase_uncertainties[i] << "\t";
			write << data_points[0][i] << "\t";
			write << data_point_errors[i] << endl;
		}
	}

	void randomize_data() {
		Normaldev norm(0.0,1.0,seed);
		for (int i=0;i!=num_rand_data_points;++i) {
			for (int j=0;j!=num_data_points;++j) {
				norm.mu=data_points[0][j];
				norm.sig=data_point_errors[j];
				data_points[i+1][j]=norm.dev();
			}
		}
	}
	
};


// Class Model_Lcs to store and manipulate individual modelled thermal light-curves


class Model_Lcs {

private:

	int numwavelengths;
	double ***model_data_points;

public:

	void create_data_points() {
		model_data_points=(double ***)malloc(numwavelengths*sizeof(double **));
		for (int i=0;i!=numwavelengths;++i) {
			model_data_points[i]=(double **)malloc(2*sizeof(double *));
			for (int j=0;j!=2;++j ) {
				model_data_points[i][j]=(double *)malloc(NTAU*sizeof(double));
			}
		}
	}

	int &rtn_numwavelengths() { return numwavelengths; }
	double &rtn_model_data_points(int n, int m, int t) { return model_data_points[n][m][t]; }

};


// Class Model_Obs to control a set of modelled thermal light-curves


class Model_Obs {

private:

	double model_TIs[num_model_TIs+1];

	Model_Lcs **MLCS;

public:

	Model_Obs() {
		MLCS=(Model_Lcs **)malloc(num_lightcurves*sizeof(Model_Lcs *));
		for (int i=0;i!=num_lightcurves;++i) {
			MLCS[i]=(Model_Lcs *)malloc((num_model_TIs+1)*sizeof(Model_Lcs));
			for (int j=0;j!=(num_model_TIs+1);++j) {
				MLCS[i][j].rtn_numwavelengths()=num_wavelengths[i];
				MLCS[i][j].create_data_points();
			}
		}
	}

	Model_Lcs*beginmlcs(int n) const { return MLCS[n]; }
	Model_Lcs*endmlcs(int n) const { return MLCS[n]+num_model_TIs+1; }

	double &rtn_model_TIs(int n) { return model_TIs[n]; }

	void read_model_data() {
		
		string smooth_data_file;
		string rough_data_file;
		int t_inertia;
		int num_wavelengths;

		for (int i=0;i!=num_lightcurves;++i) {
			for (int j=0;j!=(num_model_TIs+1);++j) {

				t_inertia=minimum_thermal_inertia+(thermal_inertia_step*j);
				stringstream inert;
				stringstream lcurve;
				inert << t_inertia;
				lcurve << i;
				model_TIs[j]=t_inertia;

				smooth_data_file=model_data_file+"_smooth_flux("+inert.str()+"_"+lcurve.str()+").txt";
				rough_data_file=model_data_file+"_rough_flux("+inert.str()+"_"+lcurve.str()+").txt";
			
				ifstream smooth_read(smooth_data_file.c_str());
				ifstream rough_read(rough_data_file.c_str());

				num_wavelengths=MLCS[i][j].rtn_numwavelengths();

				for (int step=0;step!=NTAU;++step) {
					for (int k=0;k!=num_wavelengths;++k) {
						smooth_read >> MLCS[i][j].rtn_model_data_points(k,0,step);
						rough_read >> MLCS[i][j].rtn_model_data_points(k,1,step);
					}
				}

				for (int step=0;step!=NTAU;++step) {
					MLCS[i][j].rtn_model_data_points(0,0,step) *= emissivity_multiplier;
					MLCS[i][j].rtn_model_data_points(0,1,step) *= emissivity_multiplier;
				}

				cout << (i+1) << "\t" << t_inertia << endl;
	
			}

		}

	}

};


// Class Clone to store properties of test asteroid clones


class Clone {

private:

	int thermal_inertia_number;
	int roughness_number;

	int phase_offset;

	double Hv;
	double Gv;

	double diameter;
	double roughness;
	double pv;
	double effective_albedo;
	double albedo;
	double correction_factor;
	double thermal_inertia;

	double chi_square;
	bool accepted;

public:

	Clone() { phase_offset=phase_offset0; }

	int &rtn_thermal_inertia_number() { return thermal_inertia_number; }
	int &rtn_roughness_number() { return roughness_number; }

	int &rtn_phase_offset() { return phase_offset; }

	double &rtn_Hv() { return Hv; }
	double &rtn_Gv() { return Gv; }

	double &rtn_diameter() { return diameter; }
	double &rtn_roughness() { return roughness; }
	double &rtn_pv() { return pv; }
	double &rtn_effective_albedo() { return effective_albedo; }
	double &rtn_albedo() { return albedo; }
	double &rtn_correction_factor() { return correction_factor; }
	double &rtn_thermal_inertia() { return thermal_inertia; }

	double &rtn_chi_square() { return chi_square; }
	bool &rtn_accepted() { return accepted; }

};


// Class Asteroid to control all test clones and to store common and average properties


class Asteroid {

private:

	Clone *C;
	int num_clones;
	int best_clone_number;

	int diameter_bins[num_diameter_pv_bins];
	int pv_bins[num_diameter_pv_bins];
	int TI_bins[num_model_TIs+1];
	int roughness_bins[num_roughness_steps];

	double _model_diameter;

	double min_chi_square[num_rand_data_points+1];
	double average_phase_offset[num_rand_data_points+1];
	double phase_offset_uncertainty[num_rand_data_points+1];
	double average_Hv[num_rand_data_points+1];
	double Hv_uncertainty[num_rand_data_points+1];
	double average_Gv[num_rand_data_points+1];
	double Gv_uncertainty[num_rand_data_points+1];
	double average_diameter[num_rand_data_points+1];
	double diameter_uncertainty[num_rand_data_points+1];
	double average_pv[num_rand_data_points+1];
	double pv_uncertainty[num_rand_data_points+1];
	double average_roughness[num_rand_data_points+1];
	double roughness_uncertainty[num_rand_data_points+1];
	double average_TI[num_rand_data_points+1];
	double TI_uncertainty[num_rand_data_points+1];
	double average_FCF[num_rand_data_points+1];
	double FCF_uncertainty[num_rand_data_points+1];

	double min_pv;
	double max_pv;
	double pv_step;
	double diameter_step2;

public:

	Clone*beginc() const { return C; }
	Clone*endc() const { return C+num_clones; }

	Asteroid(int n) {
		num_clones=n;
		C=new Clone[n];
		_model_diameter=model_diameter;
	}

	int &rtn_num_clones() { return num_clones; }
	int &rtn_best_clone_number() { return best_clone_number; }

	int &rtn_diameter_bins(int n) { return diameter_bins[n]; }
	int &rtn_pv_bins(int n) { return pv_bins[n]; }
	int &rtn_TI_bins(int n) { return TI_bins[n]; }
	int &rtn_roughness_bins(int n) { return roughness_bins[n]; }

	double &rtn_model_diameter() { return _model_diameter; }

	double &rtn_min_chi_square(int n) { return min_chi_square[n]; }
	double &rtn_average_phase_offset(int n) { return average_phase_offset[n]; }
	double &rtn_phase_offset_uncertainty(int n) { return phase_offset_uncertainty[n]; }
	double &rtn_average_Hv(int n) { return average_Hv[n]; }
	double &rtn_Hv_uncertainty(int n) { return Hv_uncertainty[n]; }
	double &rtn_average_Gv(int n) { return average_Gv[n]; }
	double &rtn_Gv_uncertainty(int n) { return Gv_uncertainty[n]; }
	double &rtn_average_diameter(int n) { return average_diameter[n]; }
	double &rtn_diameter_uncertainty(int n) { return diameter_uncertainty[n]; }
	double &rtn_average_pv(int n) { return average_pv[n]; }
	double &rtn_pv_uncertainty(int n) { return pv_uncertainty[n]; }
	double &rtn_average_roughness(int n) { return average_roughness[n]; }
	double &rtn_roughness_uncertainty(int n) { return roughness_uncertainty[n]; }
	double &rtn_average_TI(int n) { return average_TI[n]; }
	double &rtn_TI_uncertainty(int n) { return TI_uncertainty[n]; }
	double &rtn_average_FCF(int n) { return average_FCF[n]; }
	double &rtn_FCF_uncertainty(int n) { return FCF_uncertainty[n]; }

	double &rtn_min_pv() { return min_pv; }
	double &rtn_max_pv() { return max_pv; }
	double &rtn_pv_step() { return pv_step; }
	double &rtn_diameter_step2() { return diameter_step2; }

};


// Function to calculate the albedo and flux correction factor for each test clone


void clone_asteroid_albedo(Asteroid &A) {

	Clone *C;
	int clone_num;
	double a,b,c;

#pragma omp parallel for private(C,a,b,c) schedule(dynamic)
	for (clone_num=0;clone_num<A.rtn_num_clones();clone_num++) {

		C=A.beginc()+clone_num;
		C->rtn_pv()=pow(((1329*pow(10.0,(-1.0*C->rtn_Hv()/5.0)))/(C->rtn_diameter()/1000.0)),2.0);
		C->rtn_effective_albedo()=C->rtn_pv()*(0.29+(0.684*C->rtn_Gv()));

		if (C->rtn_roughness()!=1.0) {
			a=0.5*(C->rtn_roughness()-1.0);
			b=1.0+(0.5*(C->rtn_effective_albedo()-C->rtn_roughness()));
			c=-1.0*C->rtn_effective_albedo();
			C->rtn_albedo()=(sqrt((b*b)-(4.0*a*c))-b)/(2.0*a);
		} else {
			C->rtn_albedo()=(2*C->rtn_effective_albedo())/(1.0+C->rtn_effective_albedo());
		}

		C->rtn_correction_factor()=(1.0-C->rtn_albedo())/(1.0-model_albedo);

	}

};


// Function to calculate the chi-sqaure value for each clone


void clone_chi_square(Asteroid &A, Obs_Data &OD, Model_Obs &MO, int spectrum_number) {

	Clone *C;
	int clone_num;
	int i;
	int flux_step;
	int phase_step;
	int best_phase_offset;
	int step;
	double model_flux;
	double avg_model_flux;
	double model_flux_uncertainty;
	double min_chi_square;

#pragma omp parallel for private(C, i, step, flux_step, phase_step, best_phase_offset, model_flux, avg_model_flux, model_flux_uncertainty, min_chi_square) schedule(dynamic)
	for (clone_num=0;clone_num<A.rtn_num_clones();clone_num++) {

		C=A.beginc()+clone_num;
		C->rtn_phase_offset()=phase_offset0;
		C->rtn_chi_square()=0.0;

		for (i=0;i!=num_data_points;++i) {

			if (OD.rtn_data_point_phase_uncertainties(i)==0) {

				flux_step=OD.rtn_data_point_phases(i)+C->rtn_phase_offset();
				if (flux_step>(NTAU-1)) { flux_step=flux_step-NTAU; }
				if (flux_step<0) { flux_step=flux_step+NTAU; }

				avg_model_flux=C->rtn_roughness()*(MO.beginmlcs(OD.rtn_data_point_lightcurves(i))+C->rtn_thermal_inertia_number())->rtn_model_data_points(OD.rtn_data_point_wavelengths(i),1,flux_step);
				avg_model_flux=avg_model_flux+((1.0-C->rtn_roughness())*(MO.beginmlcs(OD.rtn_data_point_lightcurves(i))+C->rtn_thermal_inertia_number())->rtn_model_data_points(OD.rtn_data_point_wavelengths(i),0,flux_step));
				avg_model_flux=1.0e-6*C->rtn_correction_factor()*pow((C->rtn_diameter()/A.rtn_model_diameter()),2.0)*avg_model_flux;
			
				model_flux_uncertainty=0.0;

			} else {

				model_flux=0.0;
				
				for (phase_step=(-1*OD.rtn_data_point_phase_uncertainties(i));phase_step!=(OD.rtn_data_point_phase_uncertainties(i)+1);++phase_step) {

					flux_step=OD.rtn_data_point_phases(i)+C->rtn_phase_offset()+phase_step;
					if (flux_step>(NTAU-1)) { flux_step=flux_step-NTAU; }
					if (flux_step<0) { flux_step=flux_step+NTAU; }

					model_flux=model_flux+(C->rtn_roughness()*(MO.beginmlcs(OD.rtn_data_point_lightcurves(i))+C->rtn_thermal_inertia_number())->rtn_model_data_points(OD.rtn_data_point_wavelengths(i),1,flux_step));
					model_flux=model_flux+((1.0-C->rtn_roughness())*(MO.beginmlcs(OD.rtn_data_point_lightcurves(i))+C->rtn_thermal_inertia_number())->rtn_model_data_points(OD.rtn_data_point_wavelengths(i),0,flux_step));

				}

				avg_model_flux=model_flux/((2*OD.rtn_data_point_phase_uncertainties(i))+1);
				model_flux=0.0;
				model_flux_uncertainty=0.0;

				for (phase_step=(-1*OD.rtn_data_point_phase_uncertainties(i));phase_step!=(OD.rtn_data_point_phase_uncertainties(i)+1);++phase_step) {

					flux_step=OD.rtn_data_point_phases(i)+C->rtn_phase_offset()+phase_step;
					if (flux_step>(NTAU-1)) { flux_step=flux_step-NTAU; }

					model_flux=C->rtn_roughness()*(MO.beginmlcs(OD.rtn_data_point_lightcurves(i))+C->rtn_thermal_inertia_number())->rtn_model_data_points(OD.rtn_data_point_wavelengths(i),1,flux_step);
					model_flux=model_flux+((1.0-C->rtn_roughness())*(MO.beginmlcs(OD.rtn_data_point_lightcurves(i))+C->rtn_thermal_inertia_number())->rtn_model_data_points(OD.rtn_data_point_wavelengths(i),0,flux_step));

					model_flux_uncertainty=model_flux_uncertainty+pow((avg_model_flux-model_flux),2.0);
								
				}

				model_flux_uncertainty=sqrt(model_flux_uncertainty/(2*OD.rtn_data_point_phase_uncertainties(i)));

				avg_model_flux=1.0e-6*C->rtn_correction_factor()*pow((C->rtn_diameter()/A.rtn_model_diameter()),2.0)*avg_model_flux;
				model_flux_uncertainty=1.0e-6*C->rtn_correction_factor()*pow((C->rtn_diameter()/A.rtn_model_diameter()),2.0)*model_flux_uncertainty;
		
			}

			C->rtn_chi_square()=C->rtn_chi_square()+(pow((OD.rtn_data_points(spectrum_number,i)-avg_model_flux),2.0)/(pow(OD.rtn_data_point_errors(i),2.0)+pow(model_flux_uncertainty,2.0)));

		}

		best_phase_offset=C->rtn_phase_offset();
		min_chi_square=C->rtn_chi_square();

		C->rtn_phase_offset()=best_phase_offset;
		C->rtn_chi_square()=min_chi_square;

	}

};


// Function to find the minimum chi-square value and associated clone


void minimum_chi_square(Asteroid &A, int spectrum_number) {

	Clone *C;
	C=A.beginc();
	A.rtn_min_chi_square(spectrum_number)=C->rtn_chi_square();
	A.rtn_best_clone_number()=0;

	for (int i=1;i!=A.rtn_num_clones();++i) {
		C=A.beginc()+i;
		if (C->rtn_chi_square()<A.rtn_min_chi_square(spectrum_number)) { 
			A.rtn_min_chi_square(spectrum_number)=C->rtn_chi_square();
			A.rtn_best_clone_number()=i;
		}
	}

	for (C=A.beginc();C!=A.endc();++C) {
		C->rtn_accepted()=false;
		if (C->rtn_chi_square()<(chi_square_threshold+A.rtn_min_chi_square(spectrum_number))) { C->rtn_accepted()=true; }
	}

};


// Function to perform statistics on accepted clones


void statistics(Asteroid &A, int spectrum_number) {
	
	for (int i=0;i!=num_diameter_pv_bins;++i) {
		A.rtn_diameter_bins(i)=0;
		A.rtn_pv_bins(i)=0;
	}
	for (int i=0;i!=(num_model_TIs+1);++i) {
		A.rtn_TI_bins(i)=0;
	}
	for (int i=0;i!=num_roughness_steps;++i) {
		A.rtn_roughness_bins(i)=0;
	}

	A.rtn_min_pv()=pow(((1329*pow(10.0,(-1.0*H/5.0)))/(max_diameter/1000.0)),2.0);
	A.rtn_max_pv()=pow(((1329*pow(10.0,(-1.0*H/5.0)))/(min_diameter/1000.0)),2.0);
	A.rtn_pv_step()=(A.rtn_max_pv()-A.rtn_min_pv())/num_diameter_pv_bins;
	A.rtn_diameter_step2()=(max_diameter-min_diameter)/num_diameter_pv_bins;

	A.rtn_average_phase_offset(spectrum_number)=0.0;
	A.rtn_average_Hv(spectrum_number)=0.0;
	A.rtn_average_Gv(spectrum_number)=0.0;
	A.rtn_average_diameter(spectrum_number)=0.0;
	A.rtn_average_pv(spectrum_number)=0.0;
	A.rtn_average_roughness(spectrum_number)=0.0;
	A.rtn_average_TI(spectrum_number)=0.0;
	A.rtn_average_FCF(spectrum_number)=0.0;

	int num_accepted_clones;
	num_accepted_clones=0;
	double pv_lower, pv_upper, diameter_lower, diameter_upper;

	for (Clone *C=A.beginc();C!=A.endc();++C) {
		if (C->rtn_accepted()==true) {

			for (int i=0;i!=num_diameter_pv_bins;++i) {
				pv_lower=A.rtn_min_pv()+(i*A.rtn_pv_step());
				pv_upper=A.rtn_min_pv()+((i+1)*A.rtn_pv_step());
				if (C->rtn_pv()>=pv_lower && C->rtn_pv()<pv_upper) { A.rtn_pv_bins(i)=A.rtn_pv_bins(i)+1; }
				diameter_lower=min_diameter+((i-0.5)*A.rtn_diameter_step2());
				diameter_upper=min_diameter+((i+0.5)*A.rtn_diameter_step2());
				if (C->rtn_diameter()>=diameter_lower && C->rtn_diameter()<diameter_upper) { A.rtn_diameter_bins(i)=A.rtn_diameter_bins(i)+1; }
			}
			
			A.rtn_TI_bins(C->rtn_thermal_inertia_number())=A.rtn_TI_bins(C->rtn_thermal_inertia_number())+1;
			A.rtn_roughness_bins(C->rtn_roughness_number())=A.rtn_roughness_bins(C->rtn_roughness_number())+1;

			A.rtn_average_phase_offset(spectrum_number)=A.rtn_average_phase_offset(spectrum_number)+C->rtn_phase_offset();
			A.rtn_average_Hv(spectrum_number)=A.rtn_average_Hv(spectrum_number)+C->rtn_Hv();
			A.rtn_average_Gv(spectrum_number)=A.rtn_average_Gv(spectrum_number)+C->rtn_Gv();
			A.rtn_average_diameter(spectrum_number)=A.rtn_average_diameter(spectrum_number)+C->rtn_diameter();
			A.rtn_average_pv(spectrum_number)=A.rtn_average_pv(spectrum_number)+C->rtn_pv();
			A.rtn_average_roughness(spectrum_number)=A.rtn_average_roughness(spectrum_number)+C->rtn_roughness();
			A.rtn_average_TI(spectrum_number)=A.rtn_average_TI(spectrum_number)+C->rtn_thermal_inertia();
			A.rtn_average_FCF(spectrum_number)=A.rtn_average_FCF(spectrum_number)+C->rtn_correction_factor();
			num_accepted_clones=num_accepted_clones+1; 
		
		}
	}

	A.rtn_average_phase_offset(spectrum_number)=A.rtn_average_phase_offset(spectrum_number)/num_accepted_clones;
	A.rtn_average_Hv(spectrum_number)=A.rtn_average_Hv(spectrum_number)/num_accepted_clones;
	A.rtn_average_Gv(spectrum_number)=A.rtn_average_Gv(spectrum_number)/num_accepted_clones;
	A.rtn_average_diameter(spectrum_number)=A.rtn_average_diameter(spectrum_number)/num_accepted_clones;
	A.rtn_average_pv(spectrum_number)=A.rtn_average_pv(spectrum_number)/num_accepted_clones;
	A.rtn_average_roughness(spectrum_number)=A.rtn_average_roughness(spectrum_number)/num_accepted_clones;
	A.rtn_average_TI(spectrum_number)=A.rtn_average_TI(spectrum_number)/num_accepted_clones;
	A.rtn_average_FCF(spectrum_number)=A.rtn_average_FCF(spectrum_number)/num_accepted_clones;

	A.rtn_phase_offset_uncertainty(spectrum_number)=0.0;
	A.rtn_Hv_uncertainty(spectrum_number)=0.0;
	A.rtn_Gv_uncertainty(spectrum_number)=0.0;
	A.rtn_diameter_uncertainty(spectrum_number)=0.0;
	A.rtn_pv_uncertainty(spectrum_number)=0.0;
	A.rtn_roughness_uncertainty(spectrum_number)=0.0;
	A.rtn_TI_uncertainty(spectrum_number)=0.0;
	A.rtn_FCF_uncertainty(spectrum_number)=0.0;

	for (Clone *C=A.beginc();C!=A.endc();++C) {
		if (C->rtn_accepted()==true) {
			A.rtn_phase_offset_uncertainty(spectrum_number)=A.rtn_phase_offset_uncertainty(spectrum_number)+pow((C->rtn_phase_offset()-A.rtn_average_phase_offset(spectrum_number)),2.0);
			A.rtn_Hv_uncertainty(spectrum_number)=A.rtn_Hv_uncertainty(spectrum_number)+pow((C->rtn_Hv()-A.rtn_average_Hv(spectrum_number)),2.0);
			A.rtn_Gv_uncertainty(spectrum_number)=A.rtn_Gv_uncertainty(spectrum_number)+pow((C->rtn_Gv()-A.rtn_average_Gv(spectrum_number)),2.0);
			A.rtn_diameter_uncertainty(spectrum_number)=A.rtn_diameter_uncertainty(spectrum_number)+pow((C->rtn_diameter()-A.rtn_average_diameter(spectrum_number)),2.0);
			A.rtn_pv_uncertainty(spectrum_number)=A.rtn_pv_uncertainty(spectrum_number)+pow((C->rtn_pv()-A.rtn_average_pv(spectrum_number)),2.0);
			A.rtn_roughness_uncertainty(spectrum_number)=A.rtn_roughness_uncertainty(spectrum_number)+pow((C->rtn_roughness()-A.rtn_average_roughness(spectrum_number)),2.0);
			A.rtn_TI_uncertainty(spectrum_number)=A.rtn_TI_uncertainty(spectrum_number)+pow((C->rtn_thermal_inertia()-A.rtn_average_TI(spectrum_number)),2.0);
			A.rtn_FCF_uncertainty(spectrum_number)=A.rtn_FCF_uncertainty(spectrum_number)+pow((C->rtn_correction_factor()-A.rtn_average_FCF(spectrum_number)),2.0);
		}
	}

	A.rtn_phase_offset_uncertainty(spectrum_number)=sqrt(A.rtn_phase_offset_uncertainty(spectrum_number)/(num_accepted_clones-1));
	A.rtn_Hv_uncertainty(spectrum_number)=sqrt(A.rtn_Hv_uncertainty(spectrum_number)/(num_accepted_clones-1));
	A.rtn_Gv_uncertainty(spectrum_number)=sqrt(A.rtn_Gv_uncertainty(spectrum_number)/(num_accepted_clones-1));
	A.rtn_diameter_uncertainty(spectrum_number)=sqrt(A.rtn_diameter_uncertainty(spectrum_number)/(num_accepted_clones-1));
	A.rtn_pv_uncertainty(spectrum_number)=sqrt(A.rtn_pv_uncertainty(spectrum_number)/(num_accepted_clones-1));
	A.rtn_roughness_uncertainty(spectrum_number)=sqrt(A.rtn_roughness_uncertainty(spectrum_number)/(num_accepted_clones-1));
	A.rtn_TI_uncertainty(spectrum_number)=sqrt(A.rtn_TI_uncertainty(spectrum_number)/(num_accepted_clones-1));
	A.rtn_FCF_uncertainty(spectrum_number)=sqrt(A.rtn_FCF_uncertainty(spectrum_number)/(num_accepted_clones-1));

};


// Function to write out statistics of best fitting clone


void write_best_clone_stats(Asteroid &A, Obs_Data &OD, Model_Obs &MO) {

	int flux_step;
	double model_flux;
	Clone *C=A.beginc()+A.rtn_best_clone_number();

	string write6name;
	write6name=model_data_file+"_best_clone_stats.txt";
	ofstream write6(write6name.c_str());

	write6 << "Chi-Square" << "\t" << C->rtn_chi_square() << endl;
	write6 << "Diameter" << "\t" << C->rtn_diameter() << endl;
	write6 << "Pv Albedo" << "\t" << C->rtn_pv() << endl;
	write6 << "Thermal Inertia" << "\t" << C->rtn_thermal_inertia() << endl;
	write6 << "Roughness" << "\t" << C->rtn_roughness() << endl;
	write6 << "H Value" << "\t" << C->rtn_Hv() << endl;
	write6 << "G Value" << "\t" << C->rtn_Gv() << endl;
	write6 << "Phase Offset" << "\t" << C->rtn_phase_offset() << endl;
	write6 << "FCF" << "\t" << C->rtn_correction_factor() << endl;
	
	for (int i=0;i!=num_data_points;++i) {
			
		flux_step=OD.rtn_data_point_phases(i)+C->rtn_phase_offset();
		if (flux_step>(NTAU-1)) { flux_step=flux_step-NTAU; }

		model_flux=C->rtn_roughness()*(MO.beginmlcs(OD.rtn_data_point_lightcurves(i))+C->rtn_thermal_inertia_number())->rtn_model_data_points(OD.rtn_data_point_wavelengths(i),1,flux_step);
		model_flux=model_flux+((1.0-C->rtn_roughness())*(MO.beginmlcs(OD.rtn_data_point_lightcurves(i))+C->rtn_thermal_inertia_number())->rtn_model_data_points(OD.rtn_data_point_wavelengths(i),0,flux_step));
		model_flux=1.0e-6*C->rtn_correction_factor()*pow((C->rtn_diameter()/A.rtn_model_diameter()),2.0)*model_flux;

		write6 << OD.rtn_data_point_phases(i) << "\t" << model_flux << "\t" << sqrt(((OD.rtn_data_points(0,i)-model_flux)/OD.rtn_data_point_errors(i))*((OD.rtn_data_points(0,i)-model_flux)/OD.rtn_data_point_errors(i))) << endl;

	}

};


// Main program loop


int main() {
	cout << "1" << endl;
	Obs_Data OD;
	OD.read_data();
	OD.write_data();
	if (num_rand_data_points!=0) { OD.randomize_data(); }
	cout << "2" << endl;

	Model_Obs MO;
	MO.read_model_data();
	cout << "3" << endl;

	int num_size_steps;
	num_size_steps=(int)((max_diameter-min_diameter)/diameter_step)+1;
	double roughness_step;
	roughness_step=1.0/(num_roughness_steps-1.0);
	int num_clones;
	num_clones=num_rand_magnitudes*(num_model_TIs+1)*num_size_steps*num_roughness_steps;
	cout << "4" << endl;

	Asteroid A(num_clones);
	cout << "5" << endl;
	
	Clone *C;
	int clone_number;
	clone_number=0;
	cout << "6" << endl;

	double G_rand, H_rand;
	Normaldev norm1(G,G_uncertainty,seed);
	Normaldev norm2(H,H_uncertainty,seed);
	cout << "7" << endl;

	for (int m=0;m!=num_rand_magnitudes;++m) {
		
		if (m==0) {
			G_rand=G;
			H_rand=H;
		} else {
			G_rand=norm1.dev();
			H_rand=norm2.dev();
		}

		for (int j=0;j!=(num_model_TIs+1);++j) {
			for (int k=0;k!=num_size_steps;++k) {
				for (int l=0;l!=num_roughness_steps;++l) {
					C=A.beginc()+clone_number;
					C->rtn_Gv()=G_rand;
					C->rtn_Hv()=H_rand;
					C->rtn_thermal_inertia_number()=j;
					C->rtn_thermal_inertia()=MO.rtn_model_TIs(j);
					C->rtn_diameter()=min_diameter+(k*diameter_step);
					C->rtn_roughness_number()=l;
					C->rtn_roughness()=l*roughness_step;
					clone_number=clone_number+1;
				}
			}
		}

	}
	cout << "8" << endl;

	clone_asteroid_albedo(A);
	cout << "9" << endl;

	clone_chi_square(A,OD,MO,0);
	cout << "10" << endl;
	minimum_chi_square(A,0);
	cout << "11" << endl;
	statistics(A,0);
	cout << "12" << endl;

	string writename;
	writename=model_data_file+"_fitting_results.txt";
	ofstream write(writename.c_str());
	write << "Spectrum Number" << "," << "Minimum Chi-Square" << "," << "Average Diameter" << "," << "Diameter Uncertainty" << "," << "Average pv" << "," << "pv Uncertainty" << "," << "Average Thermal Inertia" << "," << "Thermal Inertia Uncertainty" << "," << "Average Roughness" << "," << "Roughness Uncertainty" << "," << "Average H" << "," << "H Uncertainty" << "," << "Average G" << "," << "G Uncertainty" << "," << "Average Phase Offset" << "," << "Phase Offset Uncertainty" << "," << "Average FCF" << "," << "FCF Uncertainty" << endl;
	write << "0" << "," << A.rtn_min_chi_square(0) << "," << A.rtn_average_diameter(0) << "," << A.rtn_diameter_uncertainty(0) << "," << A.rtn_average_pv(0) << "," << A.rtn_pv_uncertainty(0) << "," << A.rtn_average_TI(0) << "," << A.rtn_TI_uncertainty(0) << "," << A.rtn_average_roughness(0) << "," << A.rtn_roughness_uncertainty(0) << "," << A.rtn_average_Hv(0) << "," << A.rtn_Hv_uncertainty(0) << "," << A.rtn_average_Gv(0) << "," << A.rtn_Gv_uncertainty(0) << "," << A.rtn_average_phase_offset(0) << "," << A.rtn_phase_offset_uncertainty(0) << "," << A.rtn_average_FCF(0) << "," << A.rtn_FCF_uncertainty(0) << endl;

	string write2name;
	write2name=model_data_file+"_thermal_inertia_distribution.txt";
	ofstream write2(write2name.c_str());
	write2 << "Thermal Inertia" << "\t" << "Frequency" << endl;
	for (int i=0;i!=(num_model_TIs+1);++i) {
		write2 << MO.rtn_model_TIs(i) << "\t" << A.rtn_TI_bins(i) << endl;
	}

	string write3name;
	write3name=model_data_file+"_roughness_distribution.txt";
	ofstream write3(write3name.c_str());
	write3 << "Roughness" << "\t" << "Frequency" << endl;
	for (int i=0;i!=num_roughness_steps;++i) {
		write3 << (roughness_step*i) << "\t" << A.rtn_roughness_bins(i) << endl;
	}

	string write4name;
	write4name=model_data_file+"_diameter_distribution.txt";
	ofstream write4(write4name.c_str());
	write4 << "Diameter Lower Bound" << "\t" << "Diameter Upper Bound" << "\t" << "Frequency" << endl;
	for (int i=0;i!=num_diameter_pv_bins;++i) {
		write4 << (min_diameter+((i-0.5)*A.rtn_diameter_step2())) << "\t" << (min_diameter+((i+0.5)*A.rtn_diameter_step2())) << "\t" << A.rtn_diameter_bins(i) << endl;
	}

	string write5name;
	write5name=model_data_file+"_pv_distribution.txt";
	ofstream write5(write5name.c_str());
	write5 << "Pv Lower Bound" << "\t" << "Pv Upper Bound" << "\t" << "Frequency" << endl;
	for (int i=0;i!=num_diameter_pv_bins;++i) {
		write5 << (A.rtn_min_pv()+(i*A.rtn_pv_step())) << "\t" << (A.rtn_min_pv()+((i+1)*A.rtn_pv_step())) << "\t" << A.rtn_pv_bins(i) << endl;
	}

	write_best_clone_stats(A,OD,MO);

	string write7name;
	write7name=model_data_file+"_accepted_clones.txt";
	ofstream write7(write7name.c_str());
	for (Clone *C=A.beginc();C!=A.endc();++C) {
		if (C->rtn_accepted()==1) {
			write7 << C->rtn_chi_square() << "\t";
			write7 << C->rtn_diameter() << "\t";
			write7 << C->rtn_pv() << "\t";
			write7 << C->rtn_thermal_inertia() << "\t";
			write7 << C->rtn_roughness() << "\t";
			write7 << C->rtn_correction_factor() << endl;
		}
	}

	string write8name;
	write8name = model_data_file + "_all_clones.txt";
	ofstream write8(write8name.c_str());
	for (Clone* C = A.beginc(); C != A.endc(); ++C) {
		write8 << C->rtn_chi_square() << "\t";
		write8 << C->rtn_diameter() << "\t";
		write8 << C->rtn_pv() << "\t";
		write8 << C->rtn_thermal_inertia() << "\t";
		write8 << C->rtn_roughness() << "\t";
		write8 << C->rtn_correction_factor() << endl;
	}

	if (num_rand_data_points!=0) {	
		for (int i=1;i!=num_rand_data_points+1;++i) {
			clone_chi_square(A,OD,MO,i);
			minimum_chi_square(A,i);
			statistics(A,i);
			write << i << "\t" << A.rtn_min_chi_square(i) << "\t" << A.rtn_average_diameter(i) << "\t" << A.rtn_diameter_uncertainty(i) << "\t" << A.rtn_average_pv(i) << "\t" << A.rtn_pv_uncertainty(i) << "\t" << A.rtn_average_TI(i) << "\t" << A.rtn_TI_uncertainty(i) << "\t" << A.rtn_average_roughness(i) << "\t" << A.rtn_roughness_uncertainty(i) << "\t" << A.rtn_average_Hv(i) << "\t" << A.rtn_Hv_uncertainty(i) << "\t" << A.rtn_average_Gv(i) << "\t" << A.rtn_Gv_uncertainty(i) << endl;
		}
	}

};
