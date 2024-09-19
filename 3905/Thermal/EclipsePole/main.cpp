#include <omp.h>																															
#include <iostream>		
#include <string>
#include <fstream>																												
#include <sstream>																														
#include <cmath>																														
#include <time.h>
#include <stdio.h>
#include <cstdlib>
using namespace std;

const double PI = 4.0 * atan(1.0);									// Value of PI
double AU = 149.6e6;												// Astronomical unit in kilometres
double cvel = 299792458.0;											// Speed of light
float SIGMA = 5.67e-8f;												// Stefan-Boltzmann constant

// Input Parameters																																												

const int ASTNUMVTX = 3204;											// Number of shape model vertices
const int ASTNUMFACE = 6400;										// Number of shape model facets
string shape_filename = "doppler_binary";						// Shape model input file name
int load_shape_selfheating = 1;										// Load shape selfheating (= 1) or not (= 0)
const int SFCNUMVTX = 61;											// Number of surface roughness model vertices
const int SFCNUMFACE = 100;											// Number of surface roughness model facets
string roughness_filename = "new_crater2";							// Surface roughness model input file name
double polar0 = 1.977789;												// Heliocentric distance in AU
double polar1 = ((11.3684 / 180.0) * PI);								// Heliocentric ecliptic latitude in radians of asteroid position (use if polar_or_cart = 1)
double polar2 = ((110.6920 / 180.0) * PI);								// Heliocentric ecliptic longitude in radians of asteroid position (use if polar_or_cart = 1)

double semimajor_axis = 1.698;										// Semimajor axis in AU of asteroid orbit
double eccentricity = 0.508;										// Eccentricity of asteroid orbit
double inclination = 12.170;										// Inclination in degrees of asteroid orbit
double ascending_node = 356.652;									// Ascending node in degrees of asteroid orbit
double argument_periapsis = 224.658;								// Argument of periapsis in degrees of asteroid orbit

double rot0 = ((60.0 / 180.0) * PI);										// Heliocentric ecliptic latitude in radians of asteroid rotation pole
double rot1 = ((221.0 / 180.0) * PI);									// Heliocentric ecliptic longitude in radians of asteroid rotation pole

float thermal_inertia0 = 340.0f;
float r0 = 1.06f;
float alpha = -2.25f;

const int NORB = 36;												// Number of orbital steps (if running an orbital study model)
const int NOBL = 9;													// Number of obliquity steps (if running an obliquity study model)
const int NTAU = 650;												// Number of model time steps (multiples of 50)
const int NX = 56;													// Number of model depth steps
int shape_shadowing = 1;											// Perform shape shadowing calculations (= 1) or not (= 0)
int roughness_thermal_model = 1;									// Perform surface roughness thermal modelling (= 1) or not (= 0)
int multiple_scattering = 1;										// Perform multiple scattering of sunlight calculations (= 1) or not (= 0)
float TACC = 0.001f;		// 0.001								// Thermal model convergence goal in Kelvin
float scattering_accuracy = 0.0001f;								// Multiple scattering calculation accuracy in W/m^2
float EPS = 0.9f;													// Surface emissivity
float rotation_period = 50.826f * 60.0f * 60.0f;						// Asteroid rotation period in seconds
float bond_albedo = 0.137f;											// Bond albedo of a smooth flat surface
float minimum_thermal_inertia = 0.0f;								// Minimum thermal inertia value to run model for
float maximum_thermal_inertia = 500.0f;								// Maximum thermal inertia value to run model for (set to minimum thermal inertia value if you want to run for a single and fixed thermal inertia value)
float thermal_inertia_step = 10.0f;									// Thermal inertia step size
int perform_yark_yorp = 0;											// Perform Yarkovsky/YORP calculations (= 1) or not (= 0)
const int m = 10;													// Number of partial shadow divisions (splits facets into m*m smaller facets)
const int num_threads = 64;											// Number of parallel threads to use

int single_run = 1;													// Perform single model run (= 1) or not (= 0)
int orbital_run = 0;												// Perform orbital model run (= 1) or not (= 0)
int obliquity_run = 0;												// Perform obliquity model run (= 1) or (= 0)

// Output Parameters																																	

int write_shape_shadow = 1;											// Write out shape shadows (= 1) or not (= 0)
int write_shape_illumination = 1;									// Write out shape illumination (= 1) or not (= 0)
int write_shape_temperature = 1;									// Write out shape temperature (= 1) or not (= 0)
int write_roughness_temperature = 1;								// Write out surface roughness temperature (= 1) or not (= 0)
int write_yarkovsky = 0;											// Write out Yarkovsky force vectors (= 1) or not (= 0)
int write_yorp = 0;													// Write out YORP torque vectors (= 1) or not (= 0)

float thermal_inertia[98];


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}



// Class AstVertex defining position of a single vertex of the global asteroid shape model


class AstVertex {

private:

	double _equpos[3];
	double _rotequpos[NTAU][3];
	double _heliopos[NTAU][3];

public:

	friend std::istream& operator>> (std::istream& i, AstVertex& A) {
		return i >> A._equpos[0] >> A._equpos[1] >> A._equpos[2];
	}

	double& rtnequpos(int n) { return _equpos[n]; }
	double& rtnrotequpos(int t, int n) { return _rotequpos[t][n]; }
	double& rtnheliopos(int t, int n) { return _heliopos[t][n]; }

	~AstVertex() {}

};


// Class SubAstFacet to contain calculation data of surface roughness facets


class SubAstFacet {

private:

	float _shadow[NTAU];
	float _illumangle[NTAU];
	float _totalillumination[NTAU];
	float _temp[NX + 1];
	float _stemp[2][NTAU];

public:

	SubAstFacet() {
		for (int step = 0; step != NTAU; ++step) {
			_shadow[step] = 1.0f;
			_totalillumination[step] = 0.0f;
			_stemp[0][step] = 0.0f;
			_stemp[1][step] = 0.0f;
		}
	}

	void reset_subfacet_conditions() {
		for (int step = 0; step != NTAU; ++step) {
			_shadow[step] = 1.0f;
			_totalillumination[step] = 0.0f;
			_stemp[0][step] = 0.0f;
			_stemp[1][step] = 0.0f;
		}
	}

	float& rtnshadow(int t) { return _shadow[t]; }
	float& rtnillumangle(int t) { return _illumangle[t]; }
	float& rtntotalillumination(int t) { return _totalillumination[t]; }
	float& rtntemp(int n) { return _temp[n]; }
	float& rtnstemp(int n, int t) { return _stemp[n][t]; }

	~SubAstFacet() {}

};


// Class DepAstFacet to contain viewfactor data of dependent facets of global shape model


class DepAstFacet {

private:

	int _depastfacetnum;
	float _viewfactor[SFCNUMFACE + 1];
	float _viewfactor_vector[3];

public:

	int& rtndepastfacetnum() { return _depastfacetnum; }
	float& rtnviewfactor(int n) { return _viewfactor[n]; }
	float& rtnviewfactor_vector(int n) { return _viewfactor_vector[n]; }

	~DepAstFacet() {}

};


// Class AstFacet defining surface properties of a single surface facet of the global asteroid shape model


class AstFacet {

private:

	bool _operated;
	bool _astconverged;
	bool _sfcconverged;

	int _astvertices[3];
	int _astdependentfacets;
	int _shadow[NTAU];

	float _area;
	float _inert;
	float _albedo;
	float _sfcillumvector[NTAU][3];
	float _totalillumination[NTAU];
	float _temp[NX + 1];
	float _stemp[2][NTAU];

	double _midpos[NTAU][3];
	double _heliomidpos[NTAU][3];
	double _illumvector[NTAU][3];
	double _normal[NTAU][3];
	double _illumangle[NTAU];
	double _imatrix0[3][3];
	double _imatrix[NTAU][3][3];
	double _photonmomvector[3];
	double _radvector[NTAU][3];

	double _helsmoothsolarpressurevec[3];
	double _helsmoothreflectedpressurevec[3];
	double _helsmooththermalpressurevec[3];

	double _helroughsolarpressurevec[3];
	double _helroughreflectedpressurevec[3];
	double _helroughthermalpressurevec[3];

	double _smoothsolartorquevec[3];
	double _smoothreflectedtorquevec[3];
	double _smooththermaltorquevec[3];

	double _roughsolartorquevec[3];
	double _roughreflectedtorquevec[3];
	double _roughthermaltorquevec[3];

	SubAstFacet* SF;
	DepAstFacet* DAF;

public:

	SubAstFacet* beginsf() const { return SF; }
	SubAstFacet* endsf() const { return SF + SFCNUMFACE; }

	DepAstFacet* begindaf() const { return DAF; }
	DepAstFacet* enddaf() const { return DAF + _astdependentfacets; }

	AstFacet() {
		_astdependentfacets = 0;
		_inert = minimum_thermal_inertia;
		_albedo = bond_albedo;
		_photonmomvector[0] = 0.0, _photonmomvector[1] = 0.0, _photonmomvector[2] = 0.0;
	}

	void createsubfacets() {
		SF = new SubAstFacet[SFCNUMFACE];
	}

	void createdepfacets() {
		DAF = new DepAstFacet[_astdependentfacets];
	}

	void reset_facet_status() {
		_operated = false;
		_astconverged = false;
		_sfcconverged = false;
	}

	void reset_facet_conditions() {
		for (int step = 0; step != NTAU; ++step) {
			_shadow[step] = 1;
			_totalillumination[step] = 0.0f;
			_stemp[0][step] = 0.0f;
			_stemp[1][step] = 0.0f;
		}
		if (roughness_thermal_model == 1) {
			for (SubAstFacet* C = SF; C != SF + SFCNUMFACE; ++C) {
				C->reset_subfacet_conditions();
			}
		}
	}

	void reset_facet_yark_yorp() {
		_helsmoothsolarpressurevec[0] = 0.0, _helsmoothsolarpressurevec[1] = 0.0, _helsmoothsolarpressurevec[2] = 0.0;
		_helsmoothreflectedpressurevec[0] = 0.0, _helsmoothreflectedpressurevec[1] = 0.0, _helsmoothreflectedpressurevec[2] = 0.0;
		_helsmooththermalpressurevec[0] = 0.0, _helsmooththermalpressurevec[1] = 0.0, _helsmooththermalpressurevec[2] = 0.0;

		_helroughsolarpressurevec[0] = 0.0, _helroughsolarpressurevec[1] = 0.0, _helroughsolarpressurevec[2] = 0.0;
		_helroughreflectedpressurevec[0] = 0.0, _helroughreflectedpressurevec[1] = 0.0, _helroughreflectedpressurevec[2] = 0.0;
		_helroughthermalpressurevec[0] = 0.0, _helroughthermalpressurevec[1] = 0.0, _helroughthermalpressurevec[2] = 0.0;

		_smoothsolartorquevec[0] = 0.0, _smoothsolartorquevec[1] = 0.0, _smoothsolartorquevec[2] = 0.0;
		_smoothreflectedtorquevec[0] = 0.0, _smoothreflectedtorquevec[1] = 0.0, _smoothreflectedtorquevec[2] = 0.0;
		_smooththermaltorquevec[0] = 0.0, _smooththermaltorquevec[1] = 0.0, _smooththermaltorquevec[2] = 0.0;

		_roughsolartorquevec[0] = 0.0, _roughsolartorquevec[1] = 0.0, _roughsolartorquevec[2] = 0.0;
		_roughreflectedtorquevec[0] = 0.0, _roughreflectedtorquevec[1] = 0.0, _roughreflectedtorquevec[2] = 0.0;
		_roughthermaltorquevec[0] = 0.0, _roughthermaltorquevec[1] = 0.0, _roughthermaltorquevec[2] = 0.0;
	}

	friend std::istream& operator>> (std::istream& i, AstFacet& A) {
		return i >> A._astvertices[0] >> A._astvertices[1] >> A._astvertices[2];
	}

	bool& rtnoperated() { return _operated; }
	bool& rtnastconverged() { return _astconverged; }
	bool& rtnsfcconverged() { return _sfcconverged; }

	int& rtnastvertices(int n) { return _astvertices[n]; }
	int& rtnastdependentfacets() { return _astdependentfacets; }
	int& rtnshadow(int t) { return _shadow[t]; }

	float& rtnarea() { return _area; }
	float& rtninert() { return _inert; }
	float& rtnalbedo() { return _albedo; }
	float& rtnsfcillumvector(int t, int n) { return _sfcillumvector[t][n]; }
	float& rtntotalillumination(int t) { return _totalillumination[t]; }
	float& rtntemp(int n) { return _temp[n]; }
	float& rtnstemp(int n, int t) { return _stemp[n][t]; }

	double& rtnmidpos(int t, int n) { return _midpos[t][n]; }
	double& rtnheliomidpos(int t, int n) { return _heliomidpos[t][n]; }
	double& rtnillumvector(int t, int n) { return _illumvector[t][n]; }
	double& rtnnormal(int t, int n) { return _normal[t][n]; }
	double& rtnillumangle(int t) { return _illumangle[t]; }
	double& rtnimatrix0(int n, int m) { return _imatrix0[n][m]; }
	double& rtnimatrix(int t, int n, int m) { return _imatrix[t][n][m]; }
	double& rtnphotonmomvector(int n) { return _photonmomvector[n]; }
	double& rtnradvector(int t, int n) { return _radvector[t][n]; }

	double& rtnhelsmoothsolarpressurevec(int n) { return _helsmoothsolarpressurevec[n]; }
	double& rtnhelsmoothreflectedpressurevec(int n) { return _helsmoothreflectedpressurevec[n]; }
	double& rtnhelsmooththermalpressurevec(int n) { return _helsmooththermalpressurevec[n]; }

	double& rtnhelroughsolarpressurevec(int n) { return _helroughsolarpressurevec[n]; }
	double& rtnhelroughreflectedpressurevec(int n) { return _helroughreflectedpressurevec[n]; }
	double& rtnhelroughthermalpressurevec(int n) { return _helroughthermalpressurevec[n]; }

	double& rtnsmoothsolartorquevec(int n) { return _smoothsolartorquevec[n]; }
	double& rtnsmoothreflectedtorquevec(int n) { return _smoothreflectedtorquevec[n]; }
	double& rtnsmooththermaltorquevec(int n) { return _smooththermaltorquevec[n]; }

	double& rtnroughsolartorquevec(int n) { return _roughsolartorquevec[n]; }
	double& rtnroughreflectedtorquevec(int n) { return _roughreflectedtorquevec[n]; }
	double& rtnroughthermaltorquevec(int n) { return _roughthermaltorquevec[n]; }

	~AstFacet() {}

};


// Class Asteroid to control the asteroid global shape model


class Asteroid {

private:

	int NV, NF;
	int _numconverged;
	int _converged;
	int _sfcnumconverged;
	int _sfcconverged;

	float _rotperiod;
	float _solar;

	double _semimajor_axis;
	double _eccentricity;
	double _inclination;
	double _ascending_node;
	double _argument_periapsis;

	double _polcentre[3];
	double _carcentre[3];
	double _rotvector[3];
	double _imatrix[3][3];

	AstVertex* V;
	AstFacet* F;

public:

	AstVertex* beginv() const { return V; }
	AstVertex* endv() const { return V + NV; }

	AstFacet* beginf() const { return F; }
	AstFacet* endf() const { return F + NF; }

	Asteroid(int n, int m) : NV(n), NF(m) {
		V = new AstVertex[n];
		F = new AstFacet[m];
		_polcentre[0] = polar0, _polcentre[1] = polar1, _polcentre[2] = polar2;
		_rotvector[0] = rot0, _rotvector[1] = rot1;
		_semimajor_axis = semimajor_axis;
		_eccentricity = eccentricity;
		_inclination = inclination;
		_ascending_node = ascending_node;
		_argument_periapsis = argument_periapsis;
		_rotperiod = rotation_period;
	}

	void readmodel() {
		string filename;
		filename = shape_filename + ".txt";
		ifstream read(filename.c_str());
		for (AstVertex* i = V; i != V + NV; ++i) {
			read >> *i;
		}
		for (AstFacet* j = F; j != F + NF; ++j) {
			read >> *j;
		}
	}

	void invertmodel() {
		for (AstVertex* i = V; i != V + NV; ++i) {
			i->rtnequpos(2) = -1.0 * i->rtnequpos(2);
		}
	}

	void readselfheating() {
		string filename;
		filename = shape_filename + "_selfheating_list.txt";
		ifstream read(filename.c_str());
		for (AstFacet* k = F; k != F + NF; ++k) {
			read >> k->rtnastdependentfacets();
			if (k->rtnastdependentfacets() != 0) {
				k->createdepfacets();
				for (DepAstFacet* l = k->begindaf(); l != k->enddaf(); ++l) {
					read >> l->rtndepastfacetnum();
				}
				for (DepAstFacet* l = k->begindaf(); l != k->enddaf(); ++l) {
					read >> l->rtnviewfactor(0);
				}
			}
		}
	}

	void readselfheating_vectors() {
		string filename;
		filename = shape_filename + "_selfheating_vector_list.txt";
		ifstream read(filename.c_str());
		for (AstFacet* k = F; k != F + NF; ++k) {
			if (k->rtnastdependentfacets() != 0) {
				for (DepAstFacet* l = k->begindaf(); l != k->enddaf(); ++l) {
					read >> l->rtnviewfactor_vector(0);
					read >> l->rtnviewfactor_vector(1);
					read >> l->rtnviewfactor_vector(2);
				}
			}
		}
	}

	void readphotonmomvectors() {
		string filename;
		filename = shape_filename + "_photovectors.txt";
		ifstream read(filename.c_str());
		for (AstFacet* k = F; k != F + NF; ++k) {
			read >> k->rtnphotonmomvector(0);
			read >> k->rtnphotonmomvector(1);
			read >> k->rtnphotonmomvector(2);
		}
	}

	void readthermalinertia() {
		string filename;
		filename = shape_filename + "_thermal_inertia.txt";
		ifstream read(filename.c_str());
		for (AstFacet* k = F; k != F + NF; ++k) {
			read >> k->rtninert();
		}
	}

	void addsurfaceroughness() {
		for (AstFacet* i = F; i != F + NF; ++i) {
			i->createsubfacets();
		}
	}

	void reset_asteroid_status() {
		_converged = 0;
		_sfcconverged = 0;
		_numconverged = 0;
		_sfcnumconverged = 0;
		for (AstFacet* k = F; k != F + NF; ++k) {
			k->reset_facet_status();
		}
	}

	void reset_asteroid_conditions() {
		for (AstFacet* k = F; k != F + NF; ++k) {
			k->reset_facet_conditions();
		}
	}

	void reset_asteroid_yark_yorp() {
		for (AstFacet* k = F; k != F + NF; ++k) {
			k->reset_facet_yark_yorp();
		}
	}

	int& rtnnumvtx() { return NV; }
	int& rtnnumface() { return NF; }
	int& rtnnumconverged() { return _numconverged; }
	int& rtnconverged() { return _converged; }
	int& rtnsfcnumconverged() { return _sfcnumconverged; }
	int& rtnsfcconverged() { return _sfcconverged; }

	float& rtnrotperiod() { return _rotperiod; }
	float& rtnsolar() { return _solar; }

	double& rtnsemimajor_axis() { return _semimajor_axis; }
	double& rtneccentricity() { return _eccentricity; }
	double& rtninclination() { return _inclination; }
	double& rtnascending_node() { return _ascending_node; }
	double& rtnargument_periapsis() { return _argument_periapsis; }

	double& rtnpolcentre(int n) { return _polcentre[n]; }
	double& rtncarcentre(int n) { return _carcentre[n]; }
	double& rtnrotvector(int n) { return _rotvector[n]; }
	double& rtnimatrix(int n, int m) { return _imatrix[n][m]; }

};


// Class SfcVertex defining position of a single vertex of the surface roughness shape model


class SfcVertex {

private:

	float _pos[3];

public:

	friend std::istream& operator>> (std::istream& i, SfcVertex& A) {
		return i >> A._pos[0] >> A._pos[1] >> A._pos[2];
	}

	float& rtnpos(int n) { return _pos[n]; }

	~SfcVertex() {}

};


// Class DepSfcFacet to contain viewfactor data of dependent facets of surface roughness shape model


class DepSfcFacet {

private:

	int _depsfcfacetnum;
	float _viewfactor;

public:

	int& rtndepsfcfacetnum() { return _depsfcfacetnum; }
	float& rtnviewfactor() { return _viewfactor; }

	~DepSfcFacet() {}

};


// Class SfcFacet defining surface properties of a single surface facet of the surface roughness shape model


class SfcFacet {

private:

	int _sfcvertices[3];
	int _sfcdependentfacets;

	float _midpos[3];
	float _normcalc[3][3];
	float _normal[3];
	float _trivector[2][3];
	float _area;

	double _photonmomvector[3];
	double _radvector[3];

	DepSfcFacet* DSF;

public:

	DepSfcFacet* beginf() const { return DSF; }
	DepSfcFacet* endf() const { return DSF + _sfcdependentfacets; }

	friend std::istream& operator>> (std::istream& i, SfcFacet& A) {
		return i >> A._sfcvertices[0] >> A._sfcvertices[1] >> A._sfcvertices[2];
	}

	void createdepfacets() {
		DSF = new DepSfcFacet[_sfcdependentfacets];
	}

	int& rtnsfcvertices(int n) { return _sfcvertices[n]; }
	int& rtnsfcdependentfacets() { return _sfcdependentfacets; }

	float& rtnmidpos(int n) { return _midpos[n]; }
	float& rtnnormcalc(int n, int m) { return _normcalc[n][m]; }
	float& rtnnormal(int n) { return _normal[n]; }
	float& rtntrivector(int n, int m) { return _trivector[n][m]; }
	float& rtnarea() { return _area; }

	double& rtnphotonmomvector(int n) { return _photonmomvector[n]; }
	double& rtnradvector(int n) { return _radvector[n]; }

	~SfcFacet() {}

};


// Class Surface to control the surface roughness shape model


class Surface {

private:

	int NV, NF;
	float _totalarea;
	SfcVertex* SV;
	SfcFacet* SF;

public:

	Surface(int n, int m) : NV(n), NF(m) {
		SV = new SfcVertex[n];
		SF = new SfcFacet[m];
	}

	SfcVertex* beginv() const { return SV; }
	SfcVertex* endv() const { return SV + NV; }

	SfcFacet* beginf() const { return SF; }
	SfcFacet* endf() const { return SF + NF; }

	void readmodel() {
		string filename;
		filename = roughness_filename + ".txt";
		ifstream read(filename.c_str());
		for (SfcVertex* i = SV; i != SV + NV; ++i) {
			read >> *i;
		}
		for (SfcFacet* j = SF; j != SF + NF; ++j) {
			read >> *j;
		}
	}

	void readselfheating() {
		string filename;
		filename = roughness_filename + "_selfheating_list.txt";
		ifstream read(filename.c_str());
		for (SfcFacet* k = SF; k != SF + NF; ++k) {
			read >> k->rtnsfcdependentfacets();
			if (k->rtnsfcdependentfacets() != 0) {
				k->createdepfacets();
				for (DepSfcFacet* l = k->beginf(); l != k->endf(); ++l) {
					read >> l->rtndepsfcfacetnum();
				}
				for (DepSfcFacet* l = k->beginf(); l != k->endf(); ++l) {
					read >> l->rtnviewfactor();
				}
			}
		}
	}

	void readphotonmomvectors() {
		string filename;
		filename = roughness_filename + "_photovectors.txt";
		ifstream read(filename.c_str());
		for (SfcFacet* k = SF; k != SF + NF; ++k) {
			read >> k->rtnphotonmomvector(0);
			read >> k->rtnphotonmomvector(1);
			read >> k->rtnphotonmomvector(2);
		}
	}

	int& rtnnumvtx() { return NV; }
	int& rtnnumface() { return NF; }
	float& rtntotalarea() { return _totalarea; }

};


// Structure GPU_SfcVertex defining position of a single vertex of the surface roughness model


struct GPU_SfcVertex {

	float pos[3];

};


// Structure GPU_SfcFacet defining surface properties of a single surface facet of the surface roughness model


struct GPU_SfcFacet {

	int sfcvertices[3];
	float midpos[3];
	float trivector[2][3];
	float normcalc[3][3];
	float normal[3];

};


// Class GPU_Surface to control the surface roughness model


struct GPU_Surface {

	int NV, NF;
	GPU_SfcVertex* SV;
	GPU_SfcFacet* SF;

	GPU_Surface(int n, int m) : NV(n), NF(m) {
		SV = new GPU_SfcVertex[n];
		SF = new GPU_SfcFacet[m];
	}

};


// Structure Facet_Properties containing facet properties


struct Facet_Properties {

	bool sfc_converged;
	float albedo;
	float thermal_inertia;

};


// Structure Roughness_Shadow containing roughness facet shadow and illumination angle information


struct Roughness_Shadow {

	int shadow[NTAU];
	float sfc_shadow[NTAU][SFCNUMFACE];
	float sfc_illumvector[NTAU][3];
	float sfc_illumangle[NTAU][SFCNUMFACE];

};


// Structure Roughness_Illumination containing roughness facet illumination flux information


struct Roughness_Illumination {

	float sfc_illumination[NTAU][SFCNUMFACE];

};


// Structure Roughness_Globalillumination containing roughness global-illumination flux information


struct Roughness_Globalillumination {

	float sfc_globalillumination[NTAU][SFCNUMFACE];

};


// Structure Roughness_Globalheating containing roughness facet global-selfheating flux information


struct Roughness_Globalheating {

	float sfc_globalheating[NTAU][SFCNUMFACE];

};


// Structure Roughness_Temperatures containing roughness facet temperature information


struct Roughness_Temperatures {

	float sfc_stemp0[NTAU][SFCNUMFACE];
	float sfc_stemp1[NTAU][SFCNUMFACE];
	float sfc_temp[NX + 1][SFCNUMFACE];

};


// Determine orbital position geometry at a specific true anomaly step


void orbital_position_geometry(Asteroid& A, int step) {

	double TA;
	TA = 360.0 * step / NORB;

	A.rtnpolcentre(0) = (A.rtnsemimajor_axis() * AU * (1.0 - (A.rtneccentricity() * A.rtneccentricity()))) / (1.0 + (A.rtneccentricity() * cos(TA * (PI / 180.0))));
	A.rtncarcentre(0) = A.rtnpolcentre(0) * ((cos(A.rtnascending_node() * (PI / 180.0)) * cos((A.rtnargument_periapsis() + TA) * (PI / 180.0))) - (sin(A.rtnascending_node() * (PI / 180.0)) * sin((A.rtnargument_periapsis() + TA) * (PI / 180.0)) * cos(A.rtninclination() * (PI / 180.0))));
	A.rtncarcentre(1) = A.rtnpolcentre(0) * ((sin(A.rtnascending_node() * (PI / 180.0)) * cos((A.rtnargument_periapsis() + TA) * (PI / 180.0))) + (cos(A.rtnascending_node() * (PI / 180.0)) * sin((A.rtnargument_periapsis() + TA) * (PI / 180.0)) * cos(A.rtninclination() * (PI / 180.0))));
	A.rtncarcentre(2) = A.rtnpolcentre(0) * sin((A.rtnargument_periapsis() + TA) * (PI / 180.0)) * sin(A.rtninclination() * (PI / 180.0));

};


// Determine surface roughness geometry


void surfacegeometry(Surface& S, Asteroid& A) {

	float n;

	S.rtntotalarea() = 0.0;

	for (SfcFacet* C = S.beginf(); C != S.endf(); ++C) {

		C->rtnmidpos(0) = ((S.beginv() + C->rtnsfcvertices(0))->rtnpos(0) + (S.beginv() + C->rtnsfcvertices(1))->rtnpos(0) + (S.beginv() + C->rtnsfcvertices(2))->rtnpos(0)) / 3.0f;
		C->rtnmidpos(1) = ((S.beginv() + C->rtnsfcvertices(0))->rtnpos(1) + (S.beginv() + C->rtnsfcvertices(1))->rtnpos(1) + (S.beginv() + C->rtnsfcvertices(2))->rtnpos(1)) / 3.0f;
		C->rtnmidpos(2) = ((S.beginv() + C->rtnsfcvertices(0))->rtnpos(2) + (S.beginv() + C->rtnsfcvertices(1))->rtnpos(2) + (S.beginv() + C->rtnsfcvertices(2))->rtnpos(2)) / 3.0f;

		C->rtnnormcalc(0, 0) = (S.beginv() + C->rtnsfcvertices(1))->rtnpos(0) - (S.beginv() + C->rtnsfcvertices(0))->rtnpos(0);
		C->rtnnormcalc(0, 1) = (S.beginv() + C->rtnsfcvertices(1))->rtnpos(1) - (S.beginv() + C->rtnsfcvertices(0))->rtnpos(1);
		C->rtnnormcalc(0, 2) = (S.beginv() + C->rtnsfcvertices(1))->rtnpos(2) - (S.beginv() + C->rtnsfcvertices(0))->rtnpos(2);

		C->rtnnormcalc(1, 0) = (S.beginv() + C->rtnsfcvertices(2))->rtnpos(0) - (S.beginv() + C->rtnsfcvertices(0))->rtnpos(0);
		C->rtnnormcalc(1, 1) = (S.beginv() + C->rtnsfcvertices(2))->rtnpos(1) - (S.beginv() + C->rtnsfcvertices(0))->rtnpos(1);
		C->rtnnormcalc(1, 2) = (S.beginv() + C->rtnsfcvertices(2))->rtnpos(2) - (S.beginv() + C->rtnsfcvertices(0))->rtnpos(2);

		C->rtnnormcalc(2, 0) = (S.beginv() + C->rtnsfcvertices(2))->rtnpos(0) - (S.beginv() + C->rtnsfcvertices(1))->rtnpos(0);
		C->rtnnormcalc(2, 1) = (S.beginv() + C->rtnsfcvertices(2))->rtnpos(1) - (S.beginv() + C->rtnsfcvertices(1))->rtnpos(1);
		C->rtnnormcalc(2, 2) = (S.beginv() + C->rtnsfcvertices(2))->rtnpos(2) - (S.beginv() + C->rtnsfcvertices(1))->rtnpos(2);

		C->rtnnormal(0) = (C->rtnnormcalc(0, 1) * C->rtnnormcalc(1, 2)) - (C->rtnnormcalc(0, 2) * C->rtnnormcalc(1, 1));
		C->rtnnormal(1) = (C->rtnnormcalc(0, 2) * C->rtnnormcalc(1, 0)) - (C->rtnnormcalc(0, 0) * C->rtnnormcalc(1, 2));
		C->rtnnormal(2) = (C->rtnnormcalc(0, 0) * C->rtnnormcalc(1, 1)) - (C->rtnnormcalc(0, 1) * C->rtnnormcalc(1, 0));

		n = sqrt((C->rtnnormal(0) * C->rtnnormal(0)) + (C->rtnnormal(1) * C->rtnnormal(1)) + (C->rtnnormal(2) * C->rtnnormal(2)));
		C->rtnarea() = n / 2.0f;

		C->rtnnormal(0) = C->rtnnormal(0) / n;
		C->rtnnormal(1) = C->rtnnormal(1) / n;
		C->rtnnormal(2) = C->rtnnormal(2) / n;

		S.rtntotalarea() = S.rtntotalarea() + (C->rtnarea() * C->rtnnormal(2));

		C->rtntrivector(0, 0) = -1.0f * C->rtnnormcalc(1, 0) / m;
		C->rtntrivector(0, 1) = -1.0f * C->rtnnormcalc(1, 1) / m;
		C->rtntrivector(0, 2) = -1.0f * C->rtnnormcalc(1, 2) / m;

		C->rtntrivector(1, 0) = C->rtnnormcalc(0, 0) / m;
		C->rtntrivector(1, 1) = C->rtnnormcalc(0, 1) / m;
		C->rtntrivector(1, 2) = C->rtnnormcalc(0, 2) / m;

		C->rtnnormcalc(1, 0) = -1.0f * C->rtnnormcalc(1, 0);
		C->rtnnormcalc(1, 1) = -1.0f * C->rtnnormcalc(1, 1);
		C->rtnnormcalc(1, 2) = -1.0f * C->rtnnormcalc(1, 2);

		C->rtnradvector(0) = (1.5 * C->rtnphotonmomvector(0)) - C->rtnnormal(0);
		C->rtnradvector(1) = (1.5 * C->rtnphotonmomvector(1)) - C->rtnnormal(1);
		C->rtnradvector(2) = (1.5 * C->rtnphotonmomvector(2)) - C->rtnnormal(2);

	}

	int step, i;
	float _t;
	AstFacet* C;

	omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(i,_t,C)
	for (step = 0; step < NTAU; step++) {
		for (C = A.beginf(); C != A.endf(); ++C) {
			_t = sqrt((C->rtnsfcillumvector(step, 0) * C->rtnsfcillumvector(step, 0)) + (C->rtnsfcillumvector(step, 1) * C->rtnsfcillumvector(step, 1)) + (C->rtnsfcillumvector(step, 2) * C->rtnsfcillumvector(step, 2)));
			for (i = 0; i != SFCNUMFACE; ++i) {
				(C->beginsf() + i)->rtnillumangle(step) = ((C->rtnsfcillumvector(step, 0) * (S.beginf() + i)->rtnnormal(0)) + (C->rtnsfcillumvector(step, 1) * (S.beginf() + i)->rtnnormal(1)) + (C->rtnsfcillumvector(step, 2) * (S.beginf() + i)->rtnnormal(2))) / _t;
				if (C->rtnshadow(step) == 0) { (C->beginsf() + i)->rtnshadow(step) = 0.0f; }
				if ((C->beginsf() + i)->rtnillumangle(step) <= 0.0f) { (C->beginsf() + i)->rtnshadow(step) = 0.0f; }
			}
		}
	}

};


// Determine asteroid global shape model geometry


void asteroidgeometry(Asteroid& A) {

	if (single_run == 1) {
		A.rtncarcentre(0) = AU * A.rtnpolcentre(0) * cos(A.rtnpolcentre(1)) * cos(A.rtnpolcentre(2));
		A.rtncarcentre(1) = AU * A.rtnpolcentre(0) * cos(A.rtnpolcentre(1)) * sin(A.rtnpolcentre(2));
		A.rtncarcentre(2) = AU * A.rtnpolcentre(0) * sin(A.rtnpolcentre(1));
	}
	else {
		A.rtnpolcentre(0) = sqrt((A.rtncarcentre(0) * A.rtncarcentre(0)) + (A.rtncarcentre(1) * A.rtncarcentre(1)) + (A.rtncarcentre(2) * A.rtncarcentre(2))) / AU;
	}

	A.rtnsolar() = (float)(1367.0 / (A.rtnpolcentre(0) * A.rtnpolcentre(0)));

	double tmatrix[3][3], gamma;

	gamma = acos(cos(A.rtnrotvector(0)) * sin(fabs(A.rtnrotvector(1))));

	tmatrix[0][0] = sin(A.rtnrotvector(0)) / sin(gamma);
	tmatrix[0][1] = -1.0 * cos(A.rtnrotvector(0)) * cos(A.rtnrotvector(0)) * sin(A.rtnrotvector(1)) * cos(A.rtnrotvector(1)) / sin(gamma);
	tmatrix[0][2] = cos(A.rtnrotvector(0)) * cos(A.rtnrotvector(1));
	tmatrix[1][0] = 0.0;
	tmatrix[1][1] = ((sin(A.rtnrotvector(0)) * sin(A.rtnrotvector(0))) + (cos(A.rtnrotvector(0)) * cos(A.rtnrotvector(0)) * cos(A.rtnrotvector(1)) * cos(A.rtnrotvector(1)))) / sin(gamma);
	tmatrix[1][2] = cos(A.rtnrotvector(0)) * sin(A.rtnrotvector(1));
	tmatrix[2][0] = -1.0 * cos(A.rtnrotvector(0)) * cos(A.rtnrotvector(1)) / sin(gamma);
	tmatrix[2][1] = -1.0 * sin(A.rtnrotvector(0)) * cos(A.rtnrotvector(0)) * sin(A.rtnrotvector(1)) / sin(gamma);
	tmatrix[2][2] = sin(A.rtnrotvector(0));

	A.rtnimatrix(0, 0) = sin(A.rtnrotvector(0)) / sin(gamma);
	A.rtnimatrix(0, 1) = 0.0;
	A.rtnimatrix(0, 2) = -1.0 * cos(A.rtnrotvector(0)) * cos(A.rtnrotvector(1)) / sin(gamma);
	A.rtnimatrix(1, 0) = -1.0 * cos(A.rtnrotvector(0)) * cos(A.rtnrotvector(0)) * sin(A.rtnrotvector(1)) * cos(A.rtnrotvector(1)) / sin(gamma);
	A.rtnimatrix(1, 1) = ((sin(A.rtnrotvector(0)) * sin(A.rtnrotvector(0))) + (cos(A.rtnrotvector(0)) * cos(A.rtnrotvector(0)) * cos(A.rtnrotvector(1)) * cos(A.rtnrotvector(1)))) / sin(gamma);
	A.rtnimatrix(1, 2) = -1.0 * sin(A.rtnrotvector(0)) * cos(A.rtnrotvector(0)) * sin(A.rtnrotvector(1)) / sin(gamma);
	A.rtnimatrix(2, 0) = cos(A.rtnrotvector(0)) * cos(A.rtnrotvector(1));
	A.rtnimatrix(2, 1) = cos(A.rtnrotvector(0)) * sin(A.rtnrotvector(1));
	A.rtnimatrix(2, 2) = sin(A.rtnrotvector(0));

	double rotangle, eclipos[3], rotphotonmomvector[3], ecliphotonmomvector[3];
	double normcalc[2][3], y[3], _n, _t, _x;
	double imatrix[3][3];
	int step;
	AstVertex* B;
	AstFacet* C;

	omp_set_num_threads(num_threads);
#pragma	omp parallel for schedule(dynamic) private(rotangle, eclipos, rotphotonmomvector, ecliphotonmomvector, normcalc, y, _n, _t, _x, imatrix, B, C)
	for (step = 0; step < NTAU; step++) {

		rotangle = 2.0 * PI * step / NTAU;

		for (B = A.beginv(); B != A.endv(); ++B) {

			B->rtnrotequpos(step, 0) = (B->rtnequpos(0) * cos(rotangle)) - (B->rtnequpos(1) * sin(rotangle));
			B->rtnrotequpos(step, 1) = (B->rtnequpos(0) * sin(rotangle)) + (B->rtnequpos(1) * cos(rotangle));
			B->rtnrotequpos(step, 2) = B->rtnequpos(2);

			eclipos[0] = (tmatrix[0][0] * B->rtnrotequpos(step, 0)) + (tmatrix[0][1] * B->rtnrotequpos(step, 1)) + (tmatrix[0][2] * B->rtnrotequpos(step, 2));
			eclipos[1] = (tmatrix[1][0] * B->rtnrotequpos(step, 0)) + (tmatrix[1][1] * B->rtnrotequpos(step, 1)) + (tmatrix[1][2] * B->rtnrotequpos(step, 2));
			eclipos[2] = (tmatrix[2][0] * B->rtnrotequpos(step, 0)) + (tmatrix[2][1] * B->rtnrotequpos(step, 1)) + (tmatrix[2][2] * B->rtnrotequpos(step, 2));

			B->rtnheliopos(step, 0) = A.rtncarcentre(0) + eclipos[0];
			B->rtnheliopos(step, 1) = A.rtncarcentre(1) + eclipos[1];
			B->rtnheliopos(step, 2) = A.rtncarcentre(2) + eclipos[2];

		}

		for (C = A.beginf(); C != A.endf(); ++C) {

			C->rtnmidpos(step, 0) = (((A.beginv() + C->rtnastvertices(0))->rtnrotequpos(step, 0)) + ((A.beginv() + C->rtnastvertices(1))->rtnrotequpos(step, 0)) + ((A.beginv() + C->rtnastvertices(2))->rtnrotequpos(step, 0))) / 3.0;
			C->rtnmidpos(step, 1) = (((A.beginv() + C->rtnastvertices(0))->rtnrotequpos(step, 1)) + ((A.beginv() + C->rtnastvertices(1))->rtnrotequpos(step, 1)) + ((A.beginv() + C->rtnastvertices(2))->rtnrotequpos(step, 1))) / 3.0;
			C->rtnmidpos(step, 2) = (((A.beginv() + C->rtnastvertices(0))->rtnrotequpos(step, 2)) + ((A.beginv() + C->rtnastvertices(1))->rtnrotequpos(step, 2)) + ((A.beginv() + C->rtnastvertices(2))->rtnrotequpos(step, 2))) / 3.0;

			C->rtnheliomidpos(step, 0) = (((A.beginv() + C->rtnastvertices(0))->rtnheliopos(step, 0)) + ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(step, 0)) + ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(step, 0))) / 3.0;
			C->rtnheliomidpos(step, 1) = (((A.beginv() + C->rtnastvertices(0))->rtnheliopos(step, 1)) + ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(step, 1)) + ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(step, 1))) / 3.0;
			C->rtnheliomidpos(step, 2) = (((A.beginv() + C->rtnastvertices(0))->rtnheliopos(step, 2)) + ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(step, 2)) + ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(step, 2))) / 3.0;

			_t = sqrt((C->rtnheliomidpos(step, 0) * C->rtnheliomidpos(step, 0)) + (C->rtnheliomidpos(step, 1) * C->rtnheliomidpos(step, 1)) + (C->rtnheliomidpos(step, 2) * C->rtnheliomidpos(step, 2)));

			normcalc[0][0] = ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(step, 0)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(step, 0));
			normcalc[0][1] = ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(step, 1)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(step, 1));
			normcalc[0][2] = ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(step, 2)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(step, 2));

			normcalc[1][0] = ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(step, 0)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(step, 0));
			normcalc[1][1] = ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(step, 1)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(step, 1));
			normcalc[1][2] = ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(step, 2)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(step, 2));

			C->rtnnormal(step, 0) = (normcalc[0][1] * normcalc[1][2]) - (normcalc[0][2] * normcalc[1][1]);
			C->rtnnormal(step, 1) = (normcalc[0][2] * normcalc[1][0]) - (normcalc[0][0] * normcalc[1][2]);
			C->rtnnormal(step, 2) = (normcalc[0][0] * normcalc[1][1]) - (normcalc[0][1] * normcalc[1][0]);

			_n = sqrt((C->rtnnormal(step, 0) * C->rtnnormal(step, 0)) + (C->rtnnormal(step, 1) * C->rtnnormal(step, 1)) + (C->rtnnormal(step, 2) * C->rtnnormal(step, 2)));
			C->rtnarea() = (float)(_n / 2.0);

			C->rtnnormal(step, 0) = C->rtnnormal(step, 0) / _n;
			C->rtnnormal(step, 1) = C->rtnnormal(step, 1) / _n;
			C->rtnnormal(step, 2) = C->rtnnormal(step, 2) / _n;

			_x = sqrt((normcalc[0][0] * normcalc[0][0]) + (normcalc[0][1] * normcalc[0][1]) + (normcalc[0][2] * normcalc[0][2]));

			normcalc[0][0] = normcalc[0][0] / _x;
			normcalc[0][1] = normcalc[0][1] / _x;
			normcalc[0][2] = normcalc[0][2] / _x;

			y[0] = (C->rtnnormal(step, 1) * normcalc[0][2]) - (C->rtnnormal(step, 2) * normcalc[0][1]);
			y[1] = (C->rtnnormal(step, 2) * normcalc[0][0]) - (C->rtnnormal(step, 0) * normcalc[0][2]);
			y[2] = (C->rtnnormal(step, 0) * normcalc[0][1]) - (C->rtnnormal(step, 1) * normcalc[0][0]);

			imatrix[0][0] = normcalc[0][0];
			imatrix[0][1] = normcalc[0][1];
			imatrix[0][2] = normcalc[0][2];
			imatrix[1][0] = y[0];
			imatrix[1][1] = y[1];
			imatrix[1][2] = y[2];
			imatrix[2][0] = C->rtnnormal(step, 0);
			imatrix[2][1] = C->rtnnormal(step, 1);
			imatrix[2][2] = C->rtnnormal(step, 2);

			C->rtnimatrix(step, 0, 0) = normcalc[0][0];
			C->rtnimatrix(step, 0, 1) = y[0];
			C->rtnimatrix(step, 0, 2) = C->rtnnormal(step, 0);
			C->rtnimatrix(step, 1, 0) = normcalc[0][1];
			C->rtnimatrix(step, 1, 1) = y[1];
			C->rtnimatrix(step, 1, 2) = C->rtnnormal(step, 1);
			C->rtnimatrix(step, 2, 0) = normcalc[0][2];
			C->rtnimatrix(step, 2, 1) = y[2];
			C->rtnimatrix(step, 2, 2) = C->rtnnormal(step, 2);

			C->rtnillumvector(step, 0) = -1.0 * ((imatrix[0][0] * C->rtnheliomidpos(step, 0)) + (imatrix[0][1] * C->rtnheliomidpos(step, 1)) + (imatrix[0][2] * C->rtnheliomidpos(step, 2)));
			C->rtnillumvector(step, 1) = -1.0 * ((imatrix[1][0] * C->rtnheliomidpos(step, 0)) + (imatrix[1][1] * C->rtnheliomidpos(step, 1)) + (imatrix[1][2] * C->rtnheliomidpos(step, 2)));
			C->rtnillumvector(step, 2) = -1.0 * ((imatrix[2][0] * C->rtnheliomidpos(step, 0)) + (imatrix[2][1] * C->rtnheliomidpos(step, 1)) + (imatrix[2][2] * C->rtnheliomidpos(step, 2)));

			C->rtnillumangle(step) = acos(C->rtnillumvector(step, 2) / _t);
			if (C->rtnillumangle(step) > (PI / 2.0)) { C->rtnshadow(step) = 0; }

			C->rtnsfcillumvector(step, 0) = (float)(10000.0 * C->rtnillumvector(step, 0) / _t);
			C->rtnsfcillumvector(step, 1) = (float)(10000.0 * C->rtnillumvector(step, 1) / _t);
			C->rtnsfcillumvector(step, 2) = (float)(10000.0 * C->rtnillumvector(step, 2) / _t);

			rotphotonmomvector[0] = (C->rtnphotonmomvector(0) * cos(rotangle)) - (C->rtnphotonmomvector(1) * sin(rotangle));
			rotphotonmomvector[1] = (C->rtnphotonmomvector(0) * sin(rotangle)) + (C->rtnphotonmomvector(1) * cos(rotangle));
			rotphotonmomvector[2] = C->rtnphotonmomvector(2);

			ecliphotonmomvector[0] = (tmatrix[0][0] * rotphotonmomvector[0]) + (tmatrix[0][1] * rotphotonmomvector[1]) + (tmatrix[0][2] * rotphotonmomvector[2]);
			ecliphotonmomvector[1] = (tmatrix[1][0] * rotphotonmomvector[0]) + (tmatrix[1][1] * rotphotonmomvector[1]) + (tmatrix[1][2] * rotphotonmomvector[2]);
			ecliphotonmomvector[2] = (tmatrix[2][0] * rotphotonmomvector[0]) + (tmatrix[2][1] * rotphotonmomvector[1]) + (tmatrix[2][2] * rotphotonmomvector[2]);

			C->rtnradvector(step, 0) = (1.5 * ecliphotonmomvector[0]) - C->rtnnormal(step, 0);
			C->rtnradvector(step, 1) = (1.5 * ecliphotonmomvector[1]) - C->rtnnormal(step, 1);
			C->rtnradvector(step, 2) = (1.5 * ecliphotonmomvector[2]) - C->rtnnormal(step, 2);

		}

	}

};


// Code to calculate surface roughness global selfheating viewfactors


void calculate_surface_global_viewfactors(Surface& S, Asteroid& A) {

	int i;
	double normcalc[2][3], normal[3], n, x, y[3];
	AstFacet* C;

	omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(normcalc,normal,n,x,y,C)
	for (i = 0; i < ASTNUMFACE; i++) {

		C = A.beginf() + i;

		normcalc[0][0] = ((A.beginv() + C->rtnastvertices(1))->rtnequpos(0)) - ((A.beginv() + C->rtnastvertices(0))->rtnequpos(0));
		normcalc[0][1] = ((A.beginv() + C->rtnastvertices(1))->rtnequpos(1)) - ((A.beginv() + C->rtnastvertices(0))->rtnequpos(1));
		normcalc[0][2] = ((A.beginv() + C->rtnastvertices(1))->rtnequpos(2)) - ((A.beginv() + C->rtnastvertices(0))->rtnequpos(2));

		normcalc[1][0] = ((A.beginv() + C->rtnastvertices(2))->rtnequpos(0)) - ((A.beginv() + C->rtnastvertices(0))->rtnequpos(0));
		normcalc[1][1] = ((A.beginv() + C->rtnastvertices(2))->rtnequpos(1)) - ((A.beginv() + C->rtnastvertices(0))->rtnequpos(1));
		normcalc[1][2] = ((A.beginv() + C->rtnastvertices(2))->rtnequpos(2)) - ((A.beginv() + C->rtnastvertices(0))->rtnequpos(2));

		normal[0] = (normcalc[0][1] * normcalc[1][2]) - (normcalc[0][2] * normcalc[1][1]);
		normal[1] = (normcalc[0][2] * normcalc[1][0]) - (normcalc[0][0] * normcalc[1][2]);
		normal[2] = (normcalc[0][0] * normcalc[1][1]) - (normcalc[0][1] * normcalc[1][0]);

		n = sqrt((normal[0] * normal[0]) + (normal[1] * normal[1]) + (normal[2] * normal[2]));

		normal[0] = normal[0] / n;
		normal[1] = normal[1] / n;
		normal[2] = normal[2] / n;

		x = sqrt((normcalc[0][0] * normcalc[0][0]) + (normcalc[0][1] * normcalc[0][1]) + (normcalc[0][2] * normcalc[0][2]));

		normcalc[0][0] = normcalc[0][0] / x;
		normcalc[0][1] = normcalc[0][1] / x;
		normcalc[0][2] = normcalc[0][2] / x;

		y[0] = (normal[1] * normcalc[0][2]) - (normal[2] * normcalc[0][1]);
		y[1] = (normal[2] * normcalc[0][0]) - (normal[0] * normcalc[0][2]);
		y[2] = (normal[0] * normcalc[0][1]) - (normal[1] * normcalc[0][0]);

		C->rtnimatrix0(0, 0) = normcalc[0][0];
		C->rtnimatrix0(0, 1) = normcalc[0][1];
		C->rtnimatrix0(0, 2) = normcalc[0][2];
		C->rtnimatrix0(1, 0) = y[0];
		C->rtnimatrix0(1, 1) = y[1];
		C->rtnimatrix0(1, 2) = y[2];
		C->rtnimatrix0(2, 0) = normal[0];
		C->rtnimatrix0(2, 1) = normal[1];
		C->rtnimatrix0(2, 2) = normal[2];

	}

	float trivector_multiplier[m * m][2];
	int k, l, p;

	for (k = 1; k != m + 1; ++k) {
		for (l = 1; l != k + 1; ++l) {
			p = ((k - 1) * (k - 1)) + (2 * (l - 1));
			trivector_multiplier[p][0] = (float)k - (1.0f / 3.0f);
			trivector_multiplier[p][1] = (float)l - (2.0f / 3.0f);
		}
	}

	for (k = 2; k != m + 1; ++k) {
		for (l = 1; l != k; ++l) {
			p = ((k - 1) * (k - 1)) + (2 * (l - 1)) + 1;
			trivector_multiplier[p][0] = (float)k - (2.0f / 3.0f);
			trivector_multiplier[p][1] = (float)l - (1.0f / 3.0f);
		}
	}

	int j;
	float viewfactor_vector[3], visibility, angle;
	float midpos[3];
	float trans_co[2];
	float a;
	float shadow_radius;
	float crater_radius = 10.0f;
	SfcFacet* SF;
	DepAstFacet* DAF;

	omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(j,p,viewfactor_vector,visibility,angle,midpos,trans_co,a,shadow_radius,C,SF,DAF)
	for (i = 0; i < ASTNUMFACE; i++) {

		C = A.beginf() + i;

		if (C->rtnastdependentfacets() != 0) {
			for (DAF = C->begindaf(); DAF != C->enddaf(); ++DAF) {

				viewfactor_vector[0] = (float)((C->rtnimatrix0(0, 0) * DAF->rtnviewfactor_vector(0)) + ((C->rtnimatrix0(0, 1) * DAF->rtnviewfactor_vector(1))) + (C->rtnimatrix0(0, 2) * DAF->rtnviewfactor_vector(2)));
				viewfactor_vector[1] = (float)((C->rtnimatrix0(1, 0) * DAF->rtnviewfactor_vector(0)) + ((C->rtnimatrix0(1, 1) * DAF->rtnviewfactor_vector(1))) + (C->rtnimatrix0(1, 2) * DAF->rtnviewfactor_vector(2)));
				viewfactor_vector[2] = (float)((C->rtnimatrix0(2, 0) * DAF->rtnviewfactor_vector(0)) + ((C->rtnimatrix0(2, 1) * DAF->rtnviewfactor_vector(1))) + (C->rtnimatrix0(2, 2) * DAF->rtnviewfactor_vector(2)));
				a = fabs(viewfactor_vector[2]) / sqrt((viewfactor_vector[0] * viewfactor_vector[0]) + (viewfactor_vector[1] * viewfactor_vector[1]));

				for (j = 0; j != SFCNUMFACE; ++j) {

					SF = S.beginf() + j;
					visibility = 1.0f;

					angle = (viewfactor_vector[0] * SF->rtnnormal(0)) + (viewfactor_vector[1] * SF->rtnnormal(1)) + (viewfactor_vector[2] * SF->rtnnormal(2));
					if (angle < 0.0f) { visibility = 0.0f; }

					for (p = 0; p != (m * m); ++p) {

						midpos[0] = (S.beginv() + SF->rtnsfcvertices(2))->rtnpos(0) + (trivector_multiplier[p][0] * SF->rtntrivector(0, 0)) + (trivector_multiplier[p][1] * SF->rtntrivector(1, 0));
						midpos[1] = (S.beginv() + SF->rtnsfcvertices(2))->rtnpos(1) + (trivector_multiplier[p][0] * SF->rtntrivector(0, 1)) + (trivector_multiplier[p][1] * SF->rtntrivector(1, 1));
						midpos[2] = (S.beginv() + SF->rtnsfcvertices(2))->rtnpos(2) + (trivector_multiplier[p][0] * SF->rtntrivector(0, 2)) + (trivector_multiplier[p][1] * SF->rtntrivector(1, 2));

						trans_co[1] = (-1.0f * viewfactor_vector[1] * midpos[0]) + (viewfactor_vector[0] * midpos[1]);
						shadow_radius = sqrt((crater_radius * crater_radius) - (trans_co[1] * trans_co[1]));
						trans_co[0] = (viewfactor_vector[0] * midpos[0]) + (viewfactor_vector[1] * midpos[1]) - shadow_radius;

						if (midpos[2] <= (a * trans_co[0]) && visibility != 0.0f) { visibility = visibility - (1.0f / (float)(m * m)); }
						if (visibility < 0.0f) { visibility = 0.0f; }

					}

					DAF->rtnviewfactor(j + 1) = (DAF->rtnviewfactor(0) * visibility * angle) / viewfactor_vector[2];
					if (DAF->rtnviewfactor(j + 1) > 1.0f) { DAF->rtnviewfactor(j + 1) = 1.0f; }

				}

			}
		}
	}

};


// Code to test for asteroid shadows


int astshadowtest(Asteroid& A, AstFacet* G, AstFacet* H, int step) {

	double t, tcalc1, tcalc2;
	double xpoint[3];
	double normcalc1[3], normcalc2[3], normcalc3[3];
	double xvector1[3], xvector2[3], xvector3[3];
	double test1[3], test2[3], test3[3];
	double condition1, condition2, condition3;

	tcalc2 = (G->rtnheliomidpos(step, 0) * H->rtnnormal(step, 0)) + (G->rtnheliomidpos(step, 1) * H->rtnnormal(step, 1)) + (G->rtnheliomidpos(step, 2) * H->rtnnormal(step, 2));
	tcalc1 = (H->rtnheliomidpos(step, 0) * H->rtnnormal(step, 0)) + (H->rtnheliomidpos(step, 1) * H->rtnnormal(step, 1)) + (H->rtnheliomidpos(step, 2) * H->rtnnormal(step, 2));
	t = tcalc1 / tcalc2;

	xpoint[0] = t * G->rtnheliomidpos(step, 0);
	xpoint[1] = t * G->rtnheliomidpos(step, 1);
	xpoint[2] = t * G->rtnheliomidpos(step, 2);

	normcalc1[0] = ((A.beginv() + H->rtnastvertices(1))->rtnheliopos(step, 0)) - ((A.beginv() + H->rtnastvertices(0))->rtnheliopos(step, 0));
	normcalc1[1] = ((A.beginv() + H->rtnastvertices(1))->rtnheliopos(step, 1)) - ((A.beginv() + H->rtnastvertices(0))->rtnheliopos(step, 1));
	normcalc1[2] = ((A.beginv() + H->rtnastvertices(1))->rtnheliopos(step, 2)) - ((A.beginv() + H->rtnastvertices(0))->rtnheliopos(step, 2));

	normcalc2[0] = ((A.beginv() + H->rtnastvertices(0))->rtnheliopos(step, 0)) - ((A.beginv() + H->rtnastvertices(2))->rtnheliopos(step, 0));
	normcalc2[1] = ((A.beginv() + H->rtnastvertices(0))->rtnheliopos(step, 1)) - ((A.beginv() + H->rtnastvertices(2))->rtnheliopos(step, 1));
	normcalc2[2] = ((A.beginv() + H->rtnastvertices(0))->rtnheliopos(step, 2)) - ((A.beginv() + H->rtnastvertices(2))->rtnheliopos(step, 2));

	normcalc3[0] = ((A.beginv() + H->rtnastvertices(2))->rtnheliopos(step, 0)) - ((A.beginv() + H->rtnastvertices(1))->rtnheliopos(step, 0));
	normcalc3[1] = ((A.beginv() + H->rtnastvertices(2))->rtnheliopos(step, 1)) - ((A.beginv() + H->rtnastvertices(1))->rtnheliopos(step, 1));
	normcalc3[2] = ((A.beginv() + H->rtnastvertices(2))->rtnheliopos(step, 2)) - ((A.beginv() + H->rtnastvertices(1))->rtnheliopos(step, 2));

	xvector1[0] = xpoint[0] - ((A.beginv() + H->rtnastvertices(0))->rtnheliopos(step, 0));
	xvector1[1] = xpoint[1] - ((A.beginv() + H->rtnastvertices(0))->rtnheliopos(step, 1));
	xvector1[2] = xpoint[2] - ((A.beginv() + H->rtnastvertices(0))->rtnheliopos(step, 2));

	xvector3[0] = xpoint[0] - ((A.beginv() + H->rtnastvertices(1))->rtnheliopos(step, 0));
	xvector3[1] = xpoint[1] - ((A.beginv() + H->rtnastvertices(1))->rtnheliopos(step, 1));
	xvector3[2] = xpoint[2] - ((A.beginv() + H->rtnastvertices(1))->rtnheliopos(step, 2));

	xvector2[0] = xpoint[0] - ((A.beginv() + H->rtnastvertices(2))->rtnheliopos(step, 0));
	xvector2[1] = xpoint[1] - ((A.beginv() + H->rtnastvertices(2))->rtnheliopos(step, 1));
	xvector2[2] = xpoint[2] - ((A.beginv() + H->rtnastvertices(2))->rtnheliopos(step, 2));

	test1[0] = (normcalc1[1] * xvector1[2]) - (normcalc1[2] * xvector1[1]);
	test1[1] = (normcalc1[2] * xvector1[0]) - (normcalc1[0] * xvector1[2]);
	test1[2] = (normcalc1[0] * xvector1[1]) - (normcalc1[1] * xvector1[0]);

	test2[0] = (normcalc2[1] * xvector2[2]) - (normcalc2[2] * xvector2[1]);
	test2[1] = (normcalc2[2] * xvector2[0]) - (normcalc2[0] * xvector2[2]);
	test2[2] = (normcalc2[0] * xvector2[1]) - (normcalc2[1] * xvector2[0]);

	test3[0] = (normcalc3[1] * xvector3[2]) - (normcalc3[2] * xvector3[1]);
	test3[1] = (normcalc3[2] * xvector3[0]) - (normcalc3[0] * xvector3[2]);
	test3[2] = (normcalc3[0] * xvector3[1]) - (normcalc3[1] * xvector3[0]);

	condition1 = (test1[0] * H->rtnnormal(step, 0)) + (test1[1] * H->rtnnormal(step, 1)) + (test1[2] * H->rtnnormal(step, 2));
	condition2 = (test2[0] * H->rtnnormal(step, 0)) + (test2[1] * H->rtnnormal(step, 1)) + (test2[2] * H->rtnnormal(step, 2));
	condition3 = (test3[0] * H->rtnnormal(step, 0)) + (test3[1] * H->rtnnormal(step, 1)) + (test3[2] * H->rtnnormal(step, 2));

	if (condition1 >= 0.0 && condition2 >= 0.0 && condition3 >= 0.0) { return 0; }
	else { return 1; }

};


// Code to determine asteroid shadowing


//void astshadow(Asteroid &A, int step) {

//	double g, h;

//	for (AstFacet *G=A.beginf();G!=A.endf();++G) {

//		g=sqrt((G->rtnheliomidpos(step,0)*G->rtnheliomidpos(step,0))+(G->rtnheliomidpos(step,1)*G->rtnheliomidpos(step,1))+(G->rtnheliomidpos(step,2)*G->rtnheliomidpos(step,2)));

//		for (AstFacet *H=(G+1);H!=A.endf();++H) {

//			h=sqrt((H->rtnheliomidpos(step,0)*H->rtnheliomidpos(step,0))+(H->rtnheliomidpos(step,1)*H->rtnheliomidpos(step,1))+(H->rtnheliomidpos(step,2)*H->rtnheliomidpos(step,2)));

//			if (g>h) {

//				if (G->rtnshadow(step)==1) {
//					if ( astshadowtest(A,G,H,step)==0 ) { G->rtnshadow(step)=0; }
//				}

//			}	else {

//				if (H->rtnshadow(step)==1) {
//					if ( astshadowtest(A,H,G,step)==0 ) { H->rtnshadow(step)=0; }
//				}

//			}

//		}

//	}

//};


void astshadow(Asteroid& A, int step) {

	AstFacet* H;
	DepAstFacet* I;
	double g, h;

	for (AstFacet* G = A.beginf(); G != A.endf(); ++G) {

		if (G->rtnastdependentfacets() != 0 && G->rtnshadow(step) == 1) {

			g = sqrt((G->rtnheliomidpos(step, 0) * G->rtnheliomidpos(step, 0)) + (G->rtnheliomidpos(step, 1) * G->rtnheliomidpos(step, 1)) + (G->rtnheliomidpos(step, 2) * G->rtnheliomidpos(step, 2)));

			for (I = G->begindaf(); I != G->enddaf(); ++I) {

				H = A.beginf() + I->rtndepastfacetnum();
				h = sqrt((H->rtnheliomidpos(step, 0) * H->rtnheliomidpos(step, 0)) + (H->rtnheliomidpos(step, 1) * H->rtnheliomidpos(step, 1)) + (H->rtnheliomidpos(step, 2) * H->rtnheliomidpos(step, 2)));

				if (g > h) {

					if (astshadowtest(A, G, H, step) == 0) { G->rtnshadow(step) = 0; }

				}

			}

		}

	}

};


// Determine facet illumination from direct and single scattered sunlight (not including thermal radiation!)


void singlescattering(Asteroid& A) {

	AstFacet* C;
	DepAstFacet* G;
	int step;

	omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(C, G)
	for (step = 0; step < NTAU; step++) {

		for (C = A.beginf(); C != A.endf(); ++C) {

			C->rtntotalillumination(step) = A.rtnsolar() * (float)cos(C->rtnillumangle(step)) * C->rtnshadow(step);

			if (C->rtnastdependentfacets() != 0) {
				for (G = C->begindaf(); G != C->enddaf(); ++G) {
					C->rtntotalillumination(step) = C->rtntotalillumination(step) + (G->rtnviewfactor(0) * A.rtnsolar() * (A.beginf() + G->rtndepastfacetnum())->rtnalbedo() * (float)cos((A.beginf() + G->rtndepastfacetnum())->rtnillumangle(step)) * (A.beginf() + G->rtndepastfacetnum())->rtnshadow(step));
				}
			}

		}

	}

};


// Determine facet illumination from direct and multiple scattered sunlight (not including thermal radiation!)


void multiplescattering(Asteroid& A) {

	int step;
	float poweraccuracy(scattering_accuracy);

	AstFacet* C;
	DepAstFacet* G;

	bool iterate;
	int facetnumber, iteratetest;
	float beforesum, aftersum;
	float** solarpower, ** powerold, ** powernew;
	solarpower = new float* [num_threads];
	powerold = new float* [num_threads];
	powernew = new float* [num_threads];
	for (int i = 0; i != num_threads; ++i) {
		solarpower[i] = new float[ASTNUMFACE];
		powerold[i] = new float[ASTNUMFACE];
		powernew[i] = new float[ASTNUMFACE];
	}

	omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(iterate, facetnumber, iteratetest, beforesum, aftersum, C, G)
	for (step = 0; step < NTAU; ++step) {

		for (facetnumber = 0; facetnumber != ASTNUMFACE; ++facetnumber) {
			C = A.beginf() + facetnumber;
			solarpower[omp_get_thread_num()][facetnumber] = A.rtnsolar() * (float)cos(C->rtnillumangle(step)) * C->rtnshadow(step);
			powerold[omp_get_thread_num()][facetnumber] = C->rtnalbedo() * C->rtntotalillumination(step);
		}

		iterate = true;

		while (iterate == true) {

			for (facetnumber = 0; facetnumber != ASTNUMFACE; ++facetnumber) {

				C = A.beginf() + facetnumber;

				if (C->rtnastdependentfacets() != 0) {

					beforesum = 0.0f;
					aftersum = 0.0f;

					for (G = C->begindaf(); G != C->enddaf(); ++G) {
						if ((A.beginf() + G->rtndepastfacetnum()) > C) {
							aftersum = aftersum + (G->rtnviewfactor(0) * powerold[omp_get_thread_num()][G->rtndepastfacetnum()]);
						}
						else {
							beforesum = beforesum + (G->rtnviewfactor(0) * powernew[omp_get_thread_num()][G->rtndepastfacetnum()]);
						}
					}

					powernew[omp_get_thread_num()][facetnumber] = C->rtnalbedo() * (solarpower[omp_get_thread_num()][facetnumber] + aftersum + beforesum);

				}
				else {

					powernew[omp_get_thread_num()][facetnumber] = powerold[omp_get_thread_num()][facetnumber];

				}

			}

			iteratetest = 0;

			for (facetnumber = 0; facetnumber != ASTNUMFACE; ++facetnumber) {
				if (fabs(powernew[omp_get_thread_num()][facetnumber] - powerold[omp_get_thread_num()][facetnumber]) > poweraccuracy) { iteratetest = iteratetest + 1; }
				powerold[omp_get_thread_num()][facetnumber] = powernew[omp_get_thread_num()][facetnumber];
			}

			if (iteratetest == 0) { iterate = false; }

		}

		for (facetnumber = 0; facetnumber != ASTNUMFACE; ++facetnumber) {
			C = A.beginf() + facetnumber;
			C->rtntotalillumination(step) = powernew[omp_get_thread_num()][facetnumber] / C->rtnalbedo();
		}

	}

};


// Initialise temperatures to begin thermal modelling


void initialisetemperatures(Asteroid& A) {

	int facetnumber, step, i;
	float temp;
	AstFacet* C;

	omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(C, temp, step, i)
	for (facetnumber = 0; facetnumber < ASTNUMFACE; facetnumber++) {

		C = A.beginf() + facetnumber;
		temp = 0.0f;

		for (step = 0; step != NTAU; ++step) {
			temp = temp + ((1.0f / NTAU) * C->rtntotalillumination(step));
		}
		temp = pow((1.0f - C->rtnalbedo()) * temp / (EPS * SIGMA), 0.25f);

		for (i = 0; i != (NX + 1); ++i) {
			C->rtntemp(i) = temp;
		}
		for (step = 0; step != NTAU; ++step) {
			C->rtnstemp(0, step) = temp;
			C->rtnstemp(1, step) = temp;
		}

	}

};


// Finite difference method for constant thermal inertia


void inertiadifference(Asteroid& A) {

	float DZ, DTAU;
	float SBC, SBCD;

	float DT, TR[2], temp[2];
	float incidentflux;

	DZ = 8.0f / NX;
	DTAU = 2.0f * (float)PI / NTAU;

	AstFacet* C;
	DepAstFacet* G;
	int step, facetnumber, i;

	for (step = 0; step != NTAU; ++step) {

		omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(C)
		for (facetnumber = 0; facetnumber < ASTNUMFACE; facetnumber++) {
			C = A.beginf() + facetnumber;
			C->rtnoperated() = false;
		}

		omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(C, G, i, SBC, SBCD, DT, TR, temp, incidentflux)
		for (facetnumber = 0; facetnumber < ASTNUMFACE; facetnumber++) {

			C = A.beginf() + facetnumber;

			if (C->rtnastconverged() == false) {

				C->rtnstemp(0, step) = C->rtnstemp(1, step);
				C->rtnoperated() = true;

				// Finite difference scheme

				temp[0] = C->rtntemp(0);
				for (i = 1; i != NX; ++i) {
					temp[1] = C->rtntemp(i);
					C->rtntemp(i) = C->rtntemp(i) + ((DTAU / (DZ * DZ)) * (C->rtntemp(i + 1) - (2.0f * C->rtntemp(i)) + temp[0]));
					temp[0] = temp[1];
				}

				// Incident flux

				incidentflux = 0.0f;
				if (C->rtnastdependentfacets() != 0) {
					for (G = C->begindaf(); G != C->enddaf(); ++G) {
						if ((A.beginf() + G->rtndepastfacetnum())->rtnoperated() == true) {
							incidentflux = incidentflux + (EPS * SIGMA * pow((A.beginf() + G->rtndepastfacetnum())->rtnstemp(0, step), 4.0f) * G->rtnviewfactor(0));
						}
						else {
							incidentflux = incidentflux + (EPS * SIGMA * pow((A.beginf() + G->rtndepastfacetnum())->rtnstemp(1, step), 4.0f) * G->rtnviewfactor(0));
						}
					}
				}

				// Apply surface boundary condition

				DT = 1.0f;
				TR[0] = C->rtntemp(0);

				while (DT >= TACC / 10.0) {
					SBC = ((1.0f - C->rtnalbedo()) * C->rtntotalillumination(step)) + incidentflux + ((C->rtninert() * sqrt(2.0f * (float)PI / A.rtnrotperiod()) / DZ) * (C->rtntemp(1) - TR[0])) - (EPS * SIGMA * pow(TR[0], 4.0f));
					SBCD = (C->rtninert() * sqrt(2.0f * (float)PI / A.rtnrotperiod()) / DZ) + (4.0f * EPS * SIGMA * pow(TR[0], 3.0f));
					TR[1] = TR[0] + (SBC / SBCD);
					DT = fabs(TR[1] - TR[0]);
					TR[0] = TR[1];
				}

				C->rtntemp(0) = TR[1];
				C->rtnstemp(1, step) = TR[1];

				// Apply internal boundary condition

				C->rtntemp(NX) = C->rtntemp(NX - 1);

			}

		}

	}

};


// Instantaneous equilibrium thermal model


void instequilibrium(Asteroid& A) {

	AstFacet* C;
	DepAstFacet* G;
	float incidentflux;
	int step, facetnumber;

	for (step = 0; step != NTAU; ++step) {

		omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(C)
		for (facetnumber = 0; facetnumber < ASTNUMFACE; facetnumber++) {
			C = A.beginf() + facetnumber;
			C->rtnoperated() = false;
		}

		omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(C, G, incidentflux)
		for (facetnumber = 0; facetnumber < ASTNUMFACE; facetnumber++) {

			C = A.beginf() + facetnumber;

			C->rtnstemp(0, step) = C->rtnstemp(1, step);
			C->rtnoperated() = true;

			// Incident flux

			incidentflux = 0.0f;

			if (C->rtnastdependentfacets() != 0) {
				for (G = C->begindaf(); G != C->enddaf(); ++G) {
					if ((A.beginf() + G->rtndepastfacetnum())->rtnoperated() == true) {
						incidentflux = incidentflux + (EPS * SIGMA * pow((A.beginf() + G->rtndepastfacetnum())->rtnstemp(0, step), 4.0f) * G->rtnviewfactor(0));
					}
					else {
						incidentflux = incidentflux + (EPS * SIGMA * pow((A.beginf() + G->rtndepastfacetnum())->rtnstemp(1, step), 4.0f) * G->rtnviewfactor(0));
					}
				}

			}

			// Instantaneous radiative equilibrium

			C->rtnstemp(1, step) = pow(((((1.0f - C->rtnalbedo()) * C->rtntotalillumination(step)) + incidentflux) / (EPS * SIGMA)), 0.25f);

		}

	}

};


// Calculate resolved Yarkovsky forces and YORP torques


void resolved_yark_yorp(Asteroid& A, Surface& S) {

	int i, j, step;
	AstFacet* C;
	SubAstFacet* D;
	SfcFacet* E;
	DepAstFacet* DAF;

	double _t, viewfactor_vector[3];

	double subfacetsolarpressure, subfacetreflectedpressure, subfacetthermalpressure;
	double roughfacetsolarpressure[3], roughfacetreflectedpressure[3], roughfacetthermalpressure[3];
	double helioroughfacetsolarpressure[3], helioroughfacetreflectedpressure[3], helioroughfacetthermalpressure[3];
	double equroughfacetsolarpressure[3], equroughfacetreflectedpressure[3], equroughfacetthermalpressure[3];
	double roughfacetsolartorq[3], roughfacetreflectedtorq[3], roughfacetthermaltorq[3];

	double astfacetsolarpressure, astfacetreflectedpressure, astfacetthermalpressure;
	double smoothfacetsolarpressure[3], heliosmoothfacetsolarpressure[3], heliosmoothfacetreflectedpressure[3], heliosmoothfacetthermalpressure[3];
	double equsmoothfacetsolarpressure[3], equsmoothfacetreflectedpressure[3], equsmoothfacetthermalpressure[3];
	double smoothfacetsolartorq[3], smoothfacetreflectedtorq[3], smoothfacetthermaltorq[3];

	omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(j,step,C,D,E,DAF,_t,viewfactor_vector,subfacetsolarpressure,subfacetreflectedpressure,subfacetthermalpressure,roughfacetsolarpressure,roughfacetreflectedpressure,roughfacetthermalpressure,helioroughfacetsolarpressure,helioroughfacetreflectedpressure,helioroughfacetthermalpressure,equroughfacetsolarpressure,equroughfacetreflectedpressure,equroughfacetthermalpressure,roughfacetsolartorq,roughfacetreflectedtorq,roughfacetthermaltorq,astfacetsolarpressure,astfacetreflectedpressure,astfacetthermalpressure,smoothfacetsolarpressure,heliosmoothfacetsolarpressure,heliosmoothfacetreflectedpressure,heliosmoothfacetthermalpressure,equsmoothfacetsolarpressure,equsmoothfacetreflectedpressure,equsmoothfacetthermalpressure,smoothfacetsolartorq,smoothfacetreflectedtorq,smoothfacetthermaltorq)
	for (i = 0; i < ASTNUMFACE; i++) {

		C = A.beginf() + i;

		for (step = 0; step < NTAU; step++) {

			roughfacetsolarpressure[0] = 0.0, roughfacetsolarpressure[1] = 0.0, roughfacetsolarpressure[2] = 0.0;
			roughfacetreflectedpressure[0] = 0.0, roughfacetreflectedpressure[1] = 0.0, roughfacetreflectedpressure[2] = 0.0;
			roughfacetthermalpressure[0] = 0.0, roughfacetthermalpressure[1] = 0.0, roughfacetthermalpressure[2] = 0.0;

			_t = sqrt((C->rtnillumvector(step, 0) * C->rtnillumvector(step, 0)) + (C->rtnillumvector(step, 1) * C->rtnillumvector(step, 1)) + (C->rtnillumvector(step, 2) * C->rtnillumvector(step, 2)));

			if (roughness_thermal_model == 1) {

				// Roughness facets

				for (j = 0; j != SFCNUMFACE; ++j) {

					D = C->beginsf() + j;
					E = S.beginf() + j;

					subfacetsolarpressure = A.rtnsolar() * D->rtnillumangle(step) * D->rtnshadow(step) * E->rtnarea() / cvel;
					roughfacetsolarpressure[0] = roughfacetsolarpressure[0] - (subfacetsolarpressure * C->rtnillumvector(step, 0) / _t);
					roughfacetsolarpressure[1] = roughfacetsolarpressure[1] - (subfacetsolarpressure * C->rtnillumvector(step, 1) / _t);
					roughfacetsolarpressure[2] = roughfacetsolarpressure[2] - (subfacetsolarpressure * C->rtnillumvector(step, 2) / _t);

					subfacetreflectedpressure = (2.0 / 3.0) * C->rtnalbedo() * D->rtntotalillumination(step) * E->rtnarea() / cvel;
					roughfacetreflectedpressure[0] = roughfacetreflectedpressure[0] + (subfacetreflectedpressure * E->rtnradvector(0));
					roughfacetreflectedpressure[1] = roughfacetreflectedpressure[1] + (subfacetreflectedpressure * E->rtnradvector(1));
					roughfacetreflectedpressure[2] = roughfacetreflectedpressure[2] + (subfacetreflectedpressure * E->rtnradvector(2));

					subfacetthermalpressure = (2.0 / 3.0) * EPS * SIGMA * pow((double)D->rtnstemp(1, step), 4.0) * E->rtnarea() / cvel;
					roughfacetthermalpressure[0] = roughfacetthermalpressure[0] + (subfacetthermalpressure * E->rtnradvector(0));
					roughfacetthermalpressure[1] = roughfacetthermalpressure[1] + (subfacetthermalpressure * E->rtnradvector(1));
					roughfacetthermalpressure[2] = roughfacetthermalpressure[2] + (subfacetthermalpressure * E->rtnradvector(2));

					if (C->rtnastdependentfacets() != 0) {
						for (DAF = C->begindaf(); DAF != C->enddaf(); ++DAF) {

							viewfactor_vector[0] = (C->rtnimatrix0(0, 0) * DAF->rtnviewfactor_vector(0)) + ((C->rtnimatrix0(0, 1) * DAF->rtnviewfactor_vector(1))) + (C->rtnimatrix0(0, 2) * DAF->rtnviewfactor_vector(2));
							viewfactor_vector[1] = (C->rtnimatrix0(1, 0) * DAF->rtnviewfactor_vector(0)) + ((C->rtnimatrix0(1, 1) * DAF->rtnviewfactor_vector(1))) + (C->rtnimatrix0(1, 2) * DAF->rtnviewfactor_vector(2));
							viewfactor_vector[2] = (C->rtnimatrix0(2, 0) * DAF->rtnviewfactor_vector(0)) + ((C->rtnimatrix0(2, 1) * DAF->rtnviewfactor_vector(1))) + (C->rtnimatrix0(2, 2) * DAF->rtnviewfactor_vector(2));

							roughfacetreflectedpressure[0] = roughfacetreflectedpressure[0] + (subfacetreflectedpressure * DAF->rtnviewfactor(j + 1) * viewfactor_vector[0]);
							roughfacetreflectedpressure[1] = roughfacetreflectedpressure[1] + (subfacetreflectedpressure * DAF->rtnviewfactor(j + 1) * viewfactor_vector[1]);
							roughfacetreflectedpressure[2] = roughfacetreflectedpressure[2] + (subfacetreflectedpressure * DAF->rtnviewfactor(j + 1) * viewfactor_vector[2]);

							roughfacetthermalpressure[0] = roughfacetthermalpressure[0] + (subfacetthermalpressure * DAF->rtnviewfactor(j + 1) * viewfactor_vector[0]);
							roughfacetthermalpressure[1] = roughfacetthermalpressure[1] + (subfacetthermalpressure * DAF->rtnviewfactor(j + 1) * viewfactor_vector[1]);
							roughfacetthermalpressure[2] = roughfacetthermalpressure[2] + (subfacetthermalpressure * DAF->rtnviewfactor(j + 1) * viewfactor_vector[2]);

						}
					}

				}

				roughfacetsolarpressure[0] = roughfacetsolarpressure[0] * (C->rtnarea() / S.rtntotalarea());
				roughfacetsolarpressure[1] = roughfacetsolarpressure[1] * (C->rtnarea() / S.rtntotalarea());
				roughfacetsolarpressure[2] = roughfacetsolarpressure[2] * (C->rtnarea() / S.rtntotalarea());

				helioroughfacetsolarpressure[0] = (roughfacetsolarpressure[0] * C->rtnimatrix(step, 0, 0)) + (roughfacetsolarpressure[1] * C->rtnimatrix(step, 0, 1)) + (roughfacetsolarpressure[2] * C->rtnimatrix(step, 0, 2));
				helioroughfacetsolarpressure[1] = (roughfacetsolarpressure[0] * C->rtnimatrix(step, 1, 0)) + (roughfacetsolarpressure[1] * C->rtnimatrix(step, 1, 1)) + (roughfacetsolarpressure[2] * C->rtnimatrix(step, 1, 2));
				helioroughfacetsolarpressure[2] = (roughfacetsolarpressure[0] * C->rtnimatrix(step, 2, 0)) + (roughfacetsolarpressure[1] * C->rtnimatrix(step, 2, 1)) + (roughfacetsolarpressure[2] * C->rtnimatrix(step, 2, 2));

				C->rtnhelroughsolarpressurevec(0) = C->rtnhelroughsolarpressurevec(0) + ((1.0 / NTAU) * helioroughfacetsolarpressure[0]);
				C->rtnhelroughsolarpressurevec(1) = C->rtnhelroughsolarpressurevec(1) + ((1.0 / NTAU) * helioroughfacetsolarpressure[1]);
				C->rtnhelroughsolarpressurevec(2) = C->rtnhelroughsolarpressurevec(2) + ((1.0 / NTAU) * helioroughfacetsolarpressure[2]);

				equroughfacetsolarpressure[0] = (helioroughfacetsolarpressure[0] * A.rtnimatrix(0, 0)) + (helioroughfacetsolarpressure[1] * A.rtnimatrix(0, 1)) + (helioroughfacetsolarpressure[2] * A.rtnimatrix(0, 2));
				equroughfacetsolarpressure[1] = (helioroughfacetsolarpressure[0] * A.rtnimatrix(1, 0)) + (helioroughfacetsolarpressure[1] * A.rtnimatrix(1, 1)) + (helioroughfacetsolarpressure[2] * A.rtnimatrix(1, 2));
				equroughfacetsolarpressure[2] = (helioroughfacetsolarpressure[0] * A.rtnimatrix(2, 0)) + (helioroughfacetsolarpressure[1] * A.rtnimatrix(2, 1)) + (helioroughfacetsolarpressure[2] * A.rtnimatrix(2, 2));

				roughfacetsolartorq[0] = (C->rtnmidpos(step, 1) * equroughfacetsolarpressure[2]) - (C->rtnmidpos(step, 2) * equroughfacetsolarpressure[1]);
				roughfacetsolartorq[1] = (C->rtnmidpos(step, 2) * equroughfacetsolarpressure[0]) - (C->rtnmidpos(step, 0) * equroughfacetsolarpressure[2]);
				roughfacetsolartorq[2] = (C->rtnmidpos(step, 0) * equroughfacetsolarpressure[1]) - (C->rtnmidpos(step, 1) * equroughfacetsolarpressure[0]);

				C->rtnroughsolartorquevec(0) = C->rtnroughsolartorquevec(0) + ((1.0 / NTAU) * roughfacetsolartorq[0]);
				C->rtnroughsolartorquevec(1) = C->rtnroughsolartorquevec(1) + ((1.0 / NTAU) * roughfacetsolartorq[1]);
				C->rtnroughsolartorquevec(2) = C->rtnroughsolartorquevec(2) + ((1.0 / NTAU) * roughfacetsolartorq[2]);

				roughfacetreflectedpressure[0] = roughfacetreflectedpressure[0] * (C->rtnarea() / S.rtntotalarea());
				roughfacetreflectedpressure[1] = roughfacetreflectedpressure[1] * (C->rtnarea() / S.rtntotalarea());
				roughfacetreflectedpressure[2] = roughfacetreflectedpressure[2] * (C->rtnarea() / S.rtntotalarea());

				helioroughfacetreflectedpressure[0] = (roughfacetreflectedpressure[0] * C->rtnimatrix(step, 0, 0)) + (roughfacetreflectedpressure[1] * C->rtnimatrix(step, 0, 1)) + (roughfacetreflectedpressure[2] * C->rtnimatrix(step, 0, 2));
				helioroughfacetreflectedpressure[1] = (roughfacetreflectedpressure[0] * C->rtnimatrix(step, 1, 0)) + (roughfacetreflectedpressure[1] * C->rtnimatrix(step, 1, 1)) + (roughfacetreflectedpressure[2] * C->rtnimatrix(step, 1, 2));
				helioroughfacetreflectedpressure[2] = (roughfacetreflectedpressure[0] * C->rtnimatrix(step, 2, 0)) + (roughfacetreflectedpressure[1] * C->rtnimatrix(step, 2, 1)) + (roughfacetreflectedpressure[2] * C->rtnimatrix(step, 2, 2));

				C->rtnhelroughreflectedpressurevec(0) = C->rtnhelroughreflectedpressurevec(0) + ((1.0 / NTAU) * helioroughfacetreflectedpressure[0]);
				C->rtnhelroughreflectedpressurevec(1) = C->rtnhelroughreflectedpressurevec(1) + ((1.0 / NTAU) * helioroughfacetreflectedpressure[1]);
				C->rtnhelroughreflectedpressurevec(2) = C->rtnhelroughreflectedpressurevec(2) + ((1.0 / NTAU) * helioroughfacetreflectedpressure[2]);

				equroughfacetreflectedpressure[0] = (helioroughfacetreflectedpressure[0] * A.rtnimatrix(0, 0)) + (helioroughfacetreflectedpressure[1] * A.rtnimatrix(0, 1)) + (helioroughfacetreflectedpressure[2] * A.rtnimatrix(0, 2));
				equroughfacetreflectedpressure[1] = (helioroughfacetreflectedpressure[0] * A.rtnimatrix(1, 0)) + (helioroughfacetreflectedpressure[1] * A.rtnimatrix(1, 1)) + (helioroughfacetreflectedpressure[2] * A.rtnimatrix(1, 2));
				equroughfacetreflectedpressure[2] = (helioroughfacetreflectedpressure[0] * A.rtnimatrix(2, 0)) + (helioroughfacetreflectedpressure[1] * A.rtnimatrix(2, 1)) + (helioroughfacetreflectedpressure[2] * A.rtnimatrix(2, 2));

				roughfacetreflectedtorq[0] = (C->rtnmidpos(step, 1) * equroughfacetreflectedpressure[2]) - (C->rtnmidpos(step, 2) * equroughfacetreflectedpressure[1]);
				roughfacetreflectedtorq[1] = (C->rtnmidpos(step, 2) * equroughfacetreflectedpressure[0]) - (C->rtnmidpos(step, 0) * equroughfacetreflectedpressure[2]);
				roughfacetreflectedtorq[2] = (C->rtnmidpos(step, 0) * equroughfacetreflectedpressure[1]) - (C->rtnmidpos(step, 1) * equroughfacetreflectedpressure[0]);

				C->rtnroughreflectedtorquevec(0) = C->rtnroughreflectedtorquevec(0) + ((1.0 / NTAU) * roughfacetreflectedtorq[0]);
				C->rtnroughreflectedtorquevec(1) = C->rtnroughreflectedtorquevec(1) + ((1.0 / NTAU) * roughfacetreflectedtorq[1]);
				C->rtnroughreflectedtorquevec(2) = C->rtnroughreflectedtorquevec(2) + ((1.0 / NTAU) * roughfacetreflectedtorq[2]);

				roughfacetthermalpressure[0] = roughfacetthermalpressure[0] * (C->rtnarea() / S.rtntotalarea());
				roughfacetthermalpressure[1] = roughfacetthermalpressure[1] * (C->rtnarea() / S.rtntotalarea());
				roughfacetthermalpressure[2] = roughfacetthermalpressure[2] * (C->rtnarea() / S.rtntotalarea());

				helioroughfacetthermalpressure[0] = (roughfacetthermalpressure[0] * C->rtnimatrix(step, 0, 0)) + (roughfacetthermalpressure[1] * C->rtnimatrix(step, 0, 1)) + (roughfacetthermalpressure[2] * C->rtnimatrix(step, 0, 2));
				helioroughfacetthermalpressure[1] = (roughfacetthermalpressure[0] * C->rtnimatrix(step, 1, 0)) + (roughfacetthermalpressure[1] * C->rtnimatrix(step, 1, 1)) + (roughfacetthermalpressure[2] * C->rtnimatrix(step, 1, 2));
				helioroughfacetthermalpressure[2] = (roughfacetthermalpressure[0] * C->rtnimatrix(step, 2, 0)) + (roughfacetthermalpressure[1] * C->rtnimatrix(step, 2, 1)) + (roughfacetthermalpressure[2] * C->rtnimatrix(step, 2, 2));

				C->rtnhelroughthermalpressurevec(0) = C->rtnhelroughthermalpressurevec(0) + ((1.0 / NTAU) * helioroughfacetthermalpressure[0]);
				C->rtnhelroughthermalpressurevec(1) = C->rtnhelroughthermalpressurevec(1) + ((1.0 / NTAU) * helioroughfacetthermalpressure[1]);
				C->rtnhelroughthermalpressurevec(2) = C->rtnhelroughthermalpressurevec(2) + ((1.0 / NTAU) * helioroughfacetthermalpressure[2]);

				equroughfacetthermalpressure[0] = (helioroughfacetthermalpressure[0] * A.rtnimatrix(0, 0)) + (helioroughfacetthermalpressure[1] * A.rtnimatrix(0, 1)) + (helioroughfacetthermalpressure[2] * A.rtnimatrix(0, 2));
				equroughfacetthermalpressure[1] = (helioroughfacetthermalpressure[0] * A.rtnimatrix(1, 0)) + (helioroughfacetthermalpressure[1] * A.rtnimatrix(1, 1)) + (helioroughfacetthermalpressure[2] * A.rtnimatrix(1, 2));
				equroughfacetthermalpressure[2] = (helioroughfacetthermalpressure[0] * A.rtnimatrix(2, 0)) + (helioroughfacetthermalpressure[1] * A.rtnimatrix(2, 1)) + (helioroughfacetthermalpressure[2] * A.rtnimatrix(2, 2));

				roughfacetthermaltorq[0] = (C->rtnmidpos(step, 1) * equroughfacetthermalpressure[2]) - (C->rtnmidpos(step, 2) * equroughfacetthermalpressure[1]);
				roughfacetthermaltorq[1] = (C->rtnmidpos(step, 2) * equroughfacetthermalpressure[0]) - (C->rtnmidpos(step, 0) * equroughfacetthermalpressure[2]);
				roughfacetthermaltorq[2] = (C->rtnmidpos(step, 0) * equroughfacetthermalpressure[1]) - (C->rtnmidpos(step, 1) * equroughfacetthermalpressure[0]);

				C->rtnroughthermaltorquevec(0) = C->rtnroughthermaltorquevec(0) + ((1.0 / NTAU) * roughfacetthermaltorq[0]);
				C->rtnroughthermaltorquevec(1) = C->rtnroughthermaltorquevec(1) + ((1.0 / NTAU) * roughfacetthermaltorq[1]);
				C->rtnroughthermaltorquevec(2) = C->rtnroughthermaltorquevec(2) + ((1.0 / NTAU) * roughfacetthermaltorq[2]);

			}

			// Shape facets

			astfacetsolarpressure = A.rtnsolar() * cos(C->rtnillumangle(step)) * C->rtnshadow(step) * C->rtnarea() / cvel;
			smoothfacetsolarpressure[0] = astfacetsolarpressure * -1.0 * C->rtnillumvector(step, 0) / _t;
			smoothfacetsolarpressure[1] = astfacetsolarpressure * -1.0 * C->rtnillumvector(step, 1) / _t;
			smoothfacetsolarpressure[2] = astfacetsolarpressure * -1.0 * C->rtnillumvector(step, 2) / _t;

			heliosmoothfacetsolarpressure[0] = (smoothfacetsolarpressure[0] * C->rtnimatrix(step, 0, 0)) + (smoothfacetsolarpressure[1] * C->rtnimatrix(step, 0, 1)) + (smoothfacetsolarpressure[2] * C->rtnimatrix(step, 0, 2));
			heliosmoothfacetsolarpressure[1] = (smoothfacetsolarpressure[0] * C->rtnimatrix(step, 1, 0)) + (smoothfacetsolarpressure[1] * C->rtnimatrix(step, 1, 1)) + (smoothfacetsolarpressure[2] * C->rtnimatrix(step, 1, 2));
			heliosmoothfacetsolarpressure[2] = (smoothfacetsolarpressure[0] * C->rtnimatrix(step, 2, 0)) + (smoothfacetsolarpressure[1] * C->rtnimatrix(step, 2, 1)) + (smoothfacetsolarpressure[2] * C->rtnimatrix(step, 2, 2));

			C->rtnhelsmoothsolarpressurevec(0) = C->rtnhelsmoothsolarpressurevec(0) + ((1.0 / NTAU) * heliosmoothfacetsolarpressure[0]);
			C->rtnhelsmoothsolarpressurevec(1) = C->rtnhelsmoothsolarpressurevec(1) + ((1.0 / NTAU) * heliosmoothfacetsolarpressure[1]);
			C->rtnhelsmoothsolarpressurevec(2) = C->rtnhelsmoothsolarpressurevec(2) + ((1.0 / NTAU) * heliosmoothfacetsolarpressure[2]);

			equsmoothfacetsolarpressure[0] = (heliosmoothfacetsolarpressure[0] * A.rtnimatrix(0, 0)) + (heliosmoothfacetsolarpressure[1] * A.rtnimatrix(0, 1)) + (heliosmoothfacetsolarpressure[2] * A.rtnimatrix(0, 2));
			equsmoothfacetsolarpressure[1] = (heliosmoothfacetsolarpressure[0] * A.rtnimatrix(1, 0)) + (heliosmoothfacetsolarpressure[1] * A.rtnimatrix(1, 1)) + (heliosmoothfacetsolarpressure[2] * A.rtnimatrix(1, 2));
			equsmoothfacetsolarpressure[2] = (heliosmoothfacetsolarpressure[0] * A.rtnimatrix(2, 0)) + (heliosmoothfacetsolarpressure[1] * A.rtnimatrix(2, 1)) + (heliosmoothfacetsolarpressure[2] * A.rtnimatrix(2, 2));

			smoothfacetsolartorq[0] = (C->rtnmidpos(step, 1) * equsmoothfacetsolarpressure[2]) - (C->rtnmidpos(step, 2) * equsmoothfacetsolarpressure[1]);
			smoothfacetsolartorq[1] = (C->rtnmidpos(step, 2) * equsmoothfacetsolarpressure[0]) - (C->rtnmidpos(step, 0) * equsmoothfacetsolarpressure[2]);
			smoothfacetsolartorq[2] = (C->rtnmidpos(step, 0) * equsmoothfacetsolarpressure[1]) - (C->rtnmidpos(step, 1) * equsmoothfacetsolarpressure[0]);

			C->rtnsmoothsolartorquevec(0) = C->rtnsmoothsolartorquevec(0) + ((1.0 / NTAU) * smoothfacetsolartorq[0]);
			C->rtnsmoothsolartorquevec(1) = C->rtnsmoothsolartorquevec(1) + ((1.0 / NTAU) * smoothfacetsolartorq[1]);
			C->rtnsmoothsolartorquevec(2) = C->rtnsmoothsolartorquevec(2) + ((1.0 / NTAU) * smoothfacetsolartorq[2]);

			astfacetreflectedpressure = (2.0 / 3.0) * C->rtnalbedo() * C->rtntotalillumination(step) * C->rtnarea() / cvel;
			heliosmoothfacetreflectedpressure[0] = astfacetreflectedpressure * C->rtnradvector(step, 0);
			heliosmoothfacetreflectedpressure[1] = astfacetreflectedpressure * C->rtnradvector(step, 1);
			heliosmoothfacetreflectedpressure[2] = astfacetreflectedpressure * C->rtnradvector(step, 2);

			C->rtnhelsmoothreflectedpressurevec(0) = C->rtnhelsmoothreflectedpressurevec(0) + ((1.0 / NTAU) * heliosmoothfacetreflectedpressure[0]);
			C->rtnhelsmoothreflectedpressurevec(1) = C->rtnhelsmoothreflectedpressurevec(1) + ((1.0 / NTAU) * heliosmoothfacetreflectedpressure[1]);
			C->rtnhelsmoothreflectedpressurevec(2) = C->rtnhelsmoothreflectedpressurevec(2) + ((1.0 / NTAU) * heliosmoothfacetreflectedpressure[2]);

			equsmoothfacetreflectedpressure[0] = (heliosmoothfacetreflectedpressure[0] * A.rtnimatrix(0, 0)) + (heliosmoothfacetreflectedpressure[1] * A.rtnimatrix(0, 1)) + (heliosmoothfacetreflectedpressure[2] * A.rtnimatrix(0, 2));
			equsmoothfacetreflectedpressure[1] = (heliosmoothfacetreflectedpressure[0] * A.rtnimatrix(1, 0)) + (heliosmoothfacetreflectedpressure[1] * A.rtnimatrix(1, 1)) + (heliosmoothfacetreflectedpressure[2] * A.rtnimatrix(1, 2));
			equsmoothfacetreflectedpressure[2] = (heliosmoothfacetreflectedpressure[0] * A.rtnimatrix(2, 0)) + (heliosmoothfacetreflectedpressure[1] * A.rtnimatrix(2, 1)) + (heliosmoothfacetreflectedpressure[2] * A.rtnimatrix(2, 2));

			smoothfacetreflectedtorq[0] = (C->rtnmidpos(step, 1) * equsmoothfacetreflectedpressure[2]) - (C->rtnmidpos(step, 2) * equsmoothfacetreflectedpressure[1]);
			smoothfacetreflectedtorq[1] = (C->rtnmidpos(step, 2) * equsmoothfacetreflectedpressure[0]) - (C->rtnmidpos(step, 0) * equsmoothfacetreflectedpressure[2]);
			smoothfacetreflectedtorq[2] = (C->rtnmidpos(step, 0) * equsmoothfacetreflectedpressure[1]) - (C->rtnmidpos(step, 1) * equsmoothfacetreflectedpressure[0]);

			C->rtnsmoothreflectedtorquevec(0) = C->rtnsmoothreflectedtorquevec(0) + ((1.0 / NTAU) * smoothfacetreflectedtorq[0]);
			C->rtnsmoothreflectedtorquevec(1) = C->rtnsmoothreflectedtorquevec(1) + ((1.0 / NTAU) * smoothfacetreflectedtorq[1]);
			C->rtnsmoothreflectedtorquevec(2) = C->rtnsmoothreflectedtorquevec(2) + ((1.0 / NTAU) * smoothfacetreflectedtorq[2]);

			astfacetthermalpressure = (2.0 / 3.0) * EPS * SIGMA * pow((double)C->rtnstemp(1, step), 4.0) * C->rtnarea() / cvel;
			heliosmoothfacetthermalpressure[0] = astfacetthermalpressure * C->rtnradvector(step, 0);
			heliosmoothfacetthermalpressure[1] = astfacetthermalpressure * C->rtnradvector(step, 1);
			heliosmoothfacetthermalpressure[2] = astfacetthermalpressure * C->rtnradvector(step, 2);

			C->rtnhelsmooththermalpressurevec(0) = C->rtnhelsmooththermalpressurevec(0) + ((1.0 / NTAU) * heliosmoothfacetthermalpressure[0]);
			C->rtnhelsmooththermalpressurevec(1) = C->rtnhelsmooththermalpressurevec(1) + ((1.0 / NTAU) * heliosmoothfacetthermalpressure[1]);
			C->rtnhelsmooththermalpressurevec(2) = C->rtnhelsmooththermalpressurevec(2) + ((1.0 / NTAU) * heliosmoothfacetthermalpressure[2]);

			equsmoothfacetthermalpressure[0] = (heliosmoothfacetthermalpressure[0] * A.rtnimatrix(0, 0)) + (heliosmoothfacetthermalpressure[1] * A.rtnimatrix(0, 1)) + (heliosmoothfacetthermalpressure[2] * A.rtnimatrix(0, 2));
			equsmoothfacetthermalpressure[1] = (heliosmoothfacetthermalpressure[0] * A.rtnimatrix(1, 0)) + (heliosmoothfacetthermalpressure[1] * A.rtnimatrix(1, 1)) + (heliosmoothfacetthermalpressure[2] * A.rtnimatrix(1, 2));
			equsmoothfacetthermalpressure[2] = (heliosmoothfacetthermalpressure[0] * A.rtnimatrix(2, 0)) + (heliosmoothfacetthermalpressure[1] * A.rtnimatrix(2, 1)) + (heliosmoothfacetthermalpressure[2] * A.rtnimatrix(2, 2));

			smoothfacetthermaltorq[0] = (C->rtnmidpos(step, 1) * equsmoothfacetthermalpressure[2]) - (C->rtnmidpos(step, 2) * equsmoothfacetthermalpressure[1]);
			smoothfacetthermaltorq[1] = (C->rtnmidpos(step, 2) * equsmoothfacetthermalpressure[0]) - (C->rtnmidpos(step, 0) * equsmoothfacetthermalpressure[2]);
			smoothfacetthermaltorq[2] = (C->rtnmidpos(step, 0) * equsmoothfacetthermalpressure[1]) - (C->rtnmidpos(step, 1) * equsmoothfacetthermalpressure[0]);

			C->rtnsmooththermaltorquevec(0) = C->rtnsmooththermaltorquevec(0) + ((1.0 / NTAU) * smoothfacetthermaltorq[0]);
			C->rtnsmooththermaltorquevec(1) = C->rtnsmooththermaltorquevec(1) + ((1.0 / NTAU) * smoothfacetthermaltorq[1]);
			C->rtnsmooththermaltorquevec(2) = C->rtnsmooththermaltorquevec(2) + ((1.0 / NTAU) * smoothfacetthermaltorq[2]);

		}

	}

};


// Retrieve GPU facet properties


extern "C" void obtain_GPU_facet_properties(Facet_Properties * FP_H);	// Forward declaration
void retrieve_GPU_facet_properties(Asteroid& A) {

	Facet_Properties* FP_H;
	size_t size = ASTNUMFACE * sizeof(Facet_Properties);
	FP_H = (Facet_Properties*)malloc(size);

	obtain_GPU_facet_properties(FP_H);

	AstFacet* C;

	for (int i = 0; i != ASTNUMFACE; ++i) {
		C = A.beginf() + i;
		C->rtnalbedo() = FP_H[i].albedo;
		C->rtnsfcconverged() = FP_H[i].sfc_converged;
		C->rtninert() = FP_H[i].thermal_inertia;
	}

	free(FP_H);

};


// Code to test that shape facets have converged


void shape_converged(Asteroid& A) {

	float DT;
	int facetconverged;
	int convergence(1);
	A.rtnnumconverged() = 0;

	// Global facet convergence

	for (AstFacet* E = A.beginf(); E != A.endf(); ++E) {

		facetconverged = 1;

		for (int step = 0; step != NTAU; ++step) {
			DT = fabs(E->rtnstemp(1, step) - E->rtnstemp(0, step));
			if (DT >= TACC) {
				convergence = 0;
				facetconverged = 0;
			}
		}

		if (facetconverged == 1) {
			E->rtnastconverged() = true;
			A.rtnnumconverged() = A.rtnnumconverged() + 1;
		}

	}

	// Convergence output

	if (convergence == 1) { A.rtnconverged() = 1; }
	else { A.rtnconverged() = 0; }

};


// Code to test that roughness facets have converged


extern "C" void GPU_converged(float TACC);	// Forward declaration
void roughness_converged(Asteroid& A) {

	int convergence(1);
	A.rtnsfcnumconverged() = 0;

	// Roughness facet convergence

	if (roughness_thermal_model == 1) {

		GPU_converged(TACC);
		retrieve_GPU_facet_properties(A);

		for (AstFacet* E = A.beginf(); E != A.endf(); ++E) {
			if (E->rtnsfcconverged() == true) {
				A.rtnsfcnumconverged() = A.rtnsfcnumconverged() + 1;
			}
			else {
				convergence = 0;
			}
		}

	}

	// Convergence output

	if (convergence == 1) { A.rtnsfcconverged() = 1; }
	else { A.rtnsfcconverged() = 0; }

};


// Load GPU constants


extern "C" void initialise_GPU_constants(float rot_period_H, float viewfactors_H[SFCNUMFACE][SFCNUMFACE], GPU_Surface S_H);	// Forward declaration
void load_GPU_constants(Surface& S, Asteroid& A) {

	int i;
	float viewfactors[SFCNUMFACE][SFCNUMFACE];
	GPU_Surface S_H(SFCNUMVTX, SFCNUMFACE);
	SfcVertex* SV;
	SfcFacet* SF;

	for (i = 0; i != SFCNUMVTX; ++i) {
		SV = S.beginv() + i;
		S_H.SV[i].pos[0] = SV->rtnpos(0);
		S_H.SV[i].pos[1] = SV->rtnpos(1);
		S_H.SV[i].pos[2] = SV->rtnpos(2);
	}

	for (i = 0; i != SFCNUMFACE; ++i) {
		SF = S.beginf() + i;
		S_H.SF[i].sfcvertices[0] = SF->rtnsfcvertices(0);
		S_H.SF[i].sfcvertices[1] = SF->rtnsfcvertices(1);
		S_H.SF[i].sfcvertices[2] = SF->rtnsfcvertices(2);
		S_H.SF[i].midpos[0] = SF->rtnmidpos(0);
		S_H.SF[i].midpos[1] = SF->rtnmidpos(1);
		S_H.SF[i].midpos[2] = SF->rtnmidpos(2);
		S_H.SF[i].trivector[0][0] = SF->rtntrivector(0, 0);
		S_H.SF[i].trivector[0][1] = SF->rtntrivector(0, 1);
		S_H.SF[i].trivector[0][2] = SF->rtntrivector(0, 2);
		S_H.SF[i].trivector[1][0] = SF->rtntrivector(1, 0);
		S_H.SF[i].trivector[1][1] = SF->rtntrivector(1, 1);
		S_H.SF[i].trivector[1][2] = SF->rtntrivector(1, 2);
		S_H.SF[i].normcalc[0][0] = SF->rtnnormcalc(0, 0);
		S_H.SF[i].normcalc[0][1] = SF->rtnnormcalc(0, 1);
		S_H.SF[i].normcalc[0][2] = SF->rtnnormcalc(0, 2);
		S_H.SF[i].normcalc[1][0] = SF->rtnnormcalc(1, 0);
		S_H.SF[i].normcalc[1][1] = SF->rtnnormcalc(1, 1);
		S_H.SF[i].normcalc[1][2] = SF->rtnnormcalc(1, 2);
		S_H.SF[i].normcalc[2][0] = SF->rtnnormcalc(2, 0);
		S_H.SF[i].normcalc[2][1] = SF->rtnnormcalc(2, 1);
		S_H.SF[i].normcalc[2][2] = SF->rtnnormcalc(2, 2);
		S_H.SF[i].normal[0] = SF->rtnnormal(0);
		S_H.SF[i].normal[1] = SF->rtnnormal(1);
		S_H.SF[i].normal[2] = SF->rtnnormal(2);
	}

	for (i = 0; i != SFCNUMFACE; ++i) {
		SF = S.beginf() + i;
		for (DepSfcFacet* DSF = SF->beginf(); DSF != SF->endf(); ++DSF) {
			viewfactors[DSF->rtndepsfcfacetnum()][i] = DSF->rtnviewfactor();
		}
		viewfactors[i][i] = 0.0f;
	}

	initialise_GPU_constants(A.rtnrotperiod(), viewfactors, S_H);

};


// Load GPU facet properties


extern "C" void set_GPU_facet_properties(Facet_Properties * FP_H);
void load_GPU_facet_properties(Asteroid& A) {

	Facet_Properties* FP_H;
	size_t size = ASTNUMFACE * sizeof(Facet_Properties);
	FP_H = (Facet_Properties*)malloc(size);

	AstFacet* C;

	for (int i = 0; i != ASTNUMFACE; ++i) {
		C = A.beginf() + i;
		FP_H[i].albedo = C->rtnalbedo();
		FP_H[i].sfc_converged = C->rtnsfcconverged();
		FP_H[i].thermal_inertia = C->rtninert();
	}

	set_GPU_facet_properties(FP_H);
	free(FP_H);

};


// Load GPU surface roughness global illumination data


extern "C" void load_GPU_globalillumination(Roughness_Globalillumination * RGI_H);
void load_GPU_roughness_globalillumination(Asteroid& A) {

	Roughness_Globalillumination* RGI_H;
	size_t size = ASTNUMFACE * sizeof(Roughness_Globalillumination);
	RGI_H = (Roughness_Globalillumination*)malloc(size);

	int i, j, step;
	AstFacet* C;
	DepAstFacet* DAF;

	omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(j,step,C,DAF)
	for (i = 0; i < ASTNUMFACE; i++) {
		C = A.beginf() + i;
		for (step = 0; step != NTAU; ++step) {
			for (j = 0; j != SFCNUMFACE; ++j) {
				RGI_H[i].sfc_globalillumination[step][j] = 0.0f;
				if (C->rtnastdependentfacets() != 0) {
					for (DAF = C->begindaf(); DAF != C->enddaf(); ++DAF) {
						RGI_H[i].sfc_globalillumination[step][j] = RGI_H[i].sfc_globalillumination[step][j] + (DAF->rtnviewfactor(j + 1) * (A.beginf() + DAF->rtndepastfacetnum())->rtntotalillumination(step) * (A.beginf() + DAF->rtndepastfacetnum())->rtnalbedo());
					}
				}
			}
		}
	}

	load_GPU_globalillumination(RGI_H);
	free(RGI_H);

};


// Load GPU surface roughness global selfheating data


extern "C" void load_GPU_globalheating(Roughness_Globalheating * RGH_H);
void load_GPU_roughness_globalheating(Asteroid& A) {

	Roughness_Globalheating* RGH_H;
	size_t size = ASTNUMFACE * sizeof(Roughness_Globalheating);
	RGH_H = (Roughness_Globalheating*)malloc(size);

	int i, j, step;
	float temp;
	AstFacet* C;
	DepAstFacet* DAF;

	omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(j,step,temp,C,DAF)
	for (i = 0; i < ASTNUMFACE; i++) {
		C = A.beginf() + i;
		for (step = 0; step != NTAU; ++step) {
			for (j = 0; j != SFCNUMFACE; ++j) {
				RGH_H[i].sfc_globalheating[step][j] = 0.0f;
				if (C->rtnastdependentfacets() != 0) {
					for (DAF = C->begindaf(); DAF != C->enddaf(); ++DAF) {
						temp = (A.beginf() + DAF->rtndepastfacetnum())->rtnstemp(1, step);
						RGH_H[i].sfc_globalheating[step][j] = RGH_H[i].sfc_globalheating[step][j] + (DAF->rtnviewfactor(j + 1) * EPS * SIGMA * temp * temp * temp * temp);
					}
				}
			}
		}
	}

	load_GPU_globalheating(RGH_H);
	free(RGH_H);

};


// Perform GPU surface roughness shadow calculations


extern "C" void GPU_sfcpartialshadow(Roughness_Shadow * RS_H);	// Forward declaration
void perform_GPU_roughness_shadows(Asteroid& A) {

	Roughness_Shadow* RS_H;
	size_t size = ASTNUMFACE * sizeof(Roughness_Shadow);
	RS_H = (Roughness_Shadow*)malloc(size);

	AstFacet* C;
	SubAstFacet* SC;

	for (int i = 0; i != ASTNUMFACE; ++i) {
		C = A.beginf() + i;
		for (int step = 0; step != NTAU; ++step) {
			RS_H[i].shadow[step] = C->rtnshadow(step);
			RS_H[i].sfc_illumvector[step][0] = C->rtnsfcillumvector(step, 0);
			RS_H[i].sfc_illumvector[step][1] = C->rtnsfcillumvector(step, 1);
			RS_H[i].sfc_illumvector[step][2] = C->rtnsfcillumvector(step, 2);
			for (int j = 0; j != SFCNUMFACE; ++j) {
				SC = C->beginsf() + j;
				RS_H[i].sfc_shadow[step][j] = SC->rtnshadow(step);
				RS_H[i].sfc_illumangle[step][j] = SC->rtnillumangle(step);
			}
		}
	}

	GPU_sfcpartialshadow(RS_H);
	free(RS_H);

};


// Retrieve GPU surface roughness shadow data


extern "C" void obtain_GPU_roughness_shadows(Roughness_Shadow * RS_H);	// Forward declaration
void retrieve_GPU_roughness_shadows(Asteroid& A) {

	Roughness_Shadow* RS_H;
	size_t size = ASTNUMFACE * sizeof(Roughness_Shadow);
	RS_H = (Roughness_Shadow*)malloc(size);

	obtain_GPU_roughness_shadows(RS_H);

	AstFacet* C;
	SubAstFacet* SC;

	for (int i = 0; i != ASTNUMFACE; ++i) {
		C = A.beginf() + i;
		for (int step = 0; step != NTAU; ++step) {
			C->rtnshadow(step) = RS_H[i].shadow[step];
			C->rtnsfcillumvector(step, 0) = RS_H[i].sfc_illumvector[step][0];
			C->rtnsfcillumvector(step, 1) = RS_H[i].sfc_illumvector[step][1];
			C->rtnsfcillumvector(step, 2) = RS_H[i].sfc_illumvector[step][2];
			for (int j = 0; j != SFCNUMFACE; ++j) {
				SC = C->beginsf() + j;
				SC->rtnshadow(step) = RS_H[i].sfc_shadow[step][j];
				SC->rtnillumangle(step) = RS_H[i].sfc_illumangle[step][j];
			}
		}
	}

	free(RS_H);

};


// Retrieve GPU surface roughness temperature data


extern "C" void obtain_GPU_roughness_temperatures(Roughness_Temperatures * RT_H);	// Forward declaration
void retrieve_GPU_roughness_temperatures(Asteroid& A) {

	Roughness_Temperatures* RT_H;
	size_t size = ASTNUMFACE * sizeof(Roughness_Temperatures);
	RT_H = (Roughness_Temperatures*)malloc(size);

	obtain_GPU_roughness_temperatures(RT_H);

	AstFacet* C;
	SubAstFacet* SC;

	for (int i = 0; i != ASTNUMFACE; ++i) {
		C = A.beginf() + i;
		for (int j = 0; j != SFCNUMFACE; ++j) {
			SC = C->beginsf() + j;
			for (int step = 0; step != NTAU; ++step) {
				SC->rtnstemp(0, step) = RT_H[i].sfc_stemp0[step][j];
				SC->rtnstemp(1, step) = RT_H[i].sfc_stemp1[step][j];
			}
			for (int k = 0; k != (NX + 1); ++k) {
				SC->rtntemp(k) = RT_H[i].sfc_temp[k][j];
			}
		}
	}

	free(RT_H);

};


// Retrieve GPU surface roughness illumination data


extern "C" void obtain_GPU_roughness_illumination(Roughness_Illumination * RI_H);	// Forward declaration
void retrieve_GPU_roughness_illumination(Asteroid& A) {

	Roughness_Illumination* RI_H;
	size_t size = ASTNUMFACE * sizeof(Roughness_Illumination);
	RI_H = (Roughness_Illumination*)malloc(size);

	obtain_GPU_roughness_illumination(RI_H);

	AstFacet* C;
	SubAstFacet* SC;

	for (int i = 0; i != ASTNUMFACE; ++i) {
		C = A.beginf() + i;
		for (int j = 0; j != SFCNUMFACE; ++j) {
			SC = C->beginsf() + j;
			for (int step = 0; step != NTAU; ++step) {
				SC->rtntotalillumination(step) = RI_H[i].sfc_illumination[step][j];
			}
		}
	}

	free(RI_H);

};


// Write out global shadow map


void writeglobalshadow(Asteroid& A, string filename) {

	ofstream write(filename.c_str());

	for (AstVertex* B = A.beginv(); B != A.endv(); ++B) {
		write << B->rtnequpos(0) << "\t" << B->rtnequpos(1) << "\t" << B->rtnequpos(2) << endl;
	}

	for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
		write << C->rtnastvertices(0) << "\t" << C->rtnastvertices(1) << "\t" << C->rtnastvertices(2) << endl;
	}

	for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
		for (int step = 0; step != NTAU; ++step) {
			write << C->rtnshadow(step) << endl;
		}
	}

};


// Write out global illumination map


void writeglobalillumination(Asteroid& A, string filename) {

	ofstream write(filename.c_str());

	for (AstVertex* B = A.beginv(); B != A.endv(); ++B) {
		write << B->rtnequpos(0) << "\t" << B->rtnequpos(1) << "\t" << B->rtnequpos(2) << endl;
	}

	for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
		write << C->rtnastvertices(0) << "\t" << C->rtnastvertices(1) << "\t" << C->rtnastvertices(2) << endl;
	}

	for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
		for (int step = 0; step != NTAU; ++step) {
			write << C->rtntotalillumination(step) << endl;
		}
	}

};


// Write out flat global shape temperature map


void writeglobaltemperature(Asteroid& A, string filename) {

	ofstream write(filename.c_str());

	for (AstVertex* B = A.beginv(); B != A.endv(); ++B) {
		write << B->rtnequpos(0) << "\t" << B->rtnequpos(1) << "\t" << B->rtnequpos(2) << endl;
	}

	for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
		write << C->rtnastvertices(0) << "\t" << C->rtnastvertices(1) << "\t" << C->rtnastvertices(2) << endl;
	}

	for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
		for (int step = 0; step != NTAU; ++step) {
			write << C->rtnstemp(1, step) << endl;
		}
	}

};


// Binary write out surface roughness temperature maps


void binary_writesurfacetemperature(Surface& S, Asteroid& A, string filename) {

	ofstream write(filename.c_str(), ios::out | ios::binary);

	float pos[3];
	int vert[3];
	float val[NTAU];

	for (SfcVertex* B = S.beginv(); B != S.endv(); ++B) {
		pos[0] = B->rtnpos(0), pos[1] = B->rtnpos(1), pos[2] = B->rtnpos(2);
		write.write((char*)&pos, sizeof(float) * 3);
	}

	for (SfcFacet* C = S.beginf(); C != S.endf(); ++C) {
		vert[0] = C->rtnsfcvertices(0), vert[1] = C->rtnsfcvertices(1), vert[2] = C->rtnsfcvertices(2);
		write.write((char*)&vert, sizeof(int) * 3);
	}

	for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
		for (SubAstFacet* D = C->beginsf(); D != C->endsf(); ++D) {
			for (int step = 0; step != NTAU; ++step) {
				val[step] = D->rtnstemp(1, step);
			}
			write.write((char*)&val, sizeof(float) * NTAU);
		}
	}

};


// Write out Yarkovsky results


void write_yarkovsky_results(Asteroid& A, string smooth_filename, string rough_filename) {

	ofstream smooth(smooth_filename.c_str());
	for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
		smooth << C->rtnhelsmoothsolarpressurevec(0) << "\t" << C->rtnhelsmoothsolarpressurevec(1) << "\t" << C->rtnhelsmoothsolarpressurevec(2) << "\t";
		smooth << C->rtnhelsmoothreflectedpressurevec(0) << "\t" << C->rtnhelsmoothreflectedpressurevec(1) << "\t" << C->rtnhelsmoothreflectedpressurevec(2) << "\t";
		smooth << C->rtnhelsmooththermalpressurevec(0) << "\t" << C->rtnhelsmooththermalpressurevec(1) << "\t" << C->rtnhelsmooththermalpressurevec(2) << endl;
	}

	if (roughness_thermal_model == 1) {
		ofstream rough(rough_filename.c_str());
		for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
			rough << C->rtnhelroughsolarpressurevec(0) << "\t" << C->rtnhelroughsolarpressurevec(1) << "\t" << C->rtnhelroughsolarpressurevec(2) << "\t";
			rough << C->rtnhelroughreflectedpressurevec(0) << "\t" << C->rtnhelroughreflectedpressurevec(1) << "\t" << C->rtnhelroughreflectedpressurevec(2) << "\t";
			rough << C->rtnhelroughthermalpressurevec(0) << "\t" << C->rtnhelroughthermalpressurevec(1) << "\t" << C->rtnhelroughthermalpressurevec(2) << endl;
		}
	}

};


// Write out YORP results


void write_yorp_results(Asteroid& A, string smooth_filename, string rough_filename) {

	ofstream smooth(smooth_filename.c_str());
	for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
		smooth << C->rtnsmoothsolartorquevec(0) << "\t" << C->rtnsmoothsolartorquevec(1) << "\t" << C->rtnsmoothsolartorquevec(2) << "\t";
		smooth << C->rtnsmoothreflectedtorquevec(0) << "\t" << C->rtnsmoothreflectedtorquevec(1) << "\t" << C->rtnsmoothreflectedtorquevec(2) << "\t";
		smooth << C->rtnsmooththermaltorquevec(0) << "\t" << C->rtnsmooththermaltorquevec(1) << "\t" << C->rtnsmooththermaltorquevec(2) << endl;
	}

	if (roughness_thermal_model == 1) {
		ofstream rough(rough_filename.c_str());
		for (AstFacet* C = A.beginf(); C != A.endf(); ++C) {
			rough << C->rtnroughsolartorquevec(0) << "\t" << C->rtnroughsolartorquevec(1) << "\t" << C->rtnroughsolartorquevec(2) << "\t";
			rough << C->rtnroughreflectedtorquevec(0) << "\t" << C->rtnroughreflectedtorquevec(1) << "\t" << C->rtnroughreflectedtorquevec(2) << "\t";
			rough << C->rtnroughthermaltorquevec(0) << "\t" << C->rtnroughthermaltorquevec(1) << "\t" << C->rtnroughthermaltorquevec(2) << endl;
		}
	}

};


// Other forward declarations


extern "C" void set_cuda_device();
extern "C" void initialise_GPU_facet_properties_arrays();
extern "C" void initialise_GPU_shadow_arrays();
extern "C" void initialise_GPU_illumination_arrays();
extern "C" void initialise_GPU_temperature_arrays();
extern "C" void initialise_GPU_globalillumination_arrays();
extern "C" void initialise_GPU_globalheating_arrays();
extern "C" void freeup_GPU_shadow_arrays();
extern "C" void freeup_GPU_temperature_arrays();
extern "C" void freeup_GPU_globalillumination_arrays();
extern "C" void freeup_GPU_globalheating_arrays();
extern "C" void GPU_single_scattering(float solar);
extern "C" void GPU_multiple_scattering(float solar, float scattering_accuracy);
extern "C" void GPU_zero_temperatures();
extern "C" void GPU_initialise_temperatures(float EPS, float SIGMA);
extern "C" void GPU_inertia_difference(float EPS, float SIGMA, float TACC);
extern "C" void GPU_inst_equilibrium(float EPS, float SIGMA);


// Standard model for running simulation once at a specified position


void single_model(Asteroid& A, Surface& S) {

	AstFacet* C;
	int step, rev, TI_step, num_TI_steps, initialised_temp;
	float t_inertia;
	string status_filename, filename1, filename2;
	status_filename = shape_filename + "_model_status.txt";
	ofstream status(status_filename.c_str());

	status << "Start time" << "\t" << currentDateTime() << endl;

	A.reset_asteroid_conditions();
	asteroidgeometry(A);
	if (roughness_thermal_model == 1) { load_GPU_facet_properties(A); }

	if (shape_shadowing == 1) {
		omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic)
		for (step = 0; step < NTAU; step++) {
			astshadow(A, step);
			status << "Asteroid shadow step " << step << " calculated" << "\t" << currentDateTime() << endl;
		}
	}

	if (roughness_thermal_model == 1) {
		surfacegeometry(S, A);
		calculate_surface_global_viewfactors(S, A);
		load_GPU_constants(S, A);
	}

	status << "Initialised geometry" << "\t" << currentDateTime() << endl;

	if (roughness_thermal_model == 1) {
		initialise_GPU_shadow_arrays();
		initialise_GPU_globalillumination_arrays();
		perform_GPU_roughness_shadows(A);
		retrieve_GPU_roughness_shadows(A);
		status << "Roughness shadows calculated" << "\t" << currentDateTime() << endl;
	}

	singlescattering(A);
	if (multiple_scattering == 1) { multiplescattering(A); }

	status << "Shape scattering complete" << "\t" << currentDateTime() << endl;

	if (roughness_thermal_model == 1) {
		load_GPU_roughness_globalillumination(A);
		status << "Roughness global illumination data loaded" << "\t" << currentDateTime() << endl;
		GPU_single_scattering(A.rtnsolar());
		if (multiple_scattering == 1) { GPU_multiple_scattering(A.rtnsolar(), scattering_accuracy); }
		retrieve_GPU_roughness_illumination(A);
		status << "Roughness scattering complete" << "\t" << currentDateTime() << endl;
	}

	if (roughness_thermal_model == 1) {
		freeup_GPU_shadow_arrays();
		freeup_GPU_globalillumination_arrays();
		initialise_GPU_temperature_arrays();
		initialise_GPU_globalheating_arrays();
	}

	initialised_temp = 0;
	num_TI_steps = (int)((maximum_thermal_inertia - minimum_thermal_inertia) / thermal_inertia_step);

	for (TI_step = 0; TI_step != num_TI_steps + 1; ++TI_step) {
		//	for (TI_step=0;TI_step!=98;++TI_step) {

		A.reset_asteroid_status();
		A.reset_asteroid_yark_yorp();

		stringstream inert;
		t_inertia = minimum_thermal_inertia + (TI_step * thermal_inertia_step);
		//		t_inertia=thermal_inertia[TI_step];
		inert << t_inertia;
		//		inert << TI_step;
		for (C = A.beginf(); C != A.endf(); ++C) {
			C->rtninert() = t_inertia;
		}
		if (roughness_thermal_model == 1) { load_GPU_facet_properties(A); }

		if (t_inertia == 0.0f) {
			if (roughness_thermal_model == 1) {
				load_GPU_roughness_globalheating(A);
				status << "Roughness global heating data loaded" << "\t" << currentDateTime() << endl;
				GPU_zero_temperatures();
				status << "Zeroed roughness temperatures" << "\t" << currentDateTime() << endl;
			}
		}
		else {
			if (initialised_temp == 0) {
				initialisetemperatures(A);
				status << "Initialised shape temperatures" << "\t" << currentDateTime() << endl;
				if (roughness_thermal_model == 1) {
					load_GPU_roughness_globalheating(A);
					status << "Roughness global heating data loaded" << "\t" << currentDateTime() << endl;
					GPU_initialise_temperatures(EPS, SIGMA);
					status << "Initialised roughness temperatures" << "\t" << currentDateTime() << endl;
				}
				initialised_temp = 1;
			}
		}

		status << "Running model for a thermal inertia of" << "\t" << t_inertia << "\t" << "Time" << "\t" << currentDateTime() << endl;

		rev = 0;
		while (A.rtnconverged() == 0 && rev < 301) {

			if (t_inertia == 0.0f) {
				instequilibrium(A);
			}
			else {
				inertiadifference(A);
			}

			if (rev > 5) { shape_converged(A); }

			status << "Number of shape facets converged during revolution " << rev << " is" << "\t" << A.rtnnumconverged() << "\t";
			status << "Shape facet convergence completion is" << "\t" << (100.0 * A.rtnnumconverged()) / ASTNUMFACE << "\t";
			status << "Time" << "\t" << currentDateTime() << endl;

			rev = rev + 1;

		}

		if (roughness_thermal_model == 1) {

			load_GPU_roughness_globalheating(A);
			status << "Roughness global heating data loaded" << "\t" << currentDateTime() << endl;

			rev = 0;
			while (A.rtnsfcconverged() == 0 && rev < 301) {

				if (t_inertia == 0.0f) {
					GPU_inst_equilibrium(EPS, SIGMA);
				}
				else {
					GPU_inertia_difference(EPS, SIGMA, TACC);
				}

				if (rev > 5) { roughness_converged(A); }

				status << "Number of roughness facets converged during revolution " << rev << " is" << "\t" << A.rtnsfcnumconverged() << "\t";
				status << "Roughness facet convergence completion is" << "\t" << (100.0 * A.rtnsfcnumconverged()) / ASTNUMFACE << "\t";
				status << "Time" << "\t" << currentDateTime() << endl;

				rev = rev + 1;

			}

		}

		if (roughness_thermal_model == 1) { retrieve_GPU_roughness_temperatures(A); }
		if (perform_yark_yorp == 1) { resolved_yark_yorp(A, S); }

		if (write_shape_temperature == 1) {
			filename1 = shape_filename + "_temperature(" + inert.str() + ").txt";
			writeglobaltemperature(A, filename1);
		}
		if (write_roughness_temperature == 1 && roughness_thermal_model == 1) {
			filename1 = shape_filename + "_roughness_temperature(" + inert.str() + ").dat";
			binary_writesurfacetemperature(S, A, filename1);
		}
		if (write_yarkovsky == 1) {
			filename1 = "yarkovsky_smooth_" + shape_filename + "(" + inert.str() + ").txt";
			filename2 = "yarkovsky_rough_" + shape_filename + "(" + inert.str() + ").txt";
			write_yarkovsky_results(A, filename1, filename2);
		}
		if (write_yorp == 1) {
			filename1 = "yorp_smooth_" + shape_filename + "(" + inert.str() + ").txt";
			filename2 = "yorp_rough_" + shape_filename + "(" + inert.str() + ").txt";
			write_yorp_results(A, filename1, filename2);
		}

	}

	if (roughness_thermal_model == 1) {
		freeup_GPU_temperature_arrays();
		freeup_GPU_globalheating_arrays();
	}

	if (write_shape_shadow == 1) {
		filename1 = shape_filename + "_shadow.txt";
		writeglobalshadow(A, filename1);
	}
	if (write_shape_illumination == 1) {
		filename1 = shape_filename + "_illumination.txt";
		writeglobalillumination(A, filename1);
	}

};


// Orbital model for running simulation at different orbital points and thermal inertias


void orbital_model(Asteroid& A, Surface& S) {

	AstFacet* C;
	int TA, step, rev, TI_step, num_TI_steps, initialised_temp;
	float t_inertia;
	string status_filename, filename1, filename2;

	for (int i = 0; i != NORB; ++i) {
		//	for (int i=11;i!=NORB;++i) {

		stringstream pos;
		TA = 3600 * i / NORB;
		pos << TA;
		status_filename = shape_filename + "_model_status_TA_" + pos.str() + ".txt";
		ofstream status(status_filename.c_str());
		status << "Start time" << "\t" << currentDateTime() << endl;

		A.reset_asteroid_conditions();
		orbital_position_geometry(A, i);
		asteroidgeometry(A);
		if (roughness_thermal_model == 1) { load_GPU_facet_properties(A); }

		if (shape_shadowing == 1) {
			omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic)
			for (step = 0; step < NTAU; step++) {
				astshadow(A, step);
				status << "Asteroid shadow step " << step << " calculated" << "\t" << currentDateTime() << endl;
			}
		}

		if (roughness_thermal_model == 1) {
			surfacegeometry(S, A);
			calculate_surface_global_viewfactors(S, A);
			load_GPU_constants(S, A);
		}

		status << "Initialised geometry" << "\t" << currentDateTime() << endl;

		if (roughness_thermal_model == 1) {
			initialise_GPU_shadow_arrays();
			initialise_GPU_globalillumination_arrays();
			perform_GPU_roughness_shadows(A);
			retrieve_GPU_roughness_shadows(A);
			status << "Roughness shadows calculated" << "\t" << currentDateTime() << endl;
		}

		singlescattering(A);
		if (multiple_scattering == 1) { multiplescattering(A); }

		status << "Shape scattering complete" << "\t" << currentDateTime() << endl;

		if (roughness_thermal_model == 1) {
			load_GPU_roughness_globalillumination(A);
			status << "Roughness global illumination data loaded" << "\t" << currentDateTime() << endl;
			GPU_single_scattering(A.rtnsolar());
			if (multiple_scattering == 1) { GPU_multiple_scattering(A.rtnsolar(), scattering_accuracy); }
			retrieve_GPU_roughness_illumination(A);
			status << "Roughness scattering complete" << "\t" << currentDateTime() << endl;
		}

		if (roughness_thermal_model == 1) {
			freeup_GPU_shadow_arrays();
			freeup_GPU_globalillumination_arrays();
			initialise_GPU_temperature_arrays();
			initialise_GPU_globalheating_arrays();
		}

		initialised_temp = 0;
		num_TI_steps = (int)((maximum_thermal_inertia - minimum_thermal_inertia) / thermal_inertia_step);

		for (TI_step = 0; TI_step != num_TI_steps + 1; ++TI_step) {

			A.reset_asteroid_status();
			A.reset_asteroid_yark_yorp();

			stringstream inert;
			t_inertia = minimum_thermal_inertia + (TI_step * thermal_inertia_step);
			inert << t_inertia;
			//			t_inertia=thermal_inertia0*pow(((float)A.rtnpolcentre(0)/r0),alpha);
			//			t_inertia=t_inertia*pow(((float)A.rtnpolcentre(0)/r0),alpha);
			//			inert << 300;
			for (C = A.beginf(); C != A.endf(); ++C) {
				C->rtninert() = t_inertia;
			}
			if (roughness_thermal_model == 1) { load_GPU_facet_properties(A); }

			if (t_inertia == 0.0f) {
				if (roughness_thermal_model == 1) {
					load_GPU_roughness_globalheating(A);
					status << "Roughness global heating data loaded" << "\t" << currentDateTime() << endl;
					GPU_zero_temperatures();
					status << "Zeroed roughness temperatures" << "\t" << currentDateTime() << endl;
				}
			}
			else {
				if (initialised_temp == 0) {
					initialisetemperatures(A);
					status << "Initialised shape temperatures" << "\t" << currentDateTime() << endl;
					if (roughness_thermal_model == 1) {
						load_GPU_roughness_globalheating(A);
						status << "Roughness global heating data loaded" << "\t" << currentDateTime() << endl;
						GPU_initialise_temperatures(EPS, SIGMA);
						status << "Initialised roughness temperatures" << "\t" << currentDateTime() << endl;
					}
					initialised_temp = 1;
				}
			}

			status << "Running model for a thermal inertia of" << "\t" << t_inertia << "\t" << "Time" << "\t" << currentDateTime() << endl;

			rev = 0;
			while (A.rtnconverged() == 0 && rev < 301) {

				if (t_inertia == 0.0f) {
					instequilibrium(A);
				}
				else {
					inertiadifference(A);
				}

				if (rev > 5) { shape_converged(A); }

				status << "Number of shape facets converged during revolution " << rev << " is" << "\t" << A.rtnnumconverged() << "\t";
				status << "Shape facet convergence completion is" << "\t" << (100.0 * A.rtnnumconverged()) / ASTNUMFACE << "\t";
				status << "Time" << "\t" << currentDateTime() << endl;

				rev = rev + 1;

			}

			if (roughness_thermal_model == 1) {

				load_GPU_roughness_globalheating(A);
				status << "Roughness global heating data loaded" << "\t" << currentDateTime() << endl;

				rev = 0;
				while (A.rtnsfcconverged() == 0 && rev < 301) {

					if (t_inertia == 0.0f) {
						GPU_inst_equilibrium(EPS, SIGMA);
					}
					else {
						GPU_inertia_difference(EPS, SIGMA, TACC);
					}

					if (rev > 5) { roughness_converged(A); }

					status << "Number of roughness facets converged during revolution " << rev << " is" << "\t" << A.rtnsfcnumconverged() << "\t";
					status << "Roughness facet convergence completion is" << "\t" << (100.0 * A.rtnsfcnumconverged()) / ASTNUMFACE << "\t";
					status << "Time" << "\t" << currentDateTime() << endl;

					rev = rev + 1;

				}

			}

			if (roughness_thermal_model == 1) { retrieve_GPU_roughness_temperatures(A); }
			if (perform_yark_yorp == 1) { resolved_yark_yorp(A, S); }

			if (write_shape_temperature == 1) {
				filename1 = shape_filename + "_temperature_TA_" + pos.str() + "(" + inert.str() + ").txt";
				writeglobaltemperature(A, filename1);
			}
			if (write_roughness_temperature == 1 && roughness_thermal_model == 1) {
				filename1 = shape_filename + "_roughness_temperature_TA_" + pos.str() + "(" + inert.str() + ").dat";
				binary_writesurfacetemperature(S, A, filename1);
			}
			if (write_yarkovsky == 1) {
				filename1 = "yarkovsky_smooth_" + shape_filename + "_TA_" + pos.str() + "(" + inert.str() + ").txt";
				filename2 = "yarkovsky_rough_" + shape_filename + "_TA_" + pos.str() + "(" + inert.str() + ").txt";
				write_yarkovsky_results(A, filename1, filename2);
			}
			if (write_yorp == 1) {
				filename1 = "yorp_smooth_" + shape_filename + "_TA_" + pos.str() + "(" + inert.str() + ").txt";
				filename2 = "yorp_rough_" + shape_filename + "_TA_" + pos.str() + "(" + inert.str() + ").txt";
				write_yorp_results(A, filename1, filename2);
			}

		}

		if (roughness_thermal_model == 1) {
			freeup_GPU_temperature_arrays();
			freeup_GPU_globalheating_arrays();
		}

		if (write_shape_shadow == 1) {
			filename1 = shape_filename + "_shadow_TA_" + pos.str() + ".txt";
			writeglobalshadow(A, filename1);
		}
		if (write_shape_illumination == 1) {
			filename1 = shape_filename + "_illumination_TA_" + pos.str() + ".txt";
			writeglobalillumination(A, filename1);
		}

	}

};


// Obliquity model for running simulation at different obliquities and orbital points


void obliquity_model(Asteroid& A, Surface& S) {

	int TA, step, rev;
	double obliquity;
	string status_filename, filename1, filename2;

	//	for (int j=1;j!=(NOBL+1);++j) {
	for (int j = 1; j != 18; ++j) {

		obliquity = 90.0 - (j * (90.0 / NOBL));
		A.rtnrotvector(0) = obliquity * (PI / 180.0);
		A.rtnrotvector(1) = 0.0;
		stringstream obl;
		obl << (int)obliquity;

		for (int i = 0; i != NORB; ++i) {

			stringstream pos;
			TA = 3600 * i / NORB;
			pos << TA;
			status_filename = shape_filename + "_model_status_TA_" + pos.str() + "(" + obl.str() + ").txt";
			ofstream status(status_filename.c_str());
			status << "Start time" << "\t" << currentDateTime() << endl;

			A.reset_asteroid_conditions();
			orbital_position_geometry(A, i);
			asteroidgeometry(A);
			if (roughness_thermal_model == 1) { load_GPU_facet_properties(A); }

			if (shape_shadowing == 1) {
				omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic)
				for (step = 0; step < NTAU; step++) {
					astshadow(A, step);
					status << "Asteroid shadow step " << step << " calculated" << "\t" << currentDateTime() << endl;
				}
			}

			if (roughness_thermal_model == 1) {
				surfacegeometry(S, A);
				calculate_surface_global_viewfactors(S, A);
				load_GPU_constants(S, A);
			}

			status << "Initialised geometry" << "\t" << currentDateTime() << endl;

			if (roughness_thermal_model == 1) {
				initialise_GPU_shadow_arrays();
				initialise_GPU_globalillumination_arrays();
				perform_GPU_roughness_shadows(A);
				retrieve_GPU_roughness_shadows(A);
				status << "Roughness shadows calculated" << "\t" << currentDateTime() << endl;
			}

			singlescattering(A);
			if (multiple_scattering == 1) { multiplescattering(A); }

			status << "Shape scattering complete" << "\t" << currentDateTime() << endl;

			if (roughness_thermal_model == 1) {
				load_GPU_roughness_globalillumination(A);
				status << "Roughness global illumination data loaded" << "\t" << currentDateTime() << endl;
				GPU_single_scattering(A.rtnsolar());
				if (multiple_scattering == 1) { GPU_multiple_scattering(A.rtnsolar(), scattering_accuracy); }
				retrieve_GPU_roughness_illumination(A);
				status << "Roughness scattering complete" << "\t" << currentDateTime() << endl;
			}

			if (roughness_thermal_model == 1) {
				freeup_GPU_shadow_arrays();
				freeup_GPU_globalillumination_arrays();
				initialise_GPU_temperature_arrays();
				initialise_GPU_globalheating_arrays();
			}

			if (minimum_thermal_inertia != 0.0f) {
				initialisetemperatures(A);
				status << "Initialised shape temperatures" << "\t" << currentDateTime() << endl;
			}

			if (roughness_thermal_model == 1) {
				if (minimum_thermal_inertia != 0.0f) {
					load_GPU_roughness_globalheating(A);
					status << "Roughness global heating data loaded" << "\t" << currentDateTime() << endl;
					GPU_initialise_temperatures(EPS, SIGMA);
					status << "Initialised roughness temperatures" << "\t" << currentDateTime() << endl;
				}
				else {
					load_GPU_roughness_globalheating(A);
					status << "Roughness global heating data loaded" << "\t" << currentDateTime() << endl;
					GPU_zero_temperatures();
					status << "Zeroed roughness temperatures" << "\t" << currentDateTime() << endl;
				}
			}

			A.reset_asteroid_status();
			A.reset_asteroid_yark_yorp();
			if (roughness_thermal_model == 1) { load_GPU_facet_properties(A); }

			// Thermal Model Loops

			rev = 0;
			while (A.rtnconverged() == 0 && rev < 301) {

				if (minimum_thermal_inertia == 0.0f && rev < 301) {
					instequilibrium(A);
				}
				else {
					inertiadifference(A);
				}

				if (rev > 5) { shape_converged(A); }

				status << "Number of shape facets converged during revolution " << rev << " is" << "\t" << A.rtnnumconverged() << "\t";
				status << "Shape facet convergence completion is" << "\t" << (100.0 * A.rtnnumconverged()) / ASTNUMFACE << "\t";
				status << "Time" << "\t" << currentDateTime() << endl;

				rev = rev + 1;

			}

			if (roughness_thermal_model == 1) {

				load_GPU_roughness_globalheating(A);
				status << "Roughness global heating data loaded" << "\t" << currentDateTime() << endl;

				rev = 0;
				while (A.rtnsfcconverged() == 0 && rev < 501) {

					if (minimum_thermal_inertia == 0.0f) {
						GPU_inst_equilibrium(EPS, SIGMA);
					}
					else {
						GPU_inertia_difference(EPS, SIGMA, TACC);
					}

					if (rev > 5) { roughness_converged(A); }

					status << "Number of roughness facets converged during revolution " << rev << " is" << "\t" << A.rtnsfcnumconverged() << "\t";
					status << "Roughness facet convergence completion is" << "\t" << (100.0 * A.rtnsfcnumconverged()) / ASTNUMFACE << "\t";
					status << "Time" << "\t" << currentDateTime() << endl;

					rev = rev + 1;

				}

			}

			if (roughness_thermal_model == 1) { retrieve_GPU_roughness_temperatures(A); }
			if (perform_yark_yorp == 1) { resolved_yark_yorp(A, S); }

			if (write_shape_temperature == 1) {
				filename1 = shape_filename + "_temperature_TA_" + pos.str() + "(" + obl.str() + ").txt";
				writeglobaltemperature(A, filename1);
			}
			if (write_roughness_temperature == 1 && roughness_thermal_model == 1) {
				filename1 = shape_filename + "_roughness_temperature_TA_" + pos.str() + "(" + obl.str() + ").dat";
				binary_writesurfacetemperature(S, A, filename1);
			}
			if (write_yarkovsky == 1) {
				filename1 = "yarkovsky_smooth_" + shape_filename + "_TA_" + pos.str() + "(" + obl.str() + ").txt";
				filename2 = "yarkovsky_rough_" + shape_filename + "_TA_" + pos.str() + "(" + obl.str() + ").txt";
				write_yarkovsky_results(A, filename1, filename2);
			}
			if (write_yorp == 1) {
				filename1 = "yorp_smooth_" + shape_filename + "_TA_" + pos.str() + "(" + obl.str() + ").txt";
				filename2 = "yorp_rough_" + shape_filename + "_TA_" + pos.str() + "(" + obl.str() + ").txt";
				write_yorp_results(A, filename1, filename2);
			}
			if (write_shape_shadow == 1) {
				filename1 = shape_filename + "_shadow_TA_" + pos.str() + "(" + obl.str() + ").txt";
				writeglobalshadow(A, filename1);
			}
			if (write_shape_illumination == 1) {
				filename1 = shape_filename + "_illumination_TA_" + pos.str() + "(" + obl.str() + ").txt";
				writeglobalillumination(A, filename1);
			}

			if (roughness_thermal_model == 1) {
				freeup_GPU_temperature_arrays();
				freeup_GPU_globalheating_arrays();
			}

		}

	}

};


// Main program loop


int main() {

	//	ifstream read("thermal_inertia_list2.txt");
	//	for (int i=0;i!=98;++i) {
	//		read >> thermal_inertia[i];
	//	}
	
	cout << "Setting cuda device: " << currentDateTime() << endl;
	set_cuda_device();
	cout << "Initialising asteroid and surface: " << currentDateTime() << endl;
	Asteroid A(ASTNUMVTX, ASTNUMFACE);
	Surface S(SFCNUMVTX, SFCNUMFACE);

	cout << "Adding surface roughness: " << currentDateTime() << endl;

	if (roughness_thermal_model == 1) { A.addsurfaceroughness(); }

	cout << "Initialising surface roughness: " << currentDateTime() << endl;

	if (roughness_thermal_model == 1) {
		S.readmodel();
		S.readselfheating();
		S.readphotonmomvectors();
		initialise_GPU_facet_properties_arrays();
		initialise_GPU_illumination_arrays();
	}

	cout << "Readubg asteroid model: " << currentDateTime() << endl;

	A.readmodel();
	//	A.readthermalinertia();
	if (load_shape_selfheating == 1) {
		A.readselfheating();
		A.readselfheating_vectors();
		A.readphotonmomvectors();
	}
	
	cout << "Staring single model call: " << currentDateTime() << endl;

	if (single_run == 1) { single_model(A, S); }
	if (orbital_run == 1) { orbital_model(A, S); }
	if (obliquity_run == 1) { obliquity_model(A, S); }
	cout << "Finishing: " << currentDateTime() << endl;
};
