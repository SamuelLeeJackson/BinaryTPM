#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <time.h>
using namespace std;

const double PI = 4.0*atan(1.0);									// Value of PI
double AU = 149.6e9;												// Astronomical unit in metres
double cvel = 299792458.0;											// Speed of light
double hcon = 6.626068963e-34;										// Planck constant
double kcon = 1.380650424e-23;										// Boltzman constant
double EPS = 0.95;													// Asteroid emissivity
double SIGMA = 5.67e-8;												// Stefan-Boltzmann constant

// Input Parameters

const int ASTNUMVTX = 2052;											// Number of shape model vertices
const int ASTNUMFACE = 4096;										// Number of shape model facets
const int PRINUMVTX = 1026;
const int NUMPRIREV = 9;
const int NUMSECREV = 2;
const int SFCNUMVTX = 61;											// Number of surface roughness model vertices
const int SFCNUMFACE = 100;											// Number of surface roughness model facets
const int NTAU = NUMPRIREV * 650;									// Number of model time steps

string shape_filename = "1996FG3_binary_lowres";					// Shape model input file name
double hel_polar0 = 1.160088;											// Heliocentric distance in AU of asteroid position
double hel_polar1 = ((-1.9521/180.0)*PI);								// Heliocentric ecliptic latitude in radians of asteroid position
double hel_polar2 = ((198.3142/180.0)*PI);								// Heliocentric ecliptic longitude in radians of asteroid position 
double obs_polar0 = 0.583905;											// Geocentric distance in AU of asteroid position
double obs_polar1 = ((-3.882500/180.0)*PI);								// Geocentric ecliptic latitude in radians of asteroid position
double obs_polar2 = ((257.492321/180.0)*PI);								// Geocentric ecliptic longitude in radians of asteroid position
double rot0 = ((-83.0/180.0)*PI);									// Helicentric ecliptic latitude in radians of asteroid rotation pole
double rot1 = ((266.0/180.0)*PI);									// Helicentric ecliptic longitude in radians of asteroid rotation pole

int observation_set_number = 4;										// Data set number (start at 0 and increment for each unique observation geometry/set)
int minimum_thermal_inertia = 0;									// Minimum thermal inertia value to run model for
int maximum_thermal_inertia = 400;									// Maximum thermal inertia value to run model for (set to minimum thermal inertia value if you want to run for a single and fixed thermal inertia value)
int thermal_inertia_step = 10;										// Thermal inertia step size

const int num_wavelengths = 3;											// Number of wavelengths observed
double obswavelength[num_wavelengths] = {4.6028,11.0984,22.6405};		// Wavelength in microns of observations

double wise_colour_corrections[3][467];


// Class AstVertex defining position of a single vertex of the global asteroid shape model


class AstVertex {

private:

	double _equpos[3];
	double _heliopos[NTAU][3];

public:

	AstVertex() {}

	friend std::istream& operator>> (std::istream&i, AstVertex &A) {
		return i >> A._equpos[0] >> A._equpos[1] >> A._equpos[2];
	}

	double &rtnequpos(int n) { return _equpos[n]; }
	double &rtnheliopos(int t, int n) { return _heliopos[t][n]; }

	~AstVertex() {}

};


// Class SubAstFacet to contain calculation data of surface roughness facets


class SubAstFacet {

private:

//	double _illumangle[NTAU];
//	double _totalillumination[NTAU];

	double _visible[NTAU];
	double _viewangle[NTAU];

	double _stemp[NTAU];

public:

	SubAstFacet() { 
		for (int step=0;step!=NTAU;++step) {
//			_totalillumination[step]=0.0;
			_visible[step]=1.0;
		}
	}

//	double &rtnillumangle(int t) { return _illumangle[t]; }
//	double &rtntotalillumination(int t) { return _totalillumination[t]; }

	double &rtnvisible(int t) { return _visible[t]; }
	double &rtnviewangle(int t) { return _viewangle[t]; }

	double &rtnstemp(int t) { return _stemp[t]; }

	~SubAstFacet() {}

};


// Class AstFacet defining surface properties of a single surface facet of the global asteroid shape model


class AstFacet {

private:

	int _astvertices[3];
//	double _illumvector[NTAU][3];
	double _viewvector[NTAU][3];
	double _normal[NTAU][3];
	double _area;
//	double _illumangle[NTAU];
	double _viewangle[NTAU];

//	int _shadow[NTAU];
	int _visible[NTAU];

	double _albedo;
	int _surfaceroughness;

	double _stemp[NTAU];

	double _smooth_flux[num_wavelengths];
	double _rough_flux[num_wavelengths];

	SubAstFacet *SF;

public:

	SubAstFacet*beginsf() const { return SF; }
	SubAstFacet*endsf() const { return SF+SFCNUMFACE; }

	AstFacet() {
		_surfaceroughness=0;
		_albedo=0.1;
		for (int step=0;step!=NTAU;++step) {
//			_shadow[step]=1;
			_visible[step]=1;
		}
		for (int i=0;i!=num_wavelengths;++i) {
			_smooth_flux[i]=0.0;
			_rough_flux[i]=0.0;
		}
	}

	void createsubfacets() {
		SF=new SubAstFacet[SFCNUMFACE];
	}

	friend std::istream& operator>> (std::istream&i, AstFacet&A) {
		return i >> A._astvertices[0] >> A._astvertices[1] >> A._astvertices[2];
	}

	int &rtnastvertices(int n) { return _astvertices[n]; }
//	double &rtnillumvector(int t, int n) { return _illumvector[t][n]; }
	double &rtnviewvector(int t, int n) { return _viewvector[t][n]; }
	double &rtnnormal(int t, int n) { return _normal[t][n]; }
	double &rtnarea() { return _area; }
//	double &rtnillumangle(int t) { return _illumangle[t]; }
	double &rtnviewangle(int t) { return _viewangle[t]; }

	double &rtnalbedo() { return _albedo; }
	int &rtnsurfaceroughness() { return _surfaceroughness; }

//	int &rtnshadow(int t) { return _shadow[t]; }
	int &rtnvisible(int t) { return _visible[t]; }

	double &rtnstemp(int t) { return _stemp[t]; }

	double &rtnsmooth_flux(int n) { return _smooth_flux[n]; }
	double &rtnrough_flux(int n) { return _rough_flux[n]; }

	~AstFacet() {}

};


// Class Asteroid to control the asteroid global shape model


class Asteroid {

private:

	int NV, NF;
	double _polcentre[3];
	double _carcentre[3];
	double _rotvector[3];
	double _rotperiod;
	double _solar;
	double _totalarea;

	int _surfaceroughness;

	AstVertex *V;
	AstFacet *F;

public:

	AstVertex*beginv() const { return V; }
	AstVertex*endv() const { return V+NV; }

	AstFacet*beginf() const { return F; }
	AstFacet*endf() const { return F+NF; }
	
	Asteroid(int n, int m) : NV(n), NF(m) {
		V=new AstVertex[n];
		F=new AstFacet[m];
		_polcentre[0]=hel_polar0, _polcentre[1]=hel_polar1, _polcentre[2]=hel_polar2;
		_rotvector[0]=rot0, _rotvector[1]=rot1;
		_surfaceroughness=0;
	}

	void readglobalmodel(string smooth_filename) {
		ifstream read(smooth_filename.c_str());
		double temp;
		for (AstVertex *i=V;i!=V+NV;++i) {
			read >> temp;
			read >> temp;
			read >> temp;
		}
		for (AstFacet *j=F;j!=F+NF;++j) {
			read >> temp;
			read >> temp;
			read >> temp;
		}
		for (AstFacet *j=F;j!=F+NF;++j) {
			for (int step=0;step!=NTAU;++step) {
				read >> temp;
				j->rtnstemp(step)=temp;
			}
		}

		if ( read.fail()==false ) { cout << "Ok"; } else { cout << "Bad"; } 

	}

	void readvisibility(string visibility_filename) {
		ifstream read(visibility_filename.c_str());
		int temp;
		for (AstVertex *i=V;i!=V+NV;++i) {
			read >> *i;
		}
		for (AstFacet *j=F;j!=F+NF;++j) {
			read >> *j;
		}
		for (AstFacet *j=F;j!=F+NF;++j) {
			for (int step=0;step!=NTAU;++step) {
				read >> temp;
				j->rtnvisible(step)=temp;
			}
		}

		if ( read.fail()==false ) { cout << "Ok"; } else { cout << "Bad"; } 

	}

	void addsurfaceroughness() {
		for (AstFacet *i=F;i!=F+NF;++i) {
			i->createsubfacets();
			i->rtnsurfaceroughness()=1;
		}
		_surfaceroughness=1;
	}

	int &rtnnumvtx() { return NV; }
	int &rtnnumface() { return NF; }
	double &rtnpolcentre(int n) { return _polcentre[n]; }
	double &rtncarcentre(int n) { return _carcentre[n]; }
	double &rtnrotvector(int n) { return _rotvector[n]; }
	double &rtnrotperiod() { return _rotperiod; }
	double &rtnsolar() { return _solar; }
	double &rtntotalarea() { return _totalarea; }
	int &rtnsurfaceroughness() { return _surfaceroughness; }

};


// Class SfcVertex defining position of a single vertex of the surface roughness shape model


class SfcVertex {

private:

	double _pos[3];

public:

	SfcVertex() {}

	friend std::istream& operator>> (std::istream&i, SfcVertex&A) {
		return i >> A._pos[0] >> A._pos[1] >> A._pos[2];
	}

	double &rtnpos(int n) { return _pos[n]; }

	~SfcVertex() {}

};


// Class SfcFacet defining surface properties of a single surface facet of the surface roughness shape model


class SfcFacet {

private:

	int _sfcvertices[3];
	double _midpos[3];
	double _normal[3];
	double _area;

public:

	SfcFacet() {}

	friend std::istream& operator>> (std::istream&i, SfcFacet&A) {
		return i >> A._sfcvertices[0] >> A._sfcvertices[1] >> A._sfcvertices[2];
	}

	int &rtnsfcvertices(int n) { return _sfcvertices[n]; }
	double &rtnmidpos(int n) { return _midpos[n]; }
	double &rtnnormal(int n) { return _normal[n]; }
	double &rtnarea() { return _area;}

	~SfcFacet() {}

};


// Class Surface to control the surface roughness shape model


class Surface {

private:

	int NV, NF;
	double _totalarea;
	SfcVertex *SV;
	SfcFacet *SF;

public:

	Surface(int n, int m): NV(n), NF(m) {
		SV=new SfcVertex[n];
		SF=new SfcFacet[m];
	}

	SfcVertex*beginv() const { return SV; }
	SfcVertex*endv() const { return SV+NV; }

	SfcFacet*beginf() const { return SF; }
	SfcFacet*endf() const { return SF+NF; }

	int &rtnnumvtx() { return NV; }
	int &rtnnumface() { return NF; }
	double &rtntotalarea() { return _totalarea; }

};


// Binary load surface temperatures function


void binary_loadsurfacetemperatures(Surface &S, Asteroid &A, string rough_filename) {

	ifstream read(rough_filename.c_str(),ios::in|ios::binary);

	float pos[3];
	int vert[3];
	float val[NTAU];

	for (SfcVertex *B=S.beginv();B!=S.endv();++B) {
		read.read((char *)&pos,sizeof(float)*3);
	}

	for (SfcFacet *C=S.beginf();C!=S.endf();++C) {
		read.read((char *)&vert,sizeof(int)*3);
	}

	for (AstFacet *C=A.beginf();C!=A.endf();++C) {
		for (SubAstFacet *D=C->beginsf();D!=C->endsf();++D) {
			read.read((char *)&val,sizeof(float)*NTAU);	
			for (int step=0;step!=NTAU;++step) {
				D->rtnstemp(step)=val[step];
			}
		}
	}

};


// Binary load surface visibility function


void binary_loadsurfacevisibility(Surface &S, Asteroid &A, string rough_visibility_filename) {

	ifstream read(rough_visibility_filename.c_str(),ios::in|ios::binary);

	float pos[3];
	int vert[3];
	float val[NTAU];

	for (SfcVertex *B=S.beginv();B!=S.endv();++B) {
		read.read((char *)&pos,sizeof(float)*3);
		B->rtnpos(0)=pos[0], B->rtnpos(1)=pos[1], B->rtnpos(2)=pos[2];
	}

	for (SfcFacet *C=S.beginf();C!=S.endf();++C) {
		read.read((char *)&vert,sizeof(int)*3);
		C->rtnsfcvertices(0)=vert[0], C->rtnsfcvertices(1)=vert[1], C->rtnsfcvertices(2)=vert[2];
	}

	for (AstFacet *C=A.beginf();C!=A.endf();++C) {
		for (SubAstFacet *D=C->beginsf();D!=C->endsf();++D) {
			read.read((char *)&val,sizeof(float)*NTAU);	
			for (int step=0;step!=NTAU;++step) {
				D->rtnvisible(step)=val[step];
			}
		}
	}

};


// Class Observer defining properties of observer


class Observer {

private:

	double _polcentre[3];
	double _carcentre[3];

public:

	Observer() { _polcentre[0]=obs_polar0, _polcentre[1]=obs_polar1, _polcentre[2]=obs_polar2; }

	double &rtnpolcentre(int n) { return _polcentre[n]; }
	double &rtncarcentre(int n) { return _carcentre[n]; }

	~Observer() {}

};


// Determine surface roughness geometry


void surfacegeometry(Surface &S, Asteroid &A, Observer &O) {

	for (SfcFacet *C=S.beginf();C!=S.endf();++C) {

		C->rtnmidpos(0)=((S.beginv()+C->rtnsfcvertices(0))->rtnpos(0)+(S.beginv()+C->rtnsfcvertices(1))->rtnpos(0)+(S.beginv()+C->rtnsfcvertices(2))->rtnpos(0))/3.0;
		C->rtnmidpos(1)=((S.beginv()+C->rtnsfcvertices(0))->rtnpos(1)+(S.beginv()+C->rtnsfcvertices(1))->rtnpos(1)+(S.beginv()+C->rtnsfcvertices(2))->rtnpos(1))/3.0;
		C->rtnmidpos(2)=((S.beginv()+C->rtnsfcvertices(0))->rtnpos(2)+(S.beginv()+C->rtnsfcvertices(1))->rtnpos(2)+(S.beginv()+C->rtnsfcvertices(2))->rtnpos(2))/3.0;
	
	}

	double normcalc[2][3], n;

	S.rtntotalarea()=0.0;

	for (SfcFacet *C=S.beginf();C!=S.endf();++C) {

		normcalc[0][0]=(S.beginv()+C->rtnsfcvertices(1))->rtnpos(0)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(0);
		normcalc[0][1]=(S.beginv()+C->rtnsfcvertices(1))->rtnpos(1)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(1);
		normcalc[0][2]=(S.beginv()+C->rtnsfcvertices(1))->rtnpos(2)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(2);

		normcalc[1][0]=(S.beginv()+C->rtnsfcvertices(2))->rtnpos(0)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(0);
		normcalc[1][1]=(S.beginv()+C->rtnsfcvertices(2))->rtnpos(1)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(1);
		normcalc[1][2]=(S.beginv()+C->rtnsfcvertices(2))->rtnpos(2)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(2);

		C->rtnnormal(0)=(normcalc[0][1]*normcalc[1][2])-(normcalc[0][2]*normcalc[1][1]);
		C->rtnnormal(1)=(normcalc[0][2]*normcalc[1][0])-(normcalc[0][0]*normcalc[1][2]);
		C->rtnnormal(2)=(normcalc[0][0]*normcalc[1][1])-(normcalc[0][1]*normcalc[1][0]);

		n=sqrt((C->rtnnormal(0)*C->rtnnormal(0))+(C->rtnnormal(1)*C->rtnnormal(1))+(C->rtnnormal(2)*C->rtnnormal(2)));
		C->rtnarea()=n/2.0;

		C->rtnnormal(0)=C->rtnnormal(0)/n;
		C->rtnnormal(1)=C->rtnnormal(1)/n;
		C->rtnnormal(2)=C->rtnnormal(2)/n;

		S.rtntotalarea()=S.rtntotalarea()+(C->rtnarea()*C->rtnnormal(2));

	}

	double _t,_o;

	for (int step=0;step!=NTAU;++step) {

		for (AstFacet *C=A.beginf();C!=A.endf();++C) {

//			_t=sqrt((C->rtnillumvector(step,0)*C->rtnillumvector(step,0))+(C->rtnillumvector(step,1)*C->rtnillumvector(step,1))+(C->rtnillumvector(step,2)*C->rtnillumvector(step,2)));
			_o=sqrt((C->rtnviewvector(step,0)*C->rtnviewvector(step,0))+(C->rtnviewvector(step,1)*C->rtnviewvector(step,1))+(C->rtnviewvector(step,2)*C->rtnviewvector(step,2)));

			for (int i=0;i!=SFCNUMFACE;++i) {

//				(C->beginsf()+i)->rtnillumangle(step)=acos(((C->rtnillumvector(step,0)*(S.beginf()+i)->rtnnormal(0))+(C->rtnillumvector(step,1)*(S.beginf()+i)->rtnnormal(1))+(C->rtnillumvector(step,2)*(S.beginf()+i)->rtnnormal(2)))/_t);
				(C->beginsf()+i)->rtnviewangle(step)=acos(((C->rtnviewvector(step,0)*(S.beginf()+i)->rtnnormal(0))+(C->rtnviewvector(step,1)*(S.beginf()+i)->rtnnormal(1))+(C->rtnviewvector(step,2)*(S.beginf()+i)->rtnnormal(2)))/_o);

			}

		}

	}

};


// Determine asteroid global shape model geometry


void asteroidgeometry(Asteroid &A, Observer &O) {

	A.rtnsolar()=1367.0/(A.rtnpolcentre(0)*A.rtnpolcentre(0));

	A.rtncarcentre(0)=AU*A.rtnpolcentre(0)*cos(A.rtnpolcentre(1))*cos(A.rtnpolcentre(2));
	A.rtncarcentre(1)=AU*A.rtnpolcentre(0)*cos(A.rtnpolcentre(1))*sin(A.rtnpolcentre(2));
	A.rtncarcentre(2)=AU*A.rtnpolcentre(0)*sin(A.rtnpolcentre(1));

	O.rtncarcentre(0)=A.rtncarcentre(0)-(AU*O.rtnpolcentre(0)*cos(O.rtnpolcentre(1))*cos(O.rtnpolcentre(2)));
	O.rtncarcentre(1)=A.rtncarcentre(1)-(AU*O.rtnpolcentre(0)*cos(O.rtnpolcentre(1))*sin(O.rtnpolcentre(2)));
	O.rtncarcentre(2)=A.rtncarcentre(2)-(AU*O.rtnpolcentre(0)*sin(O.rtnpolcentre(1)));

	double tmatrix[3][3], gamma;

	gamma=acos(cos(A.rtnrotvector(0))*sin(fabs(A.rtnrotvector(1))));

	tmatrix[0][0]=sin(A.rtnrotvector(0))/sin(gamma);
	tmatrix[0][1]=-1.0*cos(A.rtnrotvector(0))*cos(A.rtnrotvector(0))*sin(A.rtnrotvector(1))*cos(A.rtnrotvector(1))/sin(gamma);
	tmatrix[0][2]=cos(A.rtnrotvector(0))*cos(A.rtnrotvector(1));
	tmatrix[1][0]=0.0;
	tmatrix[1][1]=((sin(A.rtnrotvector(0))*sin(A.rtnrotvector(0)))+(cos(A.rtnrotvector(0))*cos(A.rtnrotvector(0))*cos(A.rtnrotvector(1))*cos(A.rtnrotvector(1))))/sin(gamma);
	tmatrix[1][2]=cos(A.rtnrotvector(0))*sin(A.rtnrotvector(1));
	tmatrix[2][0]=-1.0*cos(A.rtnrotvector(0))*cos(A.rtnrotvector(1))/sin(gamma);
	tmatrix[2][1]=-1.0*sin(A.rtnrotvector(0))*cos(A.rtnrotvector(0))*sin(A.rtnrotvector(1))/sin(gamma);
	tmatrix[2][2]=sin(A.rtnrotvector(0));

	double rotangle, rotequpos[3], eclipos[3], heliomidpos[3], obsvector[3];
	double normcalc[2][3], y[3], _n, _t, _x, _o;
	double imatrix[3][3]; 

	for (int step=0;step!=NTAU;++step) {

	//	rotangle=2.0*PI*step/NTAU;

	//	for (AstVertex *B=A.beginv();B!=A.endv();++B) {

	//		rotequpos[0]=(B->rtnequpos(0)*cos(rotangle))-(B->rtnequpos(1)*sin(rotangle));
	//		rotequpos[1]=(B->rtnequpos(0)*sin(rotangle))+(B->rtnequpos(1)*cos(rotangle));
	//		rotequpos[2]=B->rtnequpos(2);

	//		eclipos[0]=(tmatrix[0][0]*rotequpos[0])+(tmatrix[0][1]*rotequpos[1])+(tmatrix[0][2]*rotequpos[2]);
	//		eclipos[1]=(tmatrix[1][0]*rotequpos[0])+(tmatrix[1][1]*rotequpos[1])+(tmatrix[1][2]*rotequpos[2]);
	//		eclipos[2]=(tmatrix[2][0]*rotequpos[0])+(tmatrix[2][1]*rotequpos[1])+(tmatrix[2][2]*rotequpos[2]);

	//		B->rtnheliopos(step,0)=A.rtncarcentre(0)+eclipos[0];
	//		B->rtnheliopos(step,1)=A.rtncarcentre(1)+eclipos[1];
	//		B->rtnheliopos(step,2)=A.rtncarcentre(2)+eclipos[2];

	//	}

		rotangle = 2.0*PI*step*NUMPRIREV / NTAU;

		for (AstVertex *B = A.beginv(); B != A.beginv() + PRINUMVTX; ++B) {

			rotequpos[0] = (B->rtnequpos(0)*cos(rotangle)) - (B->rtnequpos(1)*sin(rotangle));
			rotequpos[1] = (B->rtnequpos(0)*sin(rotangle)) + (B->rtnequpos(1)*cos(rotangle));
			rotequpos[2] = B->rtnequpos(2);

			eclipos[0] = (tmatrix[0][0] * rotequpos[0]) + (tmatrix[0][1] * rotequpos[1]) + (tmatrix[0][2] * rotequpos[2]);
			eclipos[1] = (tmatrix[1][0] * rotequpos[0]) + (tmatrix[1][1] * rotequpos[1]) + (tmatrix[1][2] * rotequpos[2]);
			eclipos[2] = (tmatrix[2][0] * rotequpos[0]) + (tmatrix[2][1] * rotequpos[1]) + (tmatrix[2][2] * rotequpos[2]);

			B->rtnheliopos(step, 0) = A.rtncarcentre(0) + eclipos[0];
			B->rtnheliopos(step, 1) = A.rtncarcentre(1) + eclipos[1];
			B->rtnheliopos(step, 2) = A.rtncarcentre(2) + eclipos[2];

		}

		rotangle = 2.0*PI*step*NUMSECREV / NTAU;

		for (AstVertex *B = A.beginv() + PRINUMVTX; B != A.endv(); ++B) {

			rotequpos[0] = (B->rtnequpos(0)*cos(rotangle)) - (B->rtnequpos(1)*sin(rotangle));
			rotequpos[1] = (B->rtnequpos(0)*sin(rotangle)) + (B->rtnequpos(1)*cos(rotangle));
			rotequpos[2] = B->rtnequpos(2);

			eclipos[0] = (tmatrix[0][0] * rotequpos[0]) + (tmatrix[0][1] * rotequpos[1]) + (tmatrix[0][2] * rotequpos[2]);
			eclipos[1] = (tmatrix[1][0] * rotequpos[0]) + (tmatrix[1][1] * rotequpos[1]) + (tmatrix[1][2] * rotequpos[2]);
			eclipos[2] = (tmatrix[2][0] * rotequpos[0]) + (tmatrix[2][1] * rotequpos[1]) + (tmatrix[2][2] * rotequpos[2]);

			B->rtnheliopos(step, 0) = A.rtncarcentre(0) + eclipos[0];
			B->rtnheliopos(step, 1) = A.rtncarcentre(1) + eclipos[1];
			B->rtnheliopos(step, 2) = A.rtncarcentre(2) + eclipos[2];

		}

		for (AstFacet *C=A.beginf();C!=A.endf();++C) {

			heliomidpos[0]=(((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,0))+((A.beginv()+C->rtnastvertices(1))->rtnheliopos(step,0))+((A.beginv()+C->rtnastvertices(2))->rtnheliopos(step,0)))/3.0;
			heliomidpos[1]=(((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,1))+((A.beginv()+C->rtnastvertices(1))->rtnheliopos(step,1))+((A.beginv()+C->rtnastvertices(2))->rtnheliopos(step,1)))/3.0;
			heliomidpos[2]=(((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,2))+((A.beginv()+C->rtnastvertices(1))->rtnheliopos(step,2))+((A.beginv()+C->rtnastvertices(2))->rtnheliopos(step,2)))/3.0;

			_t=sqrt((heliomidpos[0]*heliomidpos[0])+(heliomidpos[1]*heliomidpos[1])+(heliomidpos[2]*heliomidpos[2]));
			
			obsvector[0]=heliomidpos[0]-O.rtncarcentre(0);
			obsvector[1]=heliomidpos[1]-O.rtncarcentre(1);
			obsvector[2]=heliomidpos[2]-O.rtncarcentre(2);

			_o=sqrt((obsvector[0]*obsvector[0])+(obsvector[1]*obsvector[1])+(obsvector[2]*obsvector[2]));

			normcalc[0][0]=((A.beginv()+C->rtnastvertices(1))->rtnheliopos(step,0))-((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,0));
			normcalc[0][1]=((A.beginv()+C->rtnastvertices(1))->rtnheliopos(step,1))-((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,1));
			normcalc[0][2]=((A.beginv()+C->rtnastvertices(1))->rtnheliopos(step,2))-((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,2));

			normcalc[1][0]=((A.beginv()+C->rtnastvertices(2))->rtnheliopos(step,0))-((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,0));
			normcalc[1][1]=((A.beginv()+C->rtnastvertices(2))->rtnheliopos(step,1))-((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,1));	
			normcalc[1][2]=((A.beginv()+C->rtnastvertices(2))->rtnheliopos(step,2))-((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,2));

			C->rtnnormal(step,0)=(normcalc[0][1]*normcalc[1][2])-(normcalc[0][2]*normcalc[1][1]);
			C->rtnnormal(step,1)=(normcalc[0][2]*normcalc[1][0])-(normcalc[0][0]*normcalc[1][2]);
			C->rtnnormal(step,2)=(normcalc[0][0]*normcalc[1][1])-(normcalc[0][1]*normcalc[1][0]);

			_n=sqrt((C->rtnnormal(step,0)*C->rtnnormal(step,0))+(C->rtnnormal(step,1)*C->rtnnormal(step,1))+(C->rtnnormal(step,2)*C->rtnnormal(step,2)));
			C->rtnarea()=_n/2.0;

			C->rtnnormal(step,0)=C->rtnnormal(step,0)/_n;
			C->rtnnormal(step,1)=C->rtnnormal(step,1)/_n;
			C->rtnnormal(step,2)=C->rtnnormal(step,2)/_n;

			_x=sqrt((normcalc[0][0]*normcalc[0][0])+(normcalc[0][1]*normcalc[0][1])+(normcalc[0][2]*normcalc[0][2]));

			normcalc[0][0]=normcalc[0][0]/_x;
			normcalc[0][1]=normcalc[0][1]/_x;
			normcalc[0][2]=normcalc[0][2]/_x;

			y[0]=(C->rtnnormal(step,1)*normcalc[0][2])-(C->rtnnormal(step,2)*normcalc[0][1]);
			y[1]=(C->rtnnormal(step,2)*normcalc[0][0])-(C->rtnnormal(step,0)*normcalc[0][2]);
			y[2]=(C->rtnnormal(step,0)*normcalc[0][1])-(C->rtnnormal(step,1)*normcalc[0][0]);

			imatrix[0][0]=normcalc[0][0];
			imatrix[0][1]=normcalc[0][1];
			imatrix[0][2]=normcalc[0][2];
			imatrix[1][0]=y[0];
			imatrix[1][1]=y[1];
			imatrix[1][2]=y[2];
			imatrix[2][0]=C->rtnnormal(step,0);
			imatrix[2][1]=C->rtnnormal(step,1);
			imatrix[2][2]=C->rtnnormal(step,2);

//			C->rtnillumvector(step,0)=-1.0*((imatrix[0][0]*heliomidpos[0])+(imatrix[0][1]*heliomidpos[1])+(imatrix[0][2]*heliomidpos[2]));
//			C->rtnillumvector(step,1)=-1.0*((imatrix[1][0]*heliomidpos[0])+(imatrix[1][1]*heliomidpos[1])+(imatrix[1][2]*heliomidpos[2]));
//			C->rtnillumvector(step,2)=-1.0*((imatrix[2][0]*heliomidpos[0])+(imatrix[2][1]*heliomidpos[1])+(imatrix[2][2]*heliomidpos[2]));

//			C->rtnillumangle(step)=acos(C->rtnillumvector(step,2)/_t);
//			if (C->rtnillumangle(step)>(PI/2.0)) {C->rtnshadow(step)=0;}

			C->rtnviewvector(step,0)=-1.0*((imatrix[0][0]*obsvector[0])+(imatrix[0][1]*obsvector[1])+(imatrix[0][2]*obsvector[2]));
			C->rtnviewvector(step,1)=-1.0*((imatrix[1][0]*obsvector[0])+(imatrix[1][1]*obsvector[1])+(imatrix[1][2]*obsvector[2]));
			C->rtnviewvector(step,2)=-1.0*((imatrix[2][0]*obsvector[0])+(imatrix[2][1]*obsvector[1])+(imatrix[2][2]*obsvector[2]));

			C->rtnviewangle(step)=acos(C->rtnviewvector(step,2)/_o);
			if (C->rtnviewangle(step)>(PI/2.0)) {C->rtnvisible(step)=0;}

		}

	}

	A.rtntotalarea()=0.0;

	for (AstFacet *C=A.beginf();C!=A.endf();++C) {
		A.rtntotalarea()=A.rtntotalarea()+C->rtnarea();
	}

};


// Code to observe asteroid


void observe(Surface &S, Asteroid &A, string smooth_outname, string rough_outname) {

	ofstream smoothspectrum(smooth_outname.c_str());
	ofstream roughspectrum(rough_outname.c_str());

	double intensity;
	double multiplier;
	double wavelength;
	double _d;
	double astsmoothfluxtotal[num_wavelengths];
	double astroughfluxtotal[num_wavelengths];
	double sfcfluxtotal[num_wavelengths];
	double phase;

	double chisqr, chisqr_smooth, chisqr_rough; 
	int best_T_smooth, best_T_rough;
	double b2, b3, b4, avgb;
	double m2, m3, m4, avgm;

	for (int step=0;step!=NTAU;++step) {
//	for (int step=0;step!=1;++step) {

		phase=1.0*step/NTAU;

		for (int i=0;i!=num_wavelengths;++i) {
			astsmoothfluxtotal[i]=0.0;
			astroughfluxtotal[i]=0.0;
		}

		for (AstFacet *B=A.beginf();B!=A.endf();++B) {

			if ( B->rtnvisible(step)!=0 ) {

				_d=sqrt((B->rtnviewvector(step,0)*B->rtnviewvector(step,0))+(B->rtnviewvector(step,1)*B->rtnviewvector(step,1))+(B->rtnviewvector(step,2)*B->rtnviewvector(step,2)));

				for (int i=0;i!=num_wavelengths;++i) {
					sfcfluxtotal[i]=0.0;
				}

				if ( A.rtnsurfaceroughness()==1 ) {

					for (int j=0;j!=SFCNUMFACE;++j) {

						if ( (B->beginsf()+j)->rtnvisible(step)>0.0 ) {

							for (int i=0;i!=num_wavelengths;++i) {
								wavelength=obswavelength[i]*1.0e-6;
								intensity=((2.0*hcon*cvel*cvel)/pow(wavelength,5.0))*(1.0/(exp((hcon*cvel)/(wavelength*kcon*(B->beginsf()+j)->rtnstemp(step)))-1.0));
								multiplier=(EPS*(S.beginf()+j)->rtnarea()*(B->beginsf()+j)->rtnvisible(step)*cos((B->beginsf()+j)->rtnviewangle(step)))/(_d*_d);
								if (intensity!=intensity) { intensity=0.0; }
								if (multiplier!=multiplier) { multiplier=0.0; }
								sfcfluxtotal[i]=sfcfluxtotal[i]+(intensity*multiplier);
							}

						}

					}

					for (int i=0;i!=num_wavelengths;++i) {
						astroughfluxtotal[i]=astroughfluxtotal[i]+(sfcfluxtotal[i]*B->rtnarea()/S.rtntotalarea());
						B->rtnrough_flux(i)=B->rtnrough_flux(i)+((sfcfluxtotal[i]*B->rtnarea())/(NTAU*S.rtntotalarea()));
					}

				}

				for (int i=0;i!=num_wavelengths;++i) {
					wavelength=obswavelength[i]*1.0e-6;
					intensity=((2.0*hcon*cvel*cvel)/pow(wavelength,5.0))*(1.0/(exp((hcon*cvel)/(wavelength*kcon*B->rtnstemp(step)))-1.0));
					multiplier=(EPS*B->rtnarea()*B->rtnvisible(step)*cos(B->rtnviewangle(step)))/(_d*_d);
					if (intensity!=intensity) { intensity=0.0; }
					if (multiplier!=multiplier) { multiplier=0.0; }
					astsmoothfluxtotal[i]=astsmoothfluxtotal[i]+(intensity*multiplier);
					B->rtnsmooth_flux(i)=B->rtnsmooth_flux(i)+(intensity*multiplier/NTAU);
				}

			}

		}

		chisqr_smooth=100.0, chisqr_rough=100.0;
		best_T_smooth=100, best_T_rough=100;

		for (int T=100;T!=567;++T) {

			wavelength=obswavelength[0]*1.0e-6;
			b2=((2.0*hcon*cvel*cvel)/pow(wavelength,5.0))*(1.0/(exp((hcon*cvel)/(wavelength*kcon*(double)T))-1.0));
			wavelength=obswavelength[1]*1.0e-6;
			b3=((2.0*hcon*cvel*cvel)/pow(wavelength,5.0))*(1.0/(exp((hcon*cvel)/(wavelength*kcon*(double)T))-1.0));
			wavelength=obswavelength[2]*1.0e-6;
			b4=((2.0*hcon*cvel*cvel)/pow(wavelength,5.0))*(1.0/(exp((hcon*cvel)/(wavelength*kcon*(double)T))-1.0));

			avgb=(b2+b3+b4)/3.0;

			b2=b2/avgb;
			b3=b3/avgb;
			b4=b4/avgb;

			avgm=(astsmoothfluxtotal[0]+astsmoothfluxtotal[1]+astsmoothfluxtotal[2])/3.0;

			m2=astsmoothfluxtotal[0]/avgm;
			m3=astsmoothfluxtotal[1]/avgm;
			m4=astsmoothfluxtotal[2]/avgm;

			chisqr=pow((b2-m2)/m2,2.0)+((b3-m3)/m3,2.0)+((b4-m4)/m4,2.0);
			if (chisqr<chisqr_smooth) {
				chisqr_smooth=chisqr;
				best_T_smooth=T;
			}

			avgm=(astroughfluxtotal[0]+astroughfluxtotal[1]+astroughfluxtotal[2])/3.0;

			m2=astroughfluxtotal[0]/avgm;
			m3=astroughfluxtotal[1]/avgm;
			m4=astroughfluxtotal[2]/avgm;

			chisqr=pow((b2-m2)/m2,2.0)+pow((b3-m3)/m3,2.0)+pow((b4-m4)/m4,2.0);
			if (chisqr<chisqr_rough) {
				chisqr_rough=chisqr;
				best_T_rough=T;
			}

		}

		for (int i=0;i!=num_wavelengths;++i) {
			smoothspectrum << astsmoothfluxtotal[i]*wise_colour_corrections[i][best_T_smooth-100] << "\t";
			roughspectrum << astroughfluxtotal[i]*wise_colour_corrections[i][best_T_rough-100] << "\t";
		//	smoothspectrum << astsmoothfluxtotal[i] << "\t" << astsmoothfluxtotal[i]*wise_colour_corrections[i][best_T_smooth-100] << "\t";
		//	roughspectrum << astroughfluxtotal[i] << "\t" << astroughfluxtotal[i]*wise_colour_corrections[i][best_T_rough-100] << "\t";
		}

		smoothspectrum << endl;
		roughspectrum << endl;

	}

};

// Function returning current time

string get_time() {
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time (&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer,sizeof(buffer),"%Y-%m-%d %H:%M:%S",timeinfo);
	std::string str(buffer);
	return str;
}

// Main program loop


int main() {
	
	stringstream inert;
#ifdef T_INERTIA
	inert << T_INERTIA;
#else
	inert << 0;
#endif
	string status_filename=shape_filename+"_observe_status("+inert.str()+").txt";
	ofstream status(status_filename.c_str());
	string time;
	time = get_time();
	status << "Start time" << "\t" << time << endl;
	ifstream read_corrections("wise_colour_corrections.txt");
	for (int i=0;i!=467;++i) {
		read_corrections >> wise_colour_corrections[0][i];
		read_corrections >> wise_colour_corrections[1][i];
		read_corrections >> wise_colour_corrections[2][i];
	}
	status << "Read the WISE colour corrections" << endl;

	Asteroid A(ASTNUMVTX, ASTNUMFACE);
	status << get_time() << ":\tInititialised asteroid object" << endl;

	Surface S(SFCNUMVTX, SFCNUMFACE);
	A.addsurfaceroughness();
	status << get_time() << ":\tInititialised surfaces" << endl;

	Observer O;

	string visibility_filename, rough_visibility_filename;
	string smooth_filename, rough_filename;
	string smooth_outname, rough_outname;

	visibility_filename=shape_filename+"_visibility.txt";
	rough_visibility_filename=shape_filename+"_roughness_visibility.dat";

	A.readvisibility(visibility_filename);
	status << get_time() << ":\tRead shape visibility" << endl;
	binary_loadsurfacevisibility(S,A,rough_visibility_filename);
	status << get_time() << ":\tRead roughness visibility" << endl;

	asteroidgeometry(A,O);
	surfacegeometry(S,A,O);
	status << get_time() << ":\tCalculated geometries" << endl;

	stringstream lcurve;
	lcurve << observation_set_number;


	status << get_time() << ":\tRunning observe for TI = "+inert.str() << endl;
	smooth_filename=shape_filename+"_temperature("+inert.str()+").txt";
	rough_filename=shape_filename+"_roughness_temperature("+inert.str()+").dat";
	smooth_outname=shape_filename+"_smooth_flux("+inert.str()+"_"+lcurve.str()+").txt";
	rough_outname=shape_filename+"_rough_flux("+inert.str()+"_"+lcurve.str()+").txt";
//	smooth_outname="sevastopol_smooth_flux("+inert.str()+"_"+lcurve.str()+"_1).txt";
//	rough_outname="sevastopol_rough_flux("+inert.str()+"_"+lcurve.str()+"_1).txt";
	
	A.readglobalmodel(smooth_filename);
	status << get_time() << ":\tRead smooth surface temperatures" << endl;
	binary_loadsurfacetemperatures(S,A,rough_filename);
	status << get_time() << ":\tLoaded rough surface temperatures" << endl;
	observe(S,A,smooth_outname,rough_outname);
	status << get_time() << ":\tCalculated fluxes" << endl;

};