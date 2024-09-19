#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <time.h>
#include <cstdlib>
using namespace std;

// Constants

const double PI = 4.0*atan(1.0);														// Value of PI
const double AU = 149.6e6;																// Astronomical unit in kilometres

// Input Parameters

const int ASTNUMVTX = 3204;																// Number of shape model vertices
const int ASTNUMFACE = 6400;															// Number of shape model facets
string shape_filename = "doppler_binary";															// Shape model input file name
int load_shape_visibility = 0;															// Load shape visibility (= 1) or not (= 0)
const int SFCNUMVTX = 61;																// Number of surface roughness model vertices
const int SFCNUMFACE = 100;																// Number of surface roughness model facets
string roughness_filename = "new_crater2";												// Surface roughness model input file name
double rot0 = ((60.0 / 180.0) * PI);										// Heliocentric ecliptic latitude in radians of asteroid rotation pole
double rot1 = ((221.0 / 180.0) * PI);									// Heliocentric ecliptic longitude in radians of asteroid rotation pole
double ob0 = 1.631629;															// Geocentric distance in AU of asteroid
double ob1 = ((13.826432 /180.0)*PI);													// Geocentric ecliptic latitude in radians of asteroid position
double ob2 = ((80.190885 /180.0)*PI);													// Geocentric ecliptic longitude in radians of asteroid position
int shape_visibility_test = 1;															// Perform shape visibility tests (= 1) or not (= 0)
int roughness_visibility_test = 1;														// Perform surface roughness visibility tests (= 1) or not (= 0)
const int NTAU = 650;																	// Number of model time steps
const int m = 10;																		// Number of partial visibility divisions (splits facets into m*m smaller facets)
const int num_threads = 40;																// Number of parallel CPU threads to use

// Output Parameters

int write_shape_visibility = 1;															// Write out shape visibility (= 1) or not (= 0)
int write_roughness_visibility = 1;														// Write out surface roughness visibility (= 1) or not (= 0)

const std::string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
	return buf;
}

// Class AstVertex defining position of a single vertex of the global asteroid shape model


class AstVertex {

private:

	double _equpos[3];
	double _heliopos[NTAU][3];

public:

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

	float _shadow[NTAU];
	float _illumangle[NTAU];

public:

	SubAstFacet() { 
		for (int step=0;step!=NTAU;++step) {
			_shadow[step]=1.0f;
		}
	}

	float &rtnshadow(int t) { return _shadow[t]; }
	float &rtnillumangle(int t) { return _illumangle[t]; }

	~SubAstFacet() {}

};


// Class AstFacet defining surface properties of a single surface facet of the global asteroid shape model


class AstFacet {

private:

	int _astvertices[3];
	double _heliomidpos[NTAU][3];
	double _illumvector[NTAU][3];
	float _sfcillumvector[NTAU][3];
	double _normal[NTAU][3];
	double _illumangle[NTAU];
	int _shadow[NTAU];
	int _surfaceroughness;

	SubAstFacet *SF;

public:

	SubAstFacet*beginsf() const { return SF; }
	SubAstFacet*endsf() const { return SF+SFCNUMFACE; }

	AstFacet() {
		_surfaceroughness=0;
		for (int step=0;step!=NTAU;++step) {
			_shadow[step]=1;
		}
	}

	void createsubfacets() {
		SF=new SubAstFacet[SFCNUMFACE];
	}

	friend std::istream& operator>> (std::istream&i, AstFacet&A) {
		return i >> A._astvertices[0] >> A._astvertices[1] >> A._astvertices[2];
	}

	int &rtnastvertices(int n) { return _astvertices[n]; }
	double &rtnheliomidpos(int t, int n) { return _heliomidpos[t][n]; }
	double &rtnillumvector(int t, int n) { return _illumvector[t][n]; }
	float &rtnsfcillumvector(int t, int n) { return _sfcillumvector[t][n]; }
	double &rtnnormal(int t, int n) { return _normal[t][n]; }
	double &rtnillumangle(int t) { return _illumangle[t]; }
	int &rtnsurfaceroughness() { return _surfaceroughness; }
	int &rtnshadow(int t) { return _shadow[t]; }

	~AstFacet() {}

};


// Class Asteroid to control the asteroid global shape model


class Asteroid {

private:

	int NV, NF;
//	double _polcentre[3];
	double _carcentre[3];
	double _rotvector[3];
	
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
	//	if (polar_or_cart==1) {
	//		_polcentre[0]=polar0, _polcentre[1]=polar1, _polcentre[2]=polar2;
	//	} else {
	//		_carcentre[0]=cart0, _carcentre[1]=cart1, _carcentre[2]=cart2;
	//	}
		_rotvector[0]=rot0, _rotvector[1]=rot1;
	}

	void readmodel() {
		string filename;
		filename=shape_filename+".txt";
		ifstream read(filename.c_str());
		for (AstVertex *i=V;i!=V+NV;++i) {
			read >> *i;
		}
		for (AstFacet *j=F;j!=F+NF;++j) {
			read >> *j;
		}
	}

	void readshadow() {
		string shape_visibility_filename;
		shape_visibility_filename=shape_filename+"_visibility.txt";
		ifstream read(shape_visibility_filename.c_str());
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
		for (AstFacet *i=F;i!=F+NF;++i) {
			for (int step=0;step!=NTAU;++step) {
				read >> i->rtnshadow(step);
			}
		}
	}

	void addsurfaceroughness() {
		for (AstFacet *i=F;i!=F+NF;++i) {
			i->createsubfacets();
			i->rtnsurfaceroughness()=1;
		}
	}

	int &rtnnumvtx() { return NV; }
	int &rtnnumface() { return NF; }
//	double &rtnpolcentre(int n) { return _polcentre[n]; }
	double &rtncarcentre(int n) { return _carcentre[n]; }
	double &rtnrotvector(int n) { return _rotvector[n]; }
	
};


// Class SfcVertex defining position of a single vertex of the surface roughness shape model


class SfcVertex {

private:

	float _pos[3];

public:

	friend std::istream& operator>> (std::istream&i, SfcVertex&A) {
		return i >> A._pos[0] >> A._pos[1] >> A._pos[2];
	}

	float &rtnpos(int n) { return _pos[n]; }

	~SfcVertex() {}

};


// Class SfcFacet defining surface properties of a single surface facet of the surface roughness shape model


class SfcFacet {

private:

	int _sfcvertices[3];
	float _midpos[3];
	float _normcalc[3][3];
	float _normal[3];
	float _trivector[2][3];

public:

	friend std::istream& operator>> (std::istream&i, SfcFacet&A) {
		return i >> A._sfcvertices[0] >> A._sfcvertices[1] >> A._sfcvertices[2];
	}

	int &rtnsfcvertices(int n) { return _sfcvertices[n]; }
	float &rtnmidpos(int n) { return _midpos[n]; }
	float &rtnnormcalc(int n, int m) { return _normcalc[n][m]; }
	float &rtnnormal(int n) { return _normal[n]; }
	float &rtntrivector(int n, int m) { return _trivector[n][m]; }

	~SfcFacet() {}

};


// Class Surface to control the surface roughness shape model


class Surface {

private:

	int NV, NF;
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

	void readmodel() {
		string filename;
		filename=roughness_filename+".txt";
		ifstream read(filename.c_str());
		for (SfcVertex *i=SV;i!=SV+NV;++i) {
			read >> *i;
		}
		for (SfcFacet *j=SF;j!=SF+NF;++j) {
			read >> *j;
		}
	}

	int &rtnnumvtx() { return NV; }
	int &rtnnumface() { return NF; }

};


// Class Observer to specify postion of observer


class Observer {

private:

	double _polarpos[3];
	double _equpos[3];

public:

	Observer() {
		_polarpos[0]=ob0, _polarpos[1]=ob1, _polarpos[2]=ob2;
	}

	double &rtnpolarpos(int n) { return _polarpos[n]; }
	double &rtnequpos(int n) { return _equpos[n]; }

	~Observer() {}

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
	GPU_SfcVertex *SV;
	GPU_SfcFacet *SF;

	GPU_Surface(int n, int m): NV(n), NF(m) {
		SV=new GPU_SfcVertex[n];
		SF=new GPU_SfcFacet[m];
	}

};


// Structure Roughness_Shadow containing roughness facet shadow and illumination angle information


struct Roughness_Shadow {

	int shadow[NTAU];
	float sfc_shadow[NTAU][SFCNUMFACE];
	float sfc_illumvector[NTAU][3];
	float sfc_illumangle[NTAU][SFCNUMFACE];

};


// Determine surface roughness geometry


void surfacegeometry(Surface &S, Asteroid &A) {

	float n; 

	for (SfcFacet *C=S.beginf();C!=S.endf();++C) {

		C->rtnmidpos(0)=((S.beginv()+C->rtnsfcvertices(0))->rtnpos(0)+(S.beginv()+C->rtnsfcvertices(1))->rtnpos(0)+(S.beginv()+C->rtnsfcvertices(2))->rtnpos(0))/3.0f;
		C->rtnmidpos(1)=((S.beginv()+C->rtnsfcvertices(0))->rtnpos(1)+(S.beginv()+C->rtnsfcvertices(1))->rtnpos(1)+(S.beginv()+C->rtnsfcvertices(2))->rtnpos(1))/3.0f;
		C->rtnmidpos(2)=((S.beginv()+C->rtnsfcvertices(0))->rtnpos(2)+(S.beginv()+C->rtnsfcvertices(1))->rtnpos(2)+(S.beginv()+C->rtnsfcvertices(2))->rtnpos(2))/3.0f;
	
		C->rtnnormcalc(0,0)=(S.beginv()+C->rtnsfcvertices(1))->rtnpos(0)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(0);
		C->rtnnormcalc(0,1)=(S.beginv()+C->rtnsfcvertices(1))->rtnpos(1)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(1);
		C->rtnnormcalc(0,2)=(S.beginv()+C->rtnsfcvertices(1))->rtnpos(2)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(2);

		C->rtnnormcalc(1,0)=(S.beginv()+C->rtnsfcvertices(2))->rtnpos(0)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(0);
		C->rtnnormcalc(1,1)=(S.beginv()+C->rtnsfcvertices(2))->rtnpos(1)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(1);
		C->rtnnormcalc(1,2)=(S.beginv()+C->rtnsfcvertices(2))->rtnpos(2)-(S.beginv()+C->rtnsfcvertices(0))->rtnpos(2);

		C->rtnnormcalc(2,0)=(S.beginv()+C->rtnsfcvertices(2))->rtnpos(0)-(S.beginv()+C->rtnsfcvertices(1))->rtnpos(0);
		C->rtnnormcalc(2,1)=(S.beginv()+C->rtnsfcvertices(2))->rtnpos(1)-(S.beginv()+C->rtnsfcvertices(1))->rtnpos(1);
		C->rtnnormcalc(2,2)=(S.beginv()+C->rtnsfcvertices(2))->rtnpos(2)-(S.beginv()+C->rtnsfcvertices(1))->rtnpos(2);

		C->rtnnormal(0)=(C->rtnnormcalc(0,1)*C->rtnnormcalc(1,2))-(C->rtnnormcalc(0,2)*C->rtnnormcalc(1,1));
		C->rtnnormal(1)=(C->rtnnormcalc(0,2)*C->rtnnormcalc(1,0))-(C->rtnnormcalc(0,0)*C->rtnnormcalc(1,2));
		C->rtnnormal(2)=(C->rtnnormcalc(0,0)*C->rtnnormcalc(1,1))-(C->rtnnormcalc(0,1)*C->rtnnormcalc(1,0));

		n=sqrt((C->rtnnormal(0)*C->rtnnormal(0))+(C->rtnnormal(1)*C->rtnnormal(1))+(C->rtnnormal(2)*C->rtnnormal(2)));

		C->rtnnormal(0)=C->rtnnormal(0)/n;
		C->rtnnormal(1)=C->rtnnormal(1)/n;
		C->rtnnormal(2)=C->rtnnormal(2)/n;

		C->rtntrivector(0,0)=-1.0f*C->rtnnormcalc(1,0)/m;
		C->rtntrivector(0,1)=-1.0f*C->rtnnormcalc(1,1)/m;
		C->rtntrivector(0,2)=-1.0f*C->rtnnormcalc(1,2)/m;

		C->rtntrivector(1,0)=C->rtnnormcalc(0,0)/m;
		C->rtntrivector(1,1)=C->rtnnormcalc(0,1)/m;
		C->rtntrivector(1,2)=C->rtnnormcalc(0,2)/m;

		C->rtnnormcalc(1,0)=-1.0f*C->rtnnormcalc(1,0);
		C->rtnnormcalc(1,1)=-1.0f*C->rtnnormcalc(1,1);
		C->rtnnormcalc(1,2)=-1.0f*C->rtnnormcalc(1,2);

	}

	int step, i;
	float _t;
	AstFacet *C;

omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic) private(i,_t,C)
	for (step=0;step<NTAU;step++) {
		for (C=A.beginf();C!=A.endf();++C) {
			_t=sqrt((C->rtnsfcillumvector(step,0)*C->rtnsfcillumvector(step,0))+(C->rtnsfcillumvector(step,1)*C->rtnsfcillumvector(step,1))+(C->rtnsfcillumvector(step,2)*C->rtnsfcillumvector(step,2)));
			for (i=0;i!=SFCNUMFACE;++i) {
				(C->beginsf()+i)->rtnillumangle(step)=((C->rtnsfcillumvector(step,0)*(S.beginf()+i)->rtnnormal(0))+(C->rtnsfcillumvector(step,1)*(S.beginf()+i)->rtnnormal(1))+(C->rtnsfcillumvector(step,2)*(S.beginf()+i)->rtnnormal(2)))/_t;
				if ( C->rtnshadow(step)==0 ) { (C->beginsf()+i)->rtnshadow(step)=0.0f; }
				if ( (C->beginsf()+i)->rtnillumangle(step)<=0.0f ) { (C->beginsf()+i)->rtnshadow(step)=0.0f; }
			}
		}
	}

};


// Determine asteroid global shape model geometry


void asteroidgeometry(Asteroid &A, Observer &O) {

//	if (polar_or_cart==1) {
//		A.rtncarcentre(0)=AU*A.rtnpolcentre(0)*cos(A.rtnpolcentre(1))*cos(A.rtnpolcentre(2));
//		A.rtncarcentre(1)=AU*A.rtnpolcentre(0)*cos(A.rtnpolcentre(1))*sin(A.rtnpolcentre(2));
//		A.rtncarcentre(2)=AU*A.rtnpolcentre(0)*sin(A.rtnpolcentre(1));
//	}

	A.rtncarcentre(0)=AU*O.rtnpolarpos(0)*cos(O.rtnpolarpos(1))*cos(O.rtnpolarpos(2));
	A.rtncarcentre(1)=AU*O.rtnpolarpos(0)*cos(O.rtnpolarpos(1))*sin(O.rtnpolarpos(2));
	A.rtncarcentre(2)=AU*O.rtnpolarpos(0)*sin(O.rtnpolarpos(1));

//	O.rtnequpos(0)=AU*O.rtnpolarpos(0)*cos(O.rtnpolarpos(1))*cos(O.rtnpolarpos(2));
//	O.rtnequpos(1)=AU*O.rtnpolarpos(0)*cos(O.rtnpolarpos(1))*sin(O.rtnpolarpos(2));
//	O.rtnequpos(2)=AU*O.rtnpolarpos(0)*sin(O.rtnpolarpos(1));
	
//	A.rtncarcentre(0)=A.rtncarcentre(0)-O.rtnequpos(0);
//	A.rtncarcentre(1)=A.rtncarcentre(1)-O.rtnequpos(1);
//	A.rtncarcentre(2)=A.rtncarcentre(2)-O.rtnequpos(2);

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

	double rotangle, rotequpos[3], eclipos[3];
	double normcalc[2][3], y[3], _n, _t, _x;
	double imatrix[3][3]; 
	int step;
	AstVertex *B;
	AstFacet *C;

omp_set_num_threads(num_threads);
#pragma	omp parallel for schedule(dynamic) private(rotangle, rotequpos, eclipos, normcalc, y, _n, _t, _x, imatrix, B, C)
	for (step=0;step<NTAU;step++) {

		rotangle=2.0*PI*step/NTAU;

		for (B=A.beginv();B!=A.endv();++B) {

			rotequpos[0]=(B->rtnequpos(0)*cos(rotangle))-(B->rtnequpos(1)*sin(rotangle));
			rotequpos[1]=(B->rtnequpos(0)*sin(rotangle))+(B->rtnequpos(1)*cos(rotangle));
			rotequpos[2]=B->rtnequpos(2);

			eclipos[0]=(tmatrix[0][0]*rotequpos[0])+(tmatrix[0][1]*rotequpos[1])+(tmatrix[0][2]*rotequpos[2]);
			eclipos[1]=(tmatrix[1][0]*rotequpos[0])+(tmatrix[1][1]*rotequpos[1])+(tmatrix[1][2]*rotequpos[2]);
			eclipos[2]=(tmatrix[2][0]*rotequpos[0])+(tmatrix[2][1]*rotequpos[1])+(tmatrix[2][2]*rotequpos[2]);

			B->rtnheliopos(step,0)=A.rtncarcentre(0)+eclipos[0];
			B->rtnheliopos(step,1)=A.rtncarcentre(1)+eclipos[1];
			B->rtnheliopos(step,2)=A.rtncarcentre(2)+eclipos[2];

		}

		for (C=A.beginf();C!=A.endf();++C) {

			C->rtnheliomidpos(step,0)=(((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,0))+((A.beginv()+C->rtnastvertices(1))->rtnheliopos(step,0))+((A.beginv()+C->rtnastvertices(2))->rtnheliopos(step,0)))/3.0;
			C->rtnheliomidpos(step,1)=(((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,1))+((A.beginv()+C->rtnastvertices(1))->rtnheliopos(step,1))+((A.beginv()+C->rtnastvertices(2))->rtnheliopos(step,1)))/3.0;
			C->rtnheliomidpos(step,2)=(((A.beginv()+C->rtnastvertices(0))->rtnheliopos(step,2))+((A.beginv()+C->rtnastvertices(1))->rtnheliopos(step,2))+((A.beginv()+C->rtnastvertices(2))->rtnheliopos(step,2)))/3.0;

			_t=sqrt((C->rtnheliomidpos(step,0)*C->rtnheliomidpos(step,0))+(C->rtnheliomidpos(step,1)*C->rtnheliomidpos(step,1))+(C->rtnheliomidpos(step,2)*C->rtnheliomidpos(step,2)));

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

			C->rtnillumvector(step,0)=-1.0*((imatrix[0][0]*C->rtnheliomidpos(step,0))+(imatrix[0][1]*C->rtnheliomidpos(step,1))+(imatrix[0][2]*C->rtnheliomidpos(step,2)));
			C->rtnillumvector(step,1)=-1.0*((imatrix[1][0]*C->rtnheliomidpos(step,0))+(imatrix[1][1]*C->rtnheliomidpos(step,1))+(imatrix[1][2]*C->rtnheliomidpos(step,2)));
			C->rtnillumvector(step,2)=-1.0*((imatrix[2][0]*C->rtnheliomidpos(step,0))+(imatrix[2][1]*C->rtnheliomidpos(step,1))+(imatrix[2][2]*C->rtnheliomidpos(step,2)));

			C->rtnillumangle(step)=acos(C->rtnillumvector(step,2)/_t);
			if (C->rtnillumangle(step)>(PI/2.0)) {C->rtnshadow(step)=0;}

			C->rtnsfcillumvector(step,0)=(float)(10000.0*C->rtnillumvector(step,0)/_t);
			C->rtnsfcillumvector(step,1)=(float)(10000.0*C->rtnillumvector(step,1)/_t);
			C->rtnsfcillumvector(step,2)=(float)(10000.0*C->rtnillumvector(step,2)/_t);

		}

	}

};


// Code to test for asteroid shadows


int astshadowtest(Asteroid &A, AstFacet *G, AstFacet *H, int step) {

	double t, tcalc1, tcalc2;
	double xpoint[3];
	double normcalc1[3], normcalc2[3], normcalc3[3];
	double xvector1[3], xvector2[3], xvector3[3];
	double test1[3], test2[3], test3[3];
	double condition1, condition2, condition3;

	tcalc2=(G->rtnheliomidpos(step,0)*H->rtnnormal(step,0))+(G->rtnheliomidpos(step,1)*H->rtnnormal(step,1))+(G->rtnheliomidpos(step,2)*H->rtnnormal(step,2));	
	tcalc1=(H->rtnheliomidpos(step,0)*H->rtnnormal(step,0))+(H->rtnheliomidpos(step,1)*H->rtnnormal(step,1))+(H->rtnheliomidpos(step,2)*H->rtnnormal(step,2));
	t=tcalc1/tcalc2;

	xpoint[0]=t*G->rtnheliomidpos(step,0);
	xpoint[1]=t*G->rtnheliomidpos(step,1);
	xpoint[2]=t*G->rtnheliomidpos(step,2);

	normcalc1[0]=((A.beginv()+H->rtnastvertices(1))->rtnheliopos(step,0))-((A.beginv()+H->rtnastvertices(0))->rtnheliopos(step,0));
	normcalc1[1]=((A.beginv()+H->rtnastvertices(1))->rtnheliopos(step,1))-((A.beginv()+H->rtnastvertices(0))->rtnheliopos(step,1));
	normcalc1[2]=((A.beginv()+H->rtnastvertices(1))->rtnheliopos(step,2))-((A.beginv()+H->rtnastvertices(0))->rtnheliopos(step,2));

	normcalc2[0]=((A.beginv()+H->rtnastvertices(0))->rtnheliopos(step,0))-((A.beginv()+H->rtnastvertices(2))->rtnheliopos(step,0));
	normcalc2[1]=((A.beginv()+H->rtnastvertices(0))->rtnheliopos(step,1))-((A.beginv()+H->rtnastvertices(2))->rtnheliopos(step,1));
	normcalc2[2]=((A.beginv()+H->rtnastvertices(0))->rtnheliopos(step,2))-((A.beginv()+H->rtnastvertices(2))->rtnheliopos(step,2));

	normcalc3[0]=((A.beginv()+H->rtnastvertices(2))->rtnheliopos(step,0))-((A.beginv()+H->rtnastvertices(1))->rtnheliopos(step,0));
	normcalc3[1]=((A.beginv()+H->rtnastvertices(2))->rtnheliopos(step,1))-((A.beginv()+H->rtnastvertices(1))->rtnheliopos(step,1));
	normcalc3[2]=((A.beginv()+H->rtnastvertices(2))->rtnheliopos(step,2))-((A.beginv()+H->rtnastvertices(1))->rtnheliopos(step,2));

	xvector1[0]=xpoint[0]-((A.beginv()+H->rtnastvertices(0))->rtnheliopos(step,0));
	xvector1[1]=xpoint[1]-((A.beginv()+H->rtnastvertices(0))->rtnheliopos(step,1));
	xvector1[2]=xpoint[2]-((A.beginv()+H->rtnastvertices(0))->rtnheliopos(step,2));

	xvector3[0]=xpoint[0]-((A.beginv()+H->rtnastvertices(1))->rtnheliopos(step,0));
	xvector3[1]=xpoint[1]-((A.beginv()+H->rtnastvertices(1))->rtnheliopos(step,1));
	xvector3[2]=xpoint[2]-((A.beginv()+H->rtnastvertices(1))->rtnheliopos(step,2));

	xvector2[0]=xpoint[0]-((A.beginv()+H->rtnastvertices(2))->rtnheliopos(step,0));
	xvector2[1]=xpoint[1]-((A.beginv()+H->rtnastvertices(2))->rtnheliopos(step,1));
	xvector2[2]=xpoint[2]-((A.beginv()+H->rtnastvertices(2))->rtnheliopos(step,2));

	test1[0]=(normcalc1[1]*xvector1[2])-(normcalc1[2]*xvector1[1]);
    test1[1]=(normcalc1[2]*xvector1[0])-(normcalc1[0]*xvector1[2]);
	test1[2]=(normcalc1[0]*xvector1[1])-(normcalc1[1]*xvector1[0]);

	test2[0]=(normcalc2[1]*xvector2[2])-(normcalc2[2]*xvector2[1]);
	test2[1]=(normcalc2[2]*xvector2[0])-(normcalc2[0]*xvector2[2]);
	test2[2]=(normcalc2[0]*xvector2[1])-(normcalc2[1]*xvector2[0]);

	test3[0]=(normcalc3[1]*xvector3[2])-(normcalc3[2]*xvector3[1]);
	test3[1]=(normcalc3[2]*xvector3[0])-(normcalc3[0]*xvector3[2]);
	test3[2]=(normcalc3[0]*xvector3[1])-(normcalc3[1]*xvector3[0]);

	condition1=(test1[0]*H->rtnnormal(step,0))+(test1[1]*H->rtnnormal(step,1))+(test1[2]*H->rtnnormal(step,2));
	condition2=(test2[0]*H->rtnnormal(step,0))+(test2[1]*H->rtnnormal(step,1))+(test2[2]*H->rtnnormal(step,2));
	condition3=(test3[0]*H->rtnnormal(step,0))+(test3[1]*H->rtnnormal(step,1))+(test3[2]*H->rtnnormal(step,2));

	if (condition1 >= 0.0 && condition2 >= 0.0 && condition3 >= 0.0) { return 0; } else { return 1;}
	
};


// Code to determine asteroid shadowing


void astshadow(Asteroid &A, int step) {

	double g, h;

	for (AstFacet *G=A.beginf();G!=A.endf();++G) {

		g=sqrt((G->rtnheliomidpos(step,0)*G->rtnheliomidpos(step,0))+(G->rtnheliomidpos(step,1)*G->rtnheliomidpos(step,1))+(G->rtnheliomidpos(step,2)*G->rtnheliomidpos(step,2)));

		for (AstFacet *H=(G+1);H!=A.endf();++H) {

			h=sqrt((H->rtnheliomidpos(step,0)*H->rtnheliomidpos(step,0))+(H->rtnheliomidpos(step,1)*H->rtnheliomidpos(step,1))+(H->rtnheliomidpos(step,2)*H->rtnheliomidpos(step,2)));

			if (g>h) {

				if (G->rtnshadow(step)==1) {
					if ( astshadowtest(A,G,H,step)==0 ) { G->rtnshadow(step)=0; }
				}

			}	else {

				if (H->rtnshadow(step)==1) {
					if ( astshadowtest(A,H,G,step)==0 ) { H->rtnshadow(step)=0; }
				}

			}

		}

	}

};


// Load GPU constants


extern "C" void initialise_GPU_constants(GPU_Surface S_H);		// GPU code forward declaration
void load_GPU_constants(Surface &S, Asteroid &A) {

	int i;
	GPU_Surface S_H(SFCNUMVTX,SFCNUMFACE);
	SfcVertex *SV;
	SfcFacet *SF;

	for (i=0;i!=SFCNUMVTX;++i) {
		SV=S.beginv()+i;
		S_H.SV[i].pos[0]=SV->rtnpos(0);
		S_H.SV[i].pos[1]=SV->rtnpos(1);
		S_H.SV[i].pos[2]=SV->rtnpos(2);
	}

	for (i=0;i!=SFCNUMFACE;++i) {
		SF=S.beginf()+i;
		S_H.SF[i].sfcvertices[0]=SF->rtnsfcvertices(0);
		S_H.SF[i].sfcvertices[1]=SF->rtnsfcvertices(1);
		S_H.SF[i].sfcvertices[2]=SF->rtnsfcvertices(2);
		S_H.SF[i].midpos[0]=SF->rtnmidpos(0);
		S_H.SF[i].midpos[1]=SF->rtnmidpos(1);
		S_H.SF[i].midpos[2]=SF->rtnmidpos(2);
		S_H.SF[i].trivector[0][0]=SF->rtntrivector(0,0);
		S_H.SF[i].trivector[0][1]=SF->rtntrivector(0,1);
		S_H.SF[i].trivector[0][2]=SF->rtntrivector(0,2);
		S_H.SF[i].trivector[1][0]=SF->rtntrivector(1,0);
		S_H.SF[i].trivector[1][1]=SF->rtntrivector(1,1);
		S_H.SF[i].trivector[1][2]=SF->rtntrivector(1,2);
		S_H.SF[i].normcalc[0][0]=SF->rtnnormcalc(0,0);
		S_H.SF[i].normcalc[0][1]=SF->rtnnormcalc(0,1);
		S_H.SF[i].normcalc[0][2]=SF->rtnnormcalc(0,2);
		S_H.SF[i].normcalc[1][0]=SF->rtnnormcalc(1,0);
		S_H.SF[i].normcalc[1][1]=SF->rtnnormcalc(1,1);
		S_H.SF[i].normcalc[1][2]=SF->rtnnormcalc(1,2);
		S_H.SF[i].normcalc[2][0]=SF->rtnnormcalc(2,0);
		S_H.SF[i].normcalc[2][1]=SF->rtnnormcalc(2,1);
		S_H.SF[i].normcalc[2][2]=SF->rtnnormcalc(2,2);
		S_H.SF[i].normal[0]=SF->rtnnormal(0);
		S_H.SF[i].normal[1]=SF->rtnnormal(1);
		S_H.SF[i].normal[2]=SF->rtnnormal(2);
	}

	initialise_GPU_constants(S_H);

};


// Perform GPU surface roughness shadow data


extern "C" void GPU_sfcpartialshadow(Roughness_Shadow *RS_H);		// GPU code forward declaration
void perform_GPU_roughness_shadows(Asteroid &A) {

	Roughness_Shadow *RS_H;
	size_t size = ASTNUMFACE*sizeof(Roughness_Shadow);
	RS_H=(Roughness_Shadow *)malloc(size);

	AstFacet *C;
	SubAstFacet *SC;

	for (int i=0;i!=ASTNUMFACE;++i) {
		C=A.beginf()+i;
		for (int step=0;step!=NTAU;++step) {
			RS_H[i].shadow[step]=C->rtnshadow(step);
			RS_H[i].sfc_illumvector[step][0]=C->rtnsfcillumvector(step,0);
			RS_H[i].sfc_illumvector[step][1]=C->rtnsfcillumvector(step,1);
			RS_H[i].sfc_illumvector[step][2]=C->rtnsfcillumvector(step,2);
			for (int j=0;j!=SFCNUMFACE;++j) {
				SC=C->beginsf()+j;
				RS_H[i].sfc_shadow[step][j]=SC->rtnshadow(step);
				RS_H[i].sfc_illumangle[step][j]=SC->rtnillumangle(step);
			}
		}	
	}

	GPU_sfcpartialshadow(RS_H);
	free(RS_H);

};


// Retrieve GPU surface roughness shadow data


extern "C" void obtain_GPU_roughness_shadows(Roughness_Shadow *RS_H);		// GPU code forward declaration
void retrieve_GPU_roughness_shadows(Asteroid &A) {

	Roughness_Shadow *RS_H;
	size_t size = ASTNUMFACE*sizeof(Roughness_Shadow);
	RS_H=(Roughness_Shadow *)malloc(size);

	obtain_GPU_roughness_shadows(RS_H);

	AstFacet *C;
	SubAstFacet *SC;

	for (int i=0;i!=ASTNUMFACE;++i) {
		C=A.beginf()+i;
		for (int step=0;step!=NTAU;++step) {
			C->rtnshadow(step)=RS_H[i].shadow[step];
			C->rtnsfcillumvector(step,0)=RS_H[i].sfc_illumvector[step][0];
			C->rtnsfcillumvector(step,1)=RS_H[i].sfc_illumvector[step][1];
			C->rtnsfcillumvector(step,2)=RS_H[i].sfc_illumvector[step][2];
			for (int j=0;j!=SFCNUMFACE;++j) {
				SC=C->beginsf()+j;
				SC->rtnshadow(step)=RS_H[i].sfc_shadow[step][j];
				SC->rtnillumangle(step)=RS_H[i].sfc_illumangle[step][j];
			}
		}	
	}
	
	free(RS_H);

};


// Write out global shadow map


void writeglobalshadow(Asteroid &A) {

	string shape_visibility_writename;
	shape_visibility_writename=shape_filename+"_visibility.txt";
	ofstream write(shape_visibility_writename.c_str());

	for (AstVertex *B=A.beginv();B!=A.endv();++B) {
		write << B->rtnequpos(0) << "\t" << B->rtnequpos(1) << "\t" << B->rtnequpos(2) << endl;
	}

	for (AstFacet *C=A.beginf();C!=A.endf();++C) {
		write << C->rtnastvertices(0) << "\t" << C->rtnastvertices(1) << "\t" << C->rtnastvertices(2) << endl;
	}

	for (AstFacet *C=A.beginf();C!=A.endf();++C) {
		for (int step=0;step!=NTAU;++step) {
			write << C->rtnshadow(step) << endl;
		}
	}
	
};


// Binary write out surface roughness temperature maps


void binary_writesfcroughshadow(Surface &S, Asteroid &A) {

	string roughness_visibility_writename;
	roughness_visibility_writename=shape_filename+"_roughness_visibility.dat";
	ofstream write(roughness_visibility_writename.c_str(),ios::out|ios::binary);

	float pos[3];
	int vert[3];
	float val[NTAU];

	for (SfcVertex *B=S.beginv();B!=S.endv();++B) {
		pos[0]=B->rtnpos(0), pos[1]=B->rtnpos(1), pos[2]=B->rtnpos(2);
		write.write((char *)&pos,sizeof(float)*3);
	}

	for (SfcFacet *C=S.beginf();C!=S.endf();++C) {
		vert[0]=C->rtnsfcvertices(0), vert[1]=C->rtnsfcvertices(1), vert[2]=C->rtnsfcvertices(2);
		write.write((char *)&vert,sizeof(int)*3);
	}

	for (AstFacet *C=A.beginf();C!=A.endf();++C) {
		for (SubAstFacet *D=C->beginsf();D!=C->endsf();++D) {
			for (int step=0;step!=NTAU;++step) {
				val[step]=D->rtnshadow(step);
			}
			write.write((char *)&val,sizeof(float)*NTAU);
		}
	}

};


// Other GPU code forward declarations


extern "C" void set_cuda_device();
extern "C" void initialise_GPU_shadow_arrays();
extern "C" void freeup_GPU_shadow_arrays();	


// Main program loop


int main() {

	string status_filename;
	status_filename=shape_filename+"_visibility_status.txt";
	ofstream status(status_filename.c_str());

	status << "Start time" << "\t" << currentDateTime() << endl;

	Asteroid A(ASTNUMVTX, ASTNUMFACE);
	Surface S(SFCNUMVTX, SFCNUMFACE);
	Observer O;

	if (roughness_visibility_test==1) { 
		set_cuda_device();
		A.addsurfaceroughness(); 
	}

	A.readmodel();
	if (roughness_visibility_test==1) { S.readmodel(); }

	status << "Loaded asteroid models" << "\t" << currentDateTime() << endl;

	asteroidgeometry(A,O);

	if (load_shape_visibility==1) { A.readshadow(); }

	status << "Initialised asteroid geometry" << "\t" << currentDateTime() << endl;

	int step;

	if (shape_visibility_test==1) {
omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic)
		for (step=0;step<NTAU;step++) {
			astshadow(A,step);
			status << "Asteroid visibility step " << step << " calculated" << "\t" << currentDateTime() << endl;
		}
	}

	if (write_shape_visibility==1) { writeglobalshadow(A); }

	if (shape_visibility_test==1) {
		status << "Asteroid visibility calculated" << "\t" << currentDateTime() << endl;
	}

	if (roughness_visibility_test==1) {
		surfacegeometry(S,A);
		load_GPU_constants(S,A);
		status << "Initialised surface geometry" << "\t" << currentDateTime() << endl;
	}

	if (roughness_visibility_test==1) {
		initialise_GPU_shadow_arrays();
		perform_GPU_roughness_shadows(A); 
		retrieve_GPU_roughness_shadows(A);
		freeup_GPU_shadow_arrays();
	}
	
	if (write_roughness_visibility==1 && roughness_visibility_test==1) { binary_writesfcroughshadow(S,A); }

	if (roughness_visibility_test==1) {
		status << "Surface visibility calculated" << "\t" << currentDateTime() << endl;
	}

	status << "End time" << "\t" << currentDateTime() << endl;

};
