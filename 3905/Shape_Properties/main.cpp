#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

const double PI = 4.0*atan(1.0);										// Value of Pi

// Input Parameters

const int ASTNUMVTX = 1602;												// Number of shape model vertices
const int ASTNUMFACE = 3200;											// Number of shape model facets
string filename = "doppler_secondary";											// Shape model input file name


// Class AstVertex defining position of a single vertex of the global asteroid shape model


class AstVertex {

private:

	double _equpos[3];
	
public:

	AstVertex() {}

	friend std::istream& operator>> (std::istream&i, AstVertex &A) {
		return i >> A._equpos[0] >> A._equpos[1] >> A._equpos[2];
	}

	double &rtnequpos(int n) { return _equpos[n]; }
	
	~AstVertex() {}

};


// Class AstFacet defining surface properties of a single surface facet of the global asteroid shape model


class AstFacet {

private:

	int _astvertices[3];
	double _midpos[3];
	double _com[3];
	double _normal[3];
	double _radius;
	double _area;
	double _volume;
	bool _outward;

public:

	AstFacet() {
		_outward=true;
	}

	friend std::istream& operator>> (std::istream&i, AstFacet&A) {
		return i >> A._astvertices[0] >> A._astvertices[1] >> A._astvertices[2];
	}

	int &rtnastvertices(int n) { return _astvertices[n]; }
	double &rtnmidpos(int n) { return _midpos[n]; }
	double &rtncom(int n) { return _com[n]; }
	double &rtnnormal(int n) { return _normal[n]; }
	double &rtnradius() { return _radius; }
	double &rtnarea() { return _area; }
	double &rtnvolume() { return _volume; }
	bool &rtnoutward() { return _outward; }

	~AstFacet() {}

};


// Class Asteroid to control the asteroid global shape model


class Asteroid {

private:

	int NV, NF;

	AstVertex *V;
	AstFacet *F;

	double _totalarea;
	double _totalvolume;
	double _minx, _maxx;
	double _miny, _maxy;
	double _minz, _maxz;
	double _effdiameter;
	double _vertexcentroid[3];
	double _areacentroid[3];
	double _volumecentroid[3];

public:

	AstVertex*beginv() const { return V; }
	AstVertex*endv() const { return V+NV; }

	AstFacet*beginf() const { return F; }
	AstFacet*endf() const { return F+NF; }
	
	Asteroid(int n, int m) : NV(n), NF(m) {
		V=new AstVertex[n];
		F=new AstFacet[m];
		_totalarea=0.0;
		_totalvolume=0.0;
		_minx=0.0, _maxx=0.0;
		_miny=0.0, _maxy=0.0;
		_minz=0.0, _maxz=0.0;
		_vertexcentroid[0]=0.0, _vertexcentroid[1]=0.0, _vertexcentroid[2]=0.0;
		_areacentroid[0]=0.0, _areacentroid[1]=0.0, _areacentroid[2]=0.0;
		_volumecentroid[0]=0.0, _volumecentroid[1]=0.0, _volumecentroid[2]=0.0;
	}

	void readmodel() {
		string shape_filename;
		shape_filename=filename+".txt";
		ifstream read(shape_filename.c_str());
		for (AstVertex *i=V;i!=V+NV;++i) {
			read >> *i;
		}
		for (AstFacet *j=F;j!=F+NF;++j) {
			read >> *j;
		}
	}

	int &rtnnumvtx() { return NV; }
	int &rtnnumface() { return NF; }

	double &rtntotalarea() { return _totalarea; }
	double &rtntotalvolume() { return _totalvolume; }
	double &rtnminx() { return _minx; }
	double &rtnmaxx() { return _maxx; }
	double &rtnminy() { return _miny; }
	double &rtnmaxy() { return _maxy; }
	double &rtnminz() { return _minz; }
	double &rtnmaxz() { return _maxz; }
	double &rtneffdiameter() { return _effdiameter; }
	double &rtnvertexcentroid(int n) { return _vertexcentroid[n]; }
	double &rtnareacentroid(int n) { return _areacentroid[n]; }
	double &rtnvolumecentroid(int n) { return _volumecentroid[n]; }

};


// Determine asteroid global shape model geometry


void asteroidgeometry(Asteroid &A) {

	for (AstVertex *B=A.beginv();B!=A.endv();++B) {

		A.rtnvertexcentroid(0)=A.rtnvertexcentroid(0)+B->rtnequpos(0);
		A.rtnvertexcentroid(1)=A.rtnvertexcentroid(1)+B->rtnequpos(1);
		A.rtnvertexcentroid(2)=A.rtnvertexcentroid(2)+B->rtnequpos(2);

		if ( B->rtnequpos(0) < A.rtnminx() ) { A.rtnminx()=B->rtnequpos(0); }
		if ( B->rtnequpos(0) > A.rtnmaxx() ) { A.rtnmaxx()=B->rtnequpos(0); }

		if ( B->rtnequpos(1) < A.rtnminy() ) { A.rtnminy()=B->rtnequpos(1); }
		if ( B->rtnequpos(1) > A.rtnmaxy() ) { A.rtnmaxy()=B->rtnequpos(1); }

		if ( B->rtnequpos(2) < A.rtnminz() ) { A.rtnminz()=B->rtnequpos(2); }
		if ( B->rtnequpos(2) > A.rtnmaxz() ) { A.rtnmaxz()=B->rtnequpos(2); }

	}

	A.rtnvertexcentroid(0)=A.rtnvertexcentroid(0)/ASTNUMVTX;
	A.rtnvertexcentroid(1)=A.rtnvertexcentroid(1)/ASTNUMVTX;
	A.rtnvertexcentroid(2)=A.rtnvertexcentroid(2)/ASTNUMVTX;

	double normcalc[2][3], _n;
	double costheta, theta;

	for (AstFacet *C=A.beginf();C!=A.endf();++C) {
	
		C->rtnmidpos(0)=((A.beginv()+C->rtnastvertices(0))->rtnequpos(0)+(A.beginv()+C->rtnastvertices(1))->rtnequpos(0)+(A.beginv()+C->rtnastvertices(2))->rtnequpos(0))/3.0;
		C->rtnmidpos(1)=((A.beginv()+C->rtnastvertices(0))->rtnequpos(1)+(A.beginv()+C->rtnastvertices(1))->rtnequpos(1)+(A.beginv()+C->rtnastvertices(2))->rtnequpos(1))/3.0;
		C->rtnmidpos(2)=((A.beginv()+C->rtnastvertices(0))->rtnequpos(2)+(A.beginv()+C->rtnastvertices(1))->rtnequpos(2)+(A.beginv()+C->rtnastvertices(2))->rtnequpos(2))/3.0;

		C->rtnradius()=sqrt((C->rtnmidpos(0)*C->rtnmidpos(0))+(C->rtnmidpos(1)*C->rtnmidpos(1))+(C->rtnmidpos(2)*C->rtnmidpos(2)));

		C->rtncom(0)=((A.beginv()+C->rtnastvertices(0))->rtnequpos(0)+(A.beginv()+C->rtnastvertices(1))->rtnequpos(0)+(A.beginv()+C->rtnastvertices(2))->rtnequpos(0))/4.0;
		C->rtncom(1)=((A.beginv()+C->rtnastvertices(0))->rtnequpos(1)+(A.beginv()+C->rtnastvertices(1))->rtnequpos(1)+(A.beginv()+C->rtnastvertices(2))->rtnequpos(1))/4.0;
		C->rtncom(2)=((A.beginv()+C->rtnastvertices(0))->rtnequpos(2)+(A.beginv()+C->rtnastvertices(1))->rtnequpos(2)+(A.beginv()+C->rtnastvertices(2))->rtnequpos(2))/4.0;

		normcalc[0][0]=((A.beginv()+C->rtnastvertices(1))->rtnequpos(0))-((A.beginv()+C->rtnastvertices(0))->rtnequpos(0));
		normcalc[0][1]=((A.beginv()+C->rtnastvertices(1))->rtnequpos(1))-((A.beginv()+C->rtnastvertices(0))->rtnequpos(1));
		normcalc[0][2]=((A.beginv()+C->rtnastvertices(1))->rtnequpos(2))-((A.beginv()+C->rtnastvertices(0))->rtnequpos(2));

		normcalc[1][0]=((A.beginv()+C->rtnastvertices(2))->rtnequpos(0))-((A.beginv()+C->rtnastvertices(0))->rtnequpos(0));
		normcalc[1][1]=((A.beginv()+C->rtnastvertices(2))->rtnequpos(1))-((A.beginv()+C->rtnastvertices(0))->rtnequpos(1));	
		normcalc[1][2]=((A.beginv()+C->rtnastvertices(2))->rtnequpos(2))-((A.beginv()+C->rtnastvertices(0))->rtnequpos(2));

		C->rtnnormal(0)=(normcalc[0][1]*normcalc[1][2])-(normcalc[0][2]*normcalc[1][1]);
		C->rtnnormal(1)=(normcalc[0][2]*normcalc[1][0])-(normcalc[0][0]*normcalc[1][2]);
		C->rtnnormal(2)=(normcalc[0][0]*normcalc[1][1])-(normcalc[0][1]*normcalc[1][0]);

		_n=sqrt((C->rtnnormal(0)*C->rtnnormal(0))+(C->rtnnormal(1)*C->rtnnormal(1))+(C->rtnnormal(2)*C->rtnnormal(2)));
		C->rtnarea()=_n/2.0;

		if ((C->rtnradius()*_n)!=0.0) { costheta=((C->rtnmidpos(0)*C->rtnnormal(0))+(C->rtnmidpos(1)*C->rtnnormal(1))+(C->rtnmidpos(2)*C->rtnnormal(2)))/(C->rtnradius()*_n); } else { costheta=0.0; }
		theta=acos(costheta);

		if ( costheta < 0.0 ) { C->rtnoutward()=false; }

		if ( C->rtnoutward()==true ) {
		
			C->rtnvolume()=C->rtnarea()*C->rtnradius()*costheta/3.0;
		
		} else {

			C->rtnvolume()=C->rtnarea()*C->rtnradius()*cos(PI-theta)/3.0;
		
		}

	}

	for (AstFacet *C=A.beginf();C!=A.endf();++C) {

		A.rtntotalarea()=A.rtntotalarea()+C->rtnarea();

		A.rtnareacentroid(0)=A.rtnareacentroid(0)+(C->rtnarea()*C->rtnmidpos(0));
		A.rtnareacentroid(1)=A.rtnareacentroid(1)+(C->rtnarea()*C->rtnmidpos(1));
		A.rtnareacentroid(2)=A.rtnareacentroid(2)+(C->rtnarea()*C->rtnmidpos(2));

		if ( C->rtnoutward()==true ) {

			A.rtntotalvolume()=A.rtntotalvolume()+C->rtnvolume();

			A.rtnvolumecentroid(0)=A.rtnvolumecentroid(0)+(C->rtnvolume()*C->rtncom(0));
			A.rtnvolumecentroid(1)=A.rtnvolumecentroid(1)+(C->rtnvolume()*C->rtncom(1));
			A.rtnvolumecentroid(2)=A.rtnvolumecentroid(2)+(C->rtnvolume()*C->rtncom(2));

		} else {

			A.rtntotalvolume()=A.rtntotalvolume()-C->rtnvolume();

			A.rtnvolumecentroid(0)=A.rtnvolumecentroid(0)-(C->rtnvolume()*C->rtncom(0));
			A.rtnvolumecentroid(1)=A.rtnvolumecentroid(1)-(C->rtnvolume()*C->rtncom(1));
			A.rtnvolumecentroid(2)=A.rtnvolumecentroid(2)-(C->rtnvolume()*C->rtncom(2));

		}

		if (C->rtnarea() == 0.0f) {
			cout << "Help i have a zero" << endl;
		} 

	}

	A.rtnareacentroid(0)=A.rtnareacentroid(0)/A.rtntotalarea();
	A.rtnareacentroid(1)=A.rtnareacentroid(1)/A.rtntotalarea();
	A.rtnareacentroid(2)=A.rtnareacentroid(2)/A.rtntotalarea();

	A.rtnvolumecentroid(0)=A.rtnvolumecentroid(0)/A.rtntotalvolume();
	A.rtnvolumecentroid(1)=A.rtnvolumecentroid(1)/A.rtntotalvolume();
	A.rtnvolumecentroid(2)=A.rtnvolumecentroid(2)/A.rtntotalvolume();

	A.rtneffdiameter()=2.0*pow((3.0/4.0)*(1.0/PI)*A.rtntotalvolume(),(1.0/3.0));

};


// Main program loop


int main() {

	Asteroid A(ASTNUMVTX,ASTNUMFACE);
	A.readmodel();

	asteroidgeometry(A);

	string properties_filename=filename+"_properties.txt";
	ofstream write(properties_filename.c_str());

	write << "Total Area:" << "\t" << A.rtntotalarea() << endl;
	write << "Total Volume:" << "\t" << A.rtntotalvolume() << endl;
	write << "Min X:" << "\t" << A.rtnminx() << "\t" << "Max X:" << "\t" << A.rtnmaxx() << endl;
	write << "Min Y:" << "\t" << A.rtnminy() << "\t" << "Max Y:" << "\t" << A.rtnmaxy() << endl;
	write << "Min Z:" << "\t" << A.rtnminz() << "\t" << "Max Z:" << "\t" << A.rtnmaxz() << endl;
	write << "Effective Diameter:" << "\t" << A.rtneffdiameter() << endl;
	write << "Vertex Centroid:" << "\t" << A.rtnvertexcentroid(0) << "\t" << A.rtnvertexcentroid(1) << "\t" << A.rtnvertexcentroid(2) << endl;
	write << "Area Centroid:" << "\t" << A.rtnareacentroid(0) << "\t" << A.rtnareacentroid(1) << "\t" << A.rtnareacentroid(2) << endl;
	write << "Volume Centroid:" << "\t" << A.rtnvolumecentroid(0) << "\t" << A.rtnvolumecentroid(1) << "\t" << A.rtnvolumecentroid(2) << endl;

};
