#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

const float PI = 4.0f*atan(1.0f);										// Value of Pi

// Input Parameters

const int ASTNUMVTX = 3204;											// Number of model vertices
const int ASTNUMFACE = 6400;											// Number of model facets
string filename = "doppler_binary";							// Model input file name
const int n = 10;														// Number of subfacet divisions (splits facets into n*n smaller facets)
const int num_threads = 10;												// Number of parallel threads to use


// Class Vertex defining position of a single vertex


class Vertex {

private:

	float _equpos[3];

public:

	friend std::istream& operator>> (std::istream&i, Vertex &A) {
		return i >> A._equpos[0] >> A._equpos[1] >> A._equpos[2];
	}

	float &rtnequpos(int n) { return _equpos[n]; }

	~Vertex() {}

};


// Class DepAstFacet to contain viewfactor data of dependent facets of global shape model


class DepAstFacet {

private:

	int _depastfacetnum;
	float _viewfactor;
	float _viewfactorvector[3];

public:

	int &rtndepastfacetnum() { return _depastfacetnum; }
	float &rtnviewfactor() { return _viewfactor; }
	float &rtnviewfactorvector(int n) { return _viewfactorvector[n]; }

	~DepAstFacet() {}

};


// Class Facet defining physical properties of a single surface facet


class Facet {

private:

	int _vertices[3];
	float _midpos[3];
	float _normal[3];
	float _area;
//	double _sourcevector[num_threads][3];
//	float _temporary;
	float _totalviewfactor;
	float _photonmomvector[3];
	float _normcalc1[3];
	float _normcalc2[3];
	float _normcalc3[3];

	int _dependantfacets;
	DepAstFacet *DAF;

	int _dependantfacets2;
	DepAstFacet *DAF2;

public:

	DepAstFacet*begindaf() const { return DAF; }
	DepAstFacet*enddaf() const { return DAF+_dependantfacets; }

	DepAstFacet*begindaf2() const { return DAF2; }
	DepAstFacet*enddaf2() const { return DAF2+_dependantfacets2; }

	Facet() {
		_vertices[0]=0, _vertices[1]=1, _vertices[2]=2;
//		_temporary=0.0f;
		_dependantfacets=0;
		_dependantfacets2=0;
		_totalviewfactor=0.0f;
		_photonmomvector[0]=0.0f;
		_photonmomvector[1]=0.0f;
		_photonmomvector[2]=0.0f;
	}

	void createdepfacets() {
		DAF=new DepAstFacet[_dependantfacets];
	}

	void createdepfacets2() {
		DAF2=new DepAstFacet[_dependantfacets];
	}

	void deletedepfacets() {
		delete [] DAF;
	}

	friend std::istream& operator>> (std::istream&i, Facet&A) {
		return i >> A._vertices[0] >> A._vertices[1] >> A._vertices[2];
	}

	int &rtnvertices(int n) { return _vertices[n]; }
	float &rtnmidpos(int n) { return _midpos[n]; }
	float &rtnnormal(int n) { return _normal[n]; }
	float &rtnarea() { return _area; }
//	double &rtnsourcevector(int n, int m) { return _sourcevector[n][m]; }
//	float &rtntemporary() { return _temporary; }
	float &rtntotalviewfactor() { return _totalviewfactor; }
	float &rtnphotonmomvector(int n) { return _photonmomvector[n]; }
	float &rtnnormcalc1(int n) { return _normcalc1[n]; }
	float &rtnnormcalc2(int n) { return _normcalc2[n]; }
	float &rtnnormcalc3(int n) { return _normcalc3[n]; }
	int &rtndependantfacets() { return _dependantfacets; }
	int &rtndependantfacets2() { return _dependantfacets2; }

	~Facet() {}

};


// Class Asteroid to create and control N Vertices and M Facets


class Asteroid {

private:

	int NV, NF;
	Vertex *E;
	Facet *F;

public:

	Vertex*beginv() const { return E; }
	Vertex*endv() const { return E+NV; }

	Facet*beginf() const { return F; }
	Facet*endf() const { return F+NF; }

	Asteroid(int n, int m) : NV(n),NF(m) {
		E=new Vertex[n];
		F=new Facet[m];
	}	

	void readmodel() {
		string shape_filename;
		shape_filename=filename+".txt";
		ifstream read(shape_filename.c_str());
		for (Vertex *i=E;i!=E+NV;++i) {
			read >> *i;
		}
		for (Facet *j=F;j!=F+NF;++j) {
			read >> *j;
		}
	}

	int &rtnnumvtx() { return NV; }
	int &rtnnumface() { return NF; }

};


// Code to determine whether a facet is hidden


void hidden(Asteroid &A) {

	float t, tcalc1, tcalc2;
	float xpoint[3];
//	double normcalc1[3], normcalc2[3], normcalc3[3];
	float xvector1[3], xvector2[3], xvector3[3];
	float test1[3], test2[3], test3[3];
	float condition1, condition2, condition3;
	float gvector[3], svector[3];
	float g, s;

	int viewexist;
	int facetnumber;

	Facet *G;
	Facet *Z;
	Facet *S;

	Vertex *V0, *V1, *V2;

	int i, j; 
//	int tid;

omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic, 100) private(t,tcalc1,tcalc2,xpoint,xvector1,xvector2,xvector3,test1,test2,test3,condition1,condition2,condition3,g,s,viewexist,gvector,svector,G,Z,S,V0,V1,V2,i,j)
	for (facetnumber=0;facetnumber<ASTNUMFACE;facetnumber++) {

//		tid=omp_get_thread_num();
		G=A.beginf()+facetnumber;
		G->createdepfacets2();

		if (G->rtndependantfacets()!=0) {

			for (i=0;i!=G->rtndependantfacets();++i) {			

				viewexist=1;
				Z=A.beginf()+(G->begindaf()+i)->rtndepastfacetnum(); //

			//	for (Z=A.beginf();Z!=A.endf();++Z) {
			//		Z->rtnsourcevector(tid,0)=Z->rtnmidpos(0)-(A.beginf()+(G->begindaf()+i)->rtndepastfacetnum())->rtnmidpos(0);
			//		Z->rtnsourcevector(tid,1)=Z->rtnmidpos(1)-(A.beginf()+(G->begindaf()+i)->rtndepastfacetnum())->rtnmidpos(1);
			//		Z->rtnsourcevector(tid,2)=Z->rtnmidpos(2)-(A.beginf()+(G->begindaf()+i)->rtndepastfacetnum())->rtnmidpos(2);
			//	}

				gvector[0]=G->rtnmidpos(0)-Z->rtnmidpos(0); //
				gvector[1]=G->rtnmidpos(1)-Z->rtnmidpos(1); //
				gvector[2]=G->rtnmidpos(2)-Z->rtnmidpos(2); //

			//	for (S=A.beginf();S!=A.endf();++S) {
				for (j=0;j!=G->rtndependantfacets();++j) {     //

					S=A.beginf()+(G->begindaf()+j)->rtndepastfacetnum();   //
					svector[0]=S->rtnmidpos(0)-Z->rtnmidpos(0); //
					svector[1]=S->rtnmidpos(1)-Z->rtnmidpos(1); //
					svector[2]=S->rtnmidpos(2)-Z->rtnmidpos(2); //

					if (S!=G && S!=Z) {

					//	g=sqrt((G->rtnsourcevector(tid,0)*G->rtnsourcevector(tid,0))+(G->rtnsourcevector(tid,1)*G->rtnsourcevector(tid,1))+(G->rtnsourcevector(tid,2)*G->rtnsourcevector(tid,2)));
					//	s=sqrt((S->rtnsourcevector(tid,0)*S->rtnsourcevector(tid,0))+(S->rtnsourcevector(tid,1)*S->rtnsourcevector(tid,1))+(S->rtnsourcevector(tid,2)*S->rtnsourcevector(tid,2)));

						g=sqrt((gvector[0]*gvector[0])+(gvector[1]*gvector[1])+(gvector[2]*gvector[2]));
						s=sqrt((svector[0]*svector[0])+(svector[1]*svector[1])+(svector[2]*svector[2]));

						if (g > s) { 

						//	tcalc1=(S->rtnsourcevector(tid,0)*S->rtnnormal(0))+(S->rtnsourcevector(tid,1)*S->rtnnormal(1))+(S->rtnsourcevector(tid,2)*S->rtnnormal(2));	
						//	tcalc2=(G->rtnsourcevector(tid,0)*S->rtnnormal(0))+(G->rtnsourcevector(tid,1)*S->rtnnormal(1))+(G->rtnsourcevector(tid,2)*S->rtnnormal(2));

							tcalc1=(svector[0]*S->rtnnormal(0))+(svector[1]*S->rtnnormal(1))+(svector[2]*S->rtnnormal(2));	
							tcalc2=(gvector[0]*S->rtnnormal(0))+(gvector[1]*S->rtnnormal(1))+(gvector[2]*S->rtnnormal(2));
							t=tcalc1/tcalc2;

							if (t > 0.0f) {

							//	xpoint[0]=(t*G->rtnsourcevector(tid,0))+(A.beginf()+(G->begindaf()+i)->rtndepastfacetnum())->rtnmidpos(0);
							//	xpoint[1]=(t*G->rtnsourcevector(tid,1))+(A.beginf()+(G->begindaf()+i)->rtndepastfacetnum())->rtnmidpos(1);
							//	xpoint[2]=(t*G->rtnsourcevector(tid,2))+(A.beginf()+(G->begindaf()+i)->rtndepastfacetnum())->rtnmidpos(2);

								xpoint[0]=(t*gvector[0])+Z->rtnmidpos(0); //
								xpoint[1]=(t*gvector[1])+Z->rtnmidpos(1); //
								xpoint[2]=(t*gvector[2])+Z->rtnmidpos(2); //

							//	normcalc1[0]=((A.beginv()+S->rtnvertices(1))->rtnequpos(0))-((A.beginv()+S->rtnvertices(0))->rtnequpos(0));
							//	normcalc1[1]=((A.beginv()+S->rtnvertices(1))->rtnequpos(1))-((A.beginv()+S->rtnvertices(0))->rtnequpos(1));
							//	normcalc1[2]=((A.beginv()+S->rtnvertices(1))->rtnequpos(2))-((A.beginv()+S->rtnvertices(0))->rtnequpos(2));

							//	normcalc2[0]=((A.beginv()+S->rtnvertices(0))->rtnequpos(0))-((A.beginv()+S->rtnvertices(2))->rtnequpos(0));
							//	normcalc2[1]=((A.beginv()+S->rtnvertices(0))->rtnequpos(1))-((A.beginv()+S->rtnvertices(2))->rtnequpos(1));
							//	normcalc2[2]=((A.beginv()+S->rtnvertices(0))->rtnequpos(2))-((A.beginv()+S->rtnvertices(2))->rtnequpos(2));
		
							//	normcalc3[0]=((A.beginv()+S->rtnvertices(2))->rtnequpos(0))-((A.beginv()+S->rtnvertices(1))->rtnequpos(0));
							//	normcalc3[1]=((A.beginv()+S->rtnvertices(2))->rtnequpos(1))-((A.beginv()+S->rtnvertices(1))->rtnequpos(1));
							//	normcalc3[2]=((A.beginv()+S->rtnvertices(2))->rtnequpos(2))-((A.beginv()+S->rtnvertices(1))->rtnequpos(2));

							//  normcalc1[0]=S->rtnnormcalc1(0);
							//	normcalc1[1]=S->rtnnormcalc1(1);
							//	normcalc1[2]=S->rtnnormcalc1(2);

							//	normcalc2[0]=S->rtnnormcalc2(0);
							//	normcalc2[1]=S->rtnnormcalc2(1);
							//	normcalc2[2]=S->rtnnormcalc2(2);

							//	normcalc3[0]=S->rtnnormcalc3(0);
							//	normcalc3[1]=S->rtnnormcalc3(1);
							//	normcalc3[2]=S->rtnnormcalc3(2);

								V0=A.beginv()+S->rtnvertices(0);
								V1=A.beginv()+S->rtnvertices(1);
								V2=A.beginv()+S->rtnvertices(2);

								xvector1[0]=xpoint[0]-V0->rtnequpos(0);
								xvector1[1]=xpoint[1]-V0->rtnequpos(1);
								xvector1[2]=xpoint[2]-V0->rtnequpos(2);

								xvector3[0]=xpoint[0]-V1->rtnequpos(0);
								xvector3[1]=xpoint[1]-V1->rtnequpos(1);
								xvector3[2]=xpoint[2]-V1->rtnequpos(2);

								xvector2[0]=xpoint[0]-V2->rtnequpos(0);
								xvector2[1]=xpoint[1]-V2->rtnequpos(1);
								xvector2[2]=xpoint[2]-V2->rtnequpos(2);

								test1[0]=(S->rtnnormcalc1(1)*xvector1[2])-(S->rtnnormcalc1(2)*xvector1[1]);
								test1[1]=(S->rtnnormcalc1(2)*xvector1[0])-(S->rtnnormcalc1(0)*xvector1[2]);
								test1[2]=(S->rtnnormcalc1(0)*xvector1[1])-(S->rtnnormcalc1(1)*xvector1[0]);

								test2[0]=(S->rtnnormcalc2(1)*xvector2[2])-(S->rtnnormcalc2(2)*xvector2[1]);
								test2[1]=(S->rtnnormcalc2(2)*xvector2[0])-(S->rtnnormcalc2(0)*xvector2[2]);
								test2[2]=(S->rtnnormcalc2(0)*xvector2[1])-(S->rtnnormcalc2(1)*xvector2[0]);

								test3[0]=(S->rtnnormcalc3(1)*xvector3[2])-(S->rtnnormcalc3(2)*xvector3[1]);
								test3[1]=(S->rtnnormcalc3(2)*xvector3[0])-(S->rtnnormcalc3(0)*xvector3[2]);
								test3[2]=(S->rtnnormcalc3(0)*xvector3[1])-(S->rtnnormcalc3(1)*xvector3[0]);

							//	test1[0]=(normcalc1[1]*xvector1[2])-(normcalc1[2]*xvector1[1]);
							//	test1[1]=(normcalc1[2]*xvector1[0])-(normcalc1[0]*xvector1[2]);
							//	test1[2]=(normcalc1[0]*xvector1[1])-(normcalc1[1]*xvector1[0]);

							//	test2[0]=(normcalc2[1]*xvector2[2])-(normcalc2[2]*xvector2[1]);
							//	test2[1]=(normcalc2[2]*xvector2[0])-(normcalc2[0]*xvector2[2]);
							//	test2[2]=(normcalc2[0]*xvector2[1])-(normcalc2[1]*xvector2[0]);

							//	test3[0]=(normcalc3[1]*xvector3[2])-(normcalc3[2]*xvector3[1]);
							//	test3[1]=(normcalc3[2]*xvector3[0])-(normcalc3[0]*xvector3[2]);
							//	test3[2]=(normcalc3[0]*xvector3[1])-(normcalc3[1]*xvector3[0]);

								condition1=(test1[0]*S->rtnnormal(0))+(test1[1]*S->rtnnormal(1))+(test1[2]*S->rtnnormal(2));
								condition2=(test2[0]*S->rtnnormal(0))+(test2[1]*S->rtnnormal(1))+(test2[2]*S->rtnnormal(2));
								condition3=(test3[0]*S->rtnnormal(0))+(test3[1]*S->rtnnormal(1))+(test3[2]*S->rtnnormal(2));

								if (condition1 >= 0.0f && condition2 >= 0.0f && condition3 >= 0.0f) { viewexist=0; } 
							
							}

						} 
	
					}
				
				}

				if (viewexist==1) {

					(G->begindaf2()+G->rtndependantfacets2())->rtndepastfacetnum()=(G->begindaf()+i)->rtndepastfacetnum();
					(G->begindaf2()+G->rtndependantfacets2())->rtnviewfactor()=(G->begindaf()+i)->rtnviewfactor();
					(G->begindaf2()+G->rtndependantfacets2())->rtnviewfactorvector(0)=(G->begindaf()+i)->rtnviewfactorvector(0);
					(G->begindaf2()+G->rtndependantfacets2())->rtnviewfactorvector(1)=(G->begindaf()+i)->rtnviewfactorvector(1);
					(G->begindaf2()+G->rtndependantfacets2())->rtnviewfactorvector(2)=(G->begindaf()+i)->rtnviewfactorvector(2);
					G->rtndependantfacets2()=G->rtndependantfacets2()+1;
					
				}

			}

		}

		cout << facetnumber << "\t" << G->rtndependantfacets2() << endl;
		G->deletedepfacets();
		
	}

};


// Code to calculate better viewfactors


void betterviewfactors(Asteroid &A) {

	float sepvector[3];
	float sepdistance;
	float ni, nj;
	float costhetai, costhetaj;

	int const subvtx((n+1)*(n+2)/2);
	int const subface(n*n);
	float test;

	float tri1vector1[3], tri1vector2[3];
	float tri2vector1[3], tri2vector2[3];

	float tri1subarea;
	float tri2subarea;

	float tri1vpos[subvtx][3];
	int tri1facet[subface][3];
	float tri1facetmidpos[subface][3];

	float tri2vpos[subvtx][3];
	int tri2facet[subface][3];
	float tri2facetmidpos[subface][3];

	float viewfactor;

	float tri1totalviewfactor[subface];
	float tri1effviewfactor;

	float tri1totalviewfactorvector[subface][3];
	float tri1effviewfactorvector[3];

	int l, i, j, k, facetnumber;

	Facet *I;

omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic, 100) private(sepvector,sepdistance,ni,nj,costhetai,costhetaj,test,tri1vector1,tri1vector2,tri2vector1,tri2vector2,tri1subarea,tri2subarea,tri1vpos,tri1facet,tri1facetmidpos,tri2vpos,tri2facet,tri2facetmidpos,viewfactor,tri1totalviewfactor,tri1effviewfactor,tri1totalviewfactorvector,tri1effviewfactorvector,l,i,j,k,I)
	for (facetnumber=0;facetnumber<ASTNUMFACE;facetnumber++) {

		I=A.beginf()+facetnumber;

		if (I->rtndependantfacets2()!=0) {

			ni=sqrt((I->rtnnormal(0)*I->rtnnormal(0))+(I->rtnnormal(1)*I->rtnnormal(1))+(I->rtnnormal(2)*I->rtnnormal(2)));

			tri1vector1[0]=((A.beginv()+I->rtnvertices(1))->rtnequpos(0)-(A.beginv()+I->rtnvertices(0))->rtnequpos(0))/n;
			tri1vector1[1]=((A.beginv()+I->rtnvertices(1))->rtnequpos(1)-(A.beginv()+I->rtnvertices(0))->rtnequpos(1))/n;
			tri1vector1[2]=((A.beginv()+I->rtnvertices(1))->rtnequpos(2)-(A.beginv()+I->rtnvertices(0))->rtnequpos(2))/n;

			tri1vector2[0]=((A.beginv()+I->rtnvertices(2))->rtnequpos(0)-(A.beginv()+I->rtnvertices(0))->rtnequpos(0))/n;
			tri1vector2[1]=((A.beginv()+I->rtnvertices(2))->rtnequpos(1)-(A.beginv()+I->rtnvertices(0))->rtnequpos(1))/n;
			tri1vector2[2]=((A.beginv()+I->rtnvertices(2))->rtnequpos(2)-(A.beginv()+I->rtnvertices(0))->rtnequpos(2))/n;

			tri1subarea=I->rtnarea()/(n*n);

			for (i=0;i!=I->rtndependantfacets2();++i) {

				sepvector[0]=(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnmidpos(0)-I->rtnmidpos(0);
				sepvector[1]=(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnmidpos(1)-I->rtnmidpos(1);
				sepvector[2]=(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnmidpos(2)-I->rtnmidpos(2);

				sepdistance=sqrt((sepvector[0]*sepvector[0])+(sepvector[1]*sepvector[1])+(sepvector[2]*sepvector[2]));

			//	test=(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnarea()/sepdistance;
				test=(I->begindaf2()+i)->rtnviewfactor();

			//	if (test>0.1) {
				if (test>=0.05f) {

					for (j=0;j!=subface;j++) {
						tri1totalviewfactor[j]=0.0f;
						tri1totalviewfactorvector[j][0]=0.0f;
						tri1totalviewfactorvector[j][1]=0.0f;
						tri1totalviewfactorvector[j][2]=0.0f;
					}

					tri1effviewfactor=0.0f;
					tri1effviewfactorvector[0]=0.0f;
					tri1effviewfactorvector[1]=0.0f;
					tri1effviewfactorvector[2]=0.0f;

					nj=sqrt(((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnnormal(0)*(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnnormal(0))+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnnormal(1)*(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnnormal(1))+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnnormal(2)*(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnnormal(2)));

					tri2subarea=(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnarea()/(n*n);

					tri2vector1[0]=((A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(1)))->rtnequpos(0)-(A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(0)))->rtnequpos(0))/n;
					tri2vector1[1]=((A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(1)))->rtnequpos(1)-(A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(0)))->rtnequpos(1))/n;
					tri2vector1[2]=((A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(1)))->rtnequpos(2)-(A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(0)))->rtnequpos(2))/n;

					tri2vector2[0]=((A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(2)))->rtnequpos(0)-(A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(0)))->rtnequpos(0))/n;
					tri2vector2[1]=((A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(2)))->rtnequpos(1)-(A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(0)))->rtnequpos(1))/n;
					tri2vector2[2]=((A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(2)))->rtnequpos(2)-(A.beginv()+((A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(0)))->rtnequpos(2))/n;

					// Generate sub vertex tables

					l=0;

					for (j=1;j!=(n+2);j++) {
						for (k=1;k!=(j+1);k++) {

							// Triangle 1

							tri1vpos[l][0]=(A.beginv()+I->rtnvertices(2))->rtnequpos(0)+((k-1)*tri1vector1[0])-((j-1)*tri1vector2[0]);
							tri1vpos[l][1]=(A.beginv()+I->rtnvertices(2))->rtnequpos(1)+((k-1)*tri1vector1[1])-((j-1)*tri1vector2[1]);
							tri1vpos[l][2]=(A.beginv()+I->rtnvertices(2))->rtnequpos(2)+((k-1)*tri1vector1[2])-((j-1)*tri1vector2[2]);

							// Triangle 2

							tri2vpos[l][0]=(A.beginv()+(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(2))->rtnequpos(0)+((k-1)*tri2vector1[0])-((j-1)*tri2vector2[0]);
							tri2vpos[l][1]=(A.beginv()+(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(2))->rtnequpos(1)+((k-1)*tri2vector1[1])-((j-1)*tri2vector2[1]);
							tri2vpos[l][2]=(A.beginv()+(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnvertices(2))->rtnequpos(2)+((k-1)*tri2vector1[2])-((j-1)*tri2vector2[2]);

							l=l+1;

						}
					}

					//	Generate sub facet tables

					l=0;

					for (j=1;j!=(n+1);j++) {

						if (j!=1) {

							for (k=1;k!=j;k++) {

								// Triangle 1

								tri1facet[l][0]=((j*(j-1))/2)+k;
								tri1facet[l][1]=((j*(j-1))/2)+k-1;
								tri1facet[l][2]=((j*(j+1))/2)+k-1;

								// Triangle 2

								tri2facet[l][0]=((j*(j-1))/2)+k;
								tri2facet[l][1]=((j*(j-1))/2)+k-1;
								tri2facet[l][2]=((j*(j+1))/2)+k-1;

								l=l+1;

								// Triangle 1

								tri1facet[l][0]=((j*(j+1))/2)+k-1;
								tri1facet[l][1]=((j*(j+1))/2)+k;
								tri1facet[l][2]=((j*(j-1))/2)+k;

								// Triangle 2

								tri2facet[l][0]=((j*(j+1))/2)+k-1;
								tri2facet[l][1]=((j*(j+1))/2)+k;
								tri2facet[l][2]=((j*(j-1))/2)+k;

								l=l+1;

							}
	
						}

						// Triangle 1

						tri1facet[l][0]=(((j+1)*(j+2))/2)-2;
						tri1facet[l][1]=(((j+1)*(j+2))/2)-1;
						tri1facet[l][2]=((j*(j+1))/2)-1;

						// Triangle 2

						tri2facet[l][0]=(((j+1)*(j+2))/2)-2;
						tri2facet[l][1]=(((j+1)*(j+2))/2)-1;
						tri2facet[l][2]=((j*(j+1))/2)-1;

						l=l+1;

					}

					// Determine sub facet midpoints

					for (j=0;j!=subface;j++) {

						// Triangle 1

						tri1facetmidpos[j][0]=(tri1vpos[tri1facet[j][0]][0]+tri1vpos[tri1facet[j][1]][0]+tri1vpos[tri1facet[j][2]][0])/3.0f;
						tri1facetmidpos[j][1]=(tri1vpos[tri1facet[j][0]][1]+tri1vpos[tri1facet[j][1]][1]+tri1vpos[tri1facet[j][2]][1])/3.0f;
						tri1facetmidpos[j][2]=(tri1vpos[tri1facet[j][0]][2]+tri1vpos[tri1facet[j][1]][2]+tri1vpos[tri1facet[j][2]][2])/3.0f;

						// Triangle 2

						tri2facetmidpos[j][0]=(tri2vpos[tri2facet[j][0]][0]+tri2vpos[tri2facet[j][1]][0]+tri2vpos[tri2facet[j][2]][0])/3.0f;
						tri2facetmidpos[j][1]=(tri2vpos[tri2facet[j][0]][1]+tri2vpos[tri2facet[j][1]][1]+tri2vpos[tri2facet[j][2]][1])/3.0f;
						tri2facetmidpos[j][2]=(tri2vpos[tri2facet[j][0]][2]+tri2vpos[tri2facet[j][1]][2]+tri2vpos[tri2facet[j][2]][2])/3.0f;

					}

					// Determine total viewfactors for sub facets of triangle 1

					for (j=0;j!=subface;j++) {
						for (k=0;k!=subface;k++) {

							sepvector[0]=tri2facetmidpos[k][0]-tri1facetmidpos[j][0];
							sepvector[1]=tri2facetmidpos[k][1]-tri1facetmidpos[j][1];
							sepvector[2]=tri2facetmidpos[k][2]-tri1facetmidpos[j][2];

							sepdistance=sqrt((sepvector[0]*sepvector[0])+(sepvector[1]*sepvector[1])+(sepvector[2]*sepvector[2]));
	
							costhetai=((I->rtnnormal(0)*sepvector[0])+(I->rtnnormal(1)*sepvector[1])+(I->rtnnormal(2)*sepvector[2]))/(ni*sepdistance);
							costhetaj=((-1.0f*(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnnormal(0)*sepvector[0])+(-1.0f*(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnnormal(1)*sepvector[1])+(-1.0f*(A.beginf()+(I->begindaf2()+i)->rtndepastfacetnum())->rtnnormal(2)*sepvector[2]))/(nj*sepdistance);

							viewfactor=(costhetai*costhetaj*tri2subarea)/(PI*sepdistance*sepdistance);
							if ( viewfactor<0.0f ) { viewfactor=0.0f; }
							tri1totalviewfactor[j]=tri1totalviewfactor[j]+viewfactor;

							sepvector[0]=sepvector[0]/sepdistance;
							sepvector[1]=sepvector[1]/sepdistance;
							sepvector[2]=sepvector[2]/sepdistance;

							tri1totalviewfactorvector[j][0]=tri1totalviewfactorvector[j][0]+(viewfactor*sepvector[0]);
							tri1totalviewfactorvector[j][1]=tri1totalviewfactorvector[j][1]+(viewfactor*sepvector[1]);
							tri1totalviewfactorvector[j][2]=tri1totalviewfactorvector[j][2]+(viewfactor*sepvector[2]);

						}
					}

					// Determine effective viewfactor for triangle 1

					for (j=0;j!=subface;j++) {
						tri1effviewfactor=tri1effviewfactor+((tri1subarea/I->rtnarea())*tri1totalviewfactor[j]);

						tri1effviewfactorvector[0]=tri1effviewfactorvector[0]+(tri1subarea*tri1totalviewfactorvector[j][0]);
						tri1effviewfactorvector[1]=tri1effviewfactorvector[1]+(tri1subarea*tri1totalviewfactorvector[j][1]);
						tri1effviewfactorvector[2]=tri1effviewfactorvector[2]+(tri1subarea*tri1totalviewfactorvector[j][2]);
					}

					(I->begindaf2()+i)->rtnviewfactor()=tri1effviewfactor;

					(I->begindaf2()+i)->rtnviewfactorvector(0)=tri1effviewfactorvector[0];
					(I->begindaf2()+i)->rtnviewfactorvector(1)=tri1effviewfactorvector[1];
					(I->begindaf2()+i)->rtnviewfactorvector(2)=tri1effviewfactorvector[2];

				}

			}

		}

	}

};


// Code to determine asteroid selfheating


void selfheating(Asteroid &A) {

	float sepvector[3];
	float sepdistance;
	float ni, nj;
	float thetai, thetaj; 
//	float thetaI, thetaJ;
	int facetlocation;
	int facetnumber;
	Facet *I;
	Facet *J;

	// Determine maximum number of interfacing facets

omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic, 100) private(sepvector,sepdistance,ni,nj,thetai,thetaj,I,J)
	for (facetnumber=0;facetnumber<ASTNUMFACE;facetnumber++) {

		I=A.beginf()+facetnumber;

		for (J=A.beginf();J!=A.endf();++J) {

			if (I!=J) {

				sepvector[0]=J->rtnmidpos(0)-I->rtnmidpos(0);
				sepvector[1]=J->rtnmidpos(1)-I->rtnmidpos(1);
				sepvector[2]=J->rtnmidpos(2)-I->rtnmidpos(2);

				sepdistance=sqrt((sepvector[0]*sepvector[0])+(sepvector[1]*sepvector[1])+(sepvector[2]*sepvector[2]));

				ni=sqrt((I->rtnnormal(0)*I->rtnnormal(0))+(I->rtnnormal(1)*I->rtnnormal(1))+(I->rtnnormal(2)*I->rtnnormal(2)));
				nj=sqrt((J->rtnnormal(0)*J->rtnnormal(0))+(J->rtnnormal(1)*J->rtnnormal(1))+(J->rtnnormal(2)*J->rtnnormal(2)));

			//	thetaI=(I->rtnnormal(0)*sepvector[0])+(I->rtnnormal(1)*sepvector[1])+(I->rtnnormal(2)*sepvector[2]);
			//	thetaJ=(-1.0f*J->rtnnormal(0)*sepvector[0])+(-1.0f*J->rtnnormal(1)*sepvector[1])+(-1.0f*J->rtnnormal(2)*sepvector[2]);

				thetai=acos(((I->rtnnormal(0)*sepvector[0])+(I->rtnnormal(1)*sepvector[1])+(I->rtnnormal(2)*sepvector[2]))/(ni*sepdistance));
				thetaj=acos(((-1.0f*J->rtnnormal(0)*sepvector[0])+(-1.0f*J->rtnnormal(1)*sepvector[1])+(-1.0f*J->rtnnormal(2)*sepvector[2]))/(nj*sepdistance));
			
				if ( thetai <= (PI/2.0f) && thetaj <= (PI/2.0f)) { I->rtndependantfacets()=I->rtndependantfacets()+1;	}

			}

		}

		I->createdepfacets();
		cout << facetnumber << "\t" << I->rtndependantfacets() << endl;

	}

	// Determine interfacing facet properties

omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic, 100) private(sepvector,sepdistance,ni,nj,thetai,thetaj,facetlocation,I,J)
	for (facetnumber=0;facetnumber<ASTNUMFACE;facetnumber++) {

		I=A.beginf()+facetnumber;
		facetlocation=0;
		I->rtndependantfacets()=0;

		for (J=A.beginf();J!=A.endf();++J) {

			if (I!=J) {

				sepvector[0]=J->rtnmidpos(0)-I->rtnmidpos(0);
				sepvector[1]=J->rtnmidpos(1)-I->rtnmidpos(1);
				sepvector[2]=J->rtnmidpos(2)-I->rtnmidpos(2);

				sepdistance=sqrt((sepvector[0]*sepvector[0])+(sepvector[1]*sepvector[1])+(sepvector[2]*sepvector[2]));

				ni=sqrt((I->rtnnormal(0)*I->rtnnormal(0))+(I->rtnnormal(1)*I->rtnnormal(1))+(I->rtnnormal(2)*I->rtnnormal(2)));
				nj=sqrt((J->rtnnormal(0)*J->rtnnormal(0))+(J->rtnnormal(1)*J->rtnnormal(1))+(J->rtnnormal(2)*J->rtnnormal(2)));

			//	thetaI=(I->rtnnormal(0)*sepvector[0])+(I->rtnnormal(1)*sepvector[1])+(I->rtnnormal(2)*sepvector[2]);
			//	thetaJ=(-1.0f*J->rtnnormal(0)*sepvector[0])+(-1.0f*J->rtnnormal(1)*sepvector[1])+(-1.0f*J->rtnnormal(2)*sepvector[2]);

				thetai=acos(((I->rtnnormal(0)*sepvector[0])+(I->rtnnormal(1)*sepvector[1])+(I->rtnnormal(2)*sepvector[2]))/(ni*sepdistance));
				thetaj=acos(((-1.0f*J->rtnnormal(0)*sepvector[0])+(-1.0f*J->rtnnormal(1)*sepvector[1])+(-1.0f*J->rtnnormal(2)*sepvector[2]))/(nj*sepdistance));
			
				if ( thetai <= (PI/2.0f) && thetaj <= (PI/2.0f)) {

					(I->begindaf()+I->rtndependantfacets())->rtnviewfactor()=(cos(thetai)*cos(thetaj)*J->rtnarea())/(PI*sepdistance*sepdistance);
					if ((I->begindaf()+I->rtndependantfacets())->rtnviewfactor()<0.0f) { (I->begindaf()+I->rtndependantfacets())->rtnviewfactor()=0.0f; }

					sepvector[0]=sepvector[0]/sepdistance;
					sepvector[1]=sepvector[1]/sepdistance;
					sepvector[2]=sepvector[2]/sepdistance;

					(I->begindaf()+I->rtndependantfacets())->rtnviewfactorvector(0)=I->rtnarea()*(I->begindaf()+I->rtndependantfacets())->rtnviewfactor()*sepvector[0];
					(I->begindaf()+I->rtndependantfacets())->rtnviewfactorvector(1)=I->rtnarea()*(I->begindaf()+I->rtndependantfacets())->rtnviewfactor()*sepvector[1];
					(I->begindaf()+I->rtndependantfacets())->rtnviewfactorvector(2)=I->rtnarea()*(I->begindaf()+I->rtndependantfacets())->rtnviewfactor()*sepvector[2];
					
					(I->begindaf()+I->rtndependantfacets())->rtndepastfacetnum()=facetlocation;

					I->rtndependantfacets()=I->rtndependantfacets()+1;

				}

			}

			facetlocation=facetlocation+1;

		}

		cout << facetnumber << endl;

	}

	// Determine hidden interfacing facets

	hidden(A);

	// Calculate better viewfactors

	betterviewfactors(A);

};


// Determine geometry of asteroid


void initialise(Asteroid &A) {

	int facetnumber;
	Facet *C;

	// Determine midpoint vector of all facets

omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic, 100) private(C)
	for (facetnumber=0;facetnumber<ASTNUMFACE;facetnumber++) {

		C=A.beginf()+facetnumber;

		C->rtnmidpos(0)=(((A.beginv()+C->rtnvertices(0))->rtnequpos(0))+((A.beginv()+C->rtnvertices(1))->rtnequpos(0))+((A.beginv()+C->rtnvertices(2))->rtnequpos(0)))/3.0f;
		C->rtnmidpos(1)=(((A.beginv()+C->rtnvertices(0))->rtnequpos(1))+((A.beginv()+C->rtnvertices(1))->rtnequpos(1))+((A.beginv()+C->rtnvertices(2))->rtnequpos(1)))/3.0f;
		C->rtnmidpos(2)=(((A.beginv()+C->rtnvertices(0))->rtnequpos(2))+((A.beginv()+C->rtnvertices(1))->rtnequpos(2))+((A.beginv()+C->rtnvertices(2))->rtnequpos(2)))/3.0f;

	}

	// Determine normal vector of all facets

	float _normcalc[2][3];
	float _n;

omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic, 100) private(_normcalc,_n,C)
	for (facetnumber=0;facetnumber<ASTNUMFACE;facetnumber++) {

		C=A.beginf()+facetnumber;

		_normcalc[0][0]=((A.beginv()+C->rtnvertices(1))->rtnequpos(0))-((A.beginv()+C->rtnvertices(0))->rtnequpos(0));
		_normcalc[0][1]=((A.beginv()+C->rtnvertices(1))->rtnequpos(1))-((A.beginv()+C->rtnvertices(0))->rtnequpos(1));
		_normcalc[0][2]=((A.beginv()+C->rtnvertices(1))->rtnequpos(2))-((A.beginv()+C->rtnvertices(0))->rtnequpos(2));

		_normcalc[1][0]=((A.beginv()+C->rtnvertices(2))->rtnequpos(0))-((A.beginv()+C->rtnvertices(0))->rtnequpos(0));
		_normcalc[1][1]=((A.beginv()+C->rtnvertices(2))->rtnequpos(1))-((A.beginv()+C->rtnvertices(0))->rtnequpos(1));	
		_normcalc[1][2]=((A.beginv()+C->rtnvertices(2))->rtnequpos(2))-((A.beginv()+C->rtnvertices(0))->rtnequpos(2));

		C->rtnnormal(0)=(_normcalc[0][1]*_normcalc[1][2])-(_normcalc[0][2]*_normcalc[1][1]);
		C->rtnnormal(1)=(_normcalc[0][2]*_normcalc[1][0])-(_normcalc[0][0]*_normcalc[1][2]);
		C->rtnnormal(2)=(_normcalc[0][0]*_normcalc[1][1])-(_normcalc[0][1]*_normcalc[1][0]);

		_n=sqrt((C->rtnnormal(0)*C->rtnnormal(0))+(C->rtnnormal(1)*C->rtnnormal(1))+(C->rtnnormal(2)*C->rtnnormal(2)));
		C->rtnarea()=_n/2.0f;

		C->rtnnormcalc1(0)=_normcalc[0][0];
		C->rtnnormcalc1(1)=_normcalc[0][1];
		C->rtnnormcalc1(2)=_normcalc[0][2];

		C->rtnnormcalc2(0)=-1.0f*_normcalc[1][0];
		C->rtnnormcalc2(1)=-1.0f*_normcalc[1][1];
		C->rtnnormcalc2(2)=-1.0f*_normcalc[1][2];
		
		C->rtnnormcalc3(0)=((A.beginv()+C->rtnvertices(2))->rtnequpos(0))-((A.beginv()+C->rtnvertices(1))->rtnequpos(0));
		C->rtnnormcalc3(1)=((A.beginv()+C->rtnvertices(2))->rtnequpos(1))-((A.beginv()+C->rtnvertices(1))->rtnequpos(1));
		C->rtnnormcalc3(2)=((A.beginv()+C->rtnvertices(2))->rtnequpos(2))-((A.beginv()+C->rtnvertices(1))->rtnequpos(2));

	}

	selfheating(A);
	
};


// Main program loop


int main() {
	
	Asteroid A(ASTNUMVTX, ASTNUMFACE);
	A.readmodel();
	
	initialise(A);

	for (Facet *C=A.beginf();C!=A.endf();++C) {
		if (C->rtndependantfacets2()!=0) {
			for (int i=0;i!=C->rtndependantfacets2();++i) {
				C->rtntotalviewfactor()=C->rtntotalviewfactor()+(C->begindaf2()+i)->rtnviewfactor();
				C->rtnphotonmomvector(0)=C->rtnphotonmomvector(0)+((C->begindaf2()+i)->rtnviewfactorvector(0)/C->rtnarea());
				C->rtnphotonmomvector(1)=C->rtnphotonmomvector(1)+((C->begindaf2()+i)->rtnviewfactorvector(1)/C->rtnarea());
				C->rtnphotonmomvector(2)=C->rtnphotonmomvector(2)+((C->begindaf2()+i)->rtnviewfactorvector(2)/C->rtnarea());
			}
		} else {
			C->rtntotalviewfactor()=0;
		}
	}

	string selfheating_filename = filename+"_selfheating_list.txt";
	string selfheating_vectors_filename = filename+"_selfheating_vector_list.txt";
	string selfheating_map_filename = filename+"_selfheating_map.txt";
	string totalviewfactor_filename = filename+"_total_viewfactor_map.txt";
	string photovectors_filename = filename+"_photovectors.txt";

	ofstream swrite(selfheating_map_filename.c_str());
	for (Vertex *B=A.beginv();B!=A.endv();++B) {
		swrite << B->rtnequpos(0) << "\t" << B->rtnequpos(1) << "\t" << B->rtnequpos(2) << endl;
	}
	for (Facet *C=A.beginf();C!=A.endf();++C) {
		swrite << C->rtnvertices(0) << "\t" << C->rtnvertices(1) << "\t" << C->rtnvertices(2) << "\t" << C->rtndependantfacets2() << endl;
	}

	ofstream swrite1(totalviewfactor_filename.c_str());
	for (Vertex *B=A.beginv();B!=A.endv();++B) {
		swrite1 << B->rtnequpos(0) << "\t" << B->rtnequpos(1) << "\t" << B->rtnequpos(2) << endl;
	}
	for (Facet *C=A.beginf();C!=A.endf();++C) {
		swrite1 << C->rtnvertices(0) << "\t" << C->rtnvertices(1) << "\t" << C->rtnvertices(2) << "\t" << C->rtntotalviewfactor() << endl;
	}

	ofstream swrite2(selfheating_filename.c_str());
	for (Facet *C=A.beginf();C!=A.endf();++C) {
		swrite2 << C->rtndependantfacets2() << "\t";
		if (C->rtndependantfacets2()!=0) {
			for(int i=0;i!=C->rtndependantfacets2();++i) {
				swrite2 << (C->begindaf2()+i)->rtndepastfacetnum() << "\t";
			}
			for(int i=0;i!=C->rtndependantfacets2();++i) {
				swrite2 << (C->begindaf2()+i)->rtnviewfactor() << "\t";
			}
		}
		swrite2 << endl;
	}

	ofstream swrite3(selfheating_vectors_filename.c_str());
	for (Facet *C=A.beginf();C!=A.endf();++C) {
		if (C->rtndependantfacets2()!=0) {
			for(int i=0;i!=C->rtndependantfacets2();++i) {
				swrite3 << ((C->begindaf2()+i)->rtnviewfactorvector(0)/((C->begindaf2()+i)->rtnviewfactor()*C->rtnarea())) << "\t";
				swrite3 << ((C->begindaf2()+i)->rtnviewfactorvector(1)/((C->begindaf2()+i)->rtnviewfactor()*C->rtnarea())) << "\t";
				swrite3 << ((C->begindaf2()+i)->rtnviewfactorvector(2)/((C->begindaf2()+i)->rtnviewfactor()*C->rtnarea())) << "\t";
			}
		}
		swrite3 << endl;
	}

	ofstream swrite4(photovectors_filename.c_str());
	for (Facet *C=A.beginf();C!=A.endf();++C) {
		swrite4 << C->rtnphotonmomvector(0) << "\t" << C->rtnphotonmomvector(1) << "\t" << C->rtnphotonmomvector(2) << endl;
	}

};