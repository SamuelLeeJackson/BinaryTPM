#include "cuda_runtime.h"
#include "device_launch_parameters.h"

const int ASTNUMFACE = 6400;											// Number of shape model facets
const int SFCNUMVTX = 61;												// Number of surface roughness model vertices
const int SFCNUMFACE = 100;												// Number of surface roughness model facets
const int NTAU = 650;													// Number of model time steps
const int m = 10;														// Number of partial shadow divisions (splits facets into m*m smaller facets)
const int cuda_device = 0;												// CUDA device number that is used


// Surface shadowing values stored in GPU constant memory


__constant__ float sfc_vertex[SFCNUMVTX][3];
__constant__ int sfc_facet[SFCNUMFACE][3];
__constant__ float sfc_midpos[SFCNUMFACE][3];
__constant__ float sfc_trivector[SFCNUMFACE][2][3];
__constant__ float sfc_normcalc[SFCNUMFACE][3][3];
__constant__ float sfc_normal[SFCNUMFACE][3];
__constant__ float trivector_multiplier[m*m][2];


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

};


// Structure Roughness_Shadow containing roughness facet shadow and illumination angle information


struct Roughness_Shadow {

	int shadow[NTAU];
	float sfc_shadow[NTAU][SFCNUMFACE];
	float sfc_illumvector[NTAU][3];
	float sfc_illumangle[NTAU][SFCNUMFACE];

};


// Pointers to GPU device memory


Roughness_Shadow *RS_D;


// Code to determine surface partial shadowing


__global__ void sfcpartialshadow(Roughness_Shadow *RS_D, int step) {

	__shared__ int shadow;
	__shared__ float sunvector[3];
	__shared__ float sfcshadow[SFCNUMFACE];
	
	int p;
	float midpos[3];
	float trans_co[2];
	float a, s;
	float shadow_radius; 
	float crater_radius = 10.0f;

	if (threadIdx.x==0) {
		shadow=RS_D[blockIdx.x].shadow[step];
	}
	__syncthreads();

	if (shadow==1) {

		if (threadIdx.x<4) {
			sunvector[threadIdx.x]=RS_D[blockIdx.x].sfc_illumvector[step][threadIdx.x];
		}
		sfcshadow[threadIdx.x]=RS_D[blockIdx.x].sfc_shadow[step][threadIdx.x];
		__syncthreads();

		if (threadIdx.x==0) {
			s=sqrt((sunvector[0]*sunvector[0])+(sunvector[1]*sunvector[1]));
			sunvector[0]=sunvector[0]/s;
			sunvector[1]=sunvector[1]/s;
			sunvector[2]=sunvector[2]/s;
		}
		__syncthreads();

		a=fabs(sunvector[2])/sqrt((sunvector[0]*sunvector[0])+(sunvector[1]*sunvector[1]));

		for (p=0;p!=(m*m);++p) {
		
			midpos[0]=sfc_vertex[sfc_facet[threadIdx.x][2]][0]+(trivector_multiplier[p][0]*sfc_trivector[threadIdx.x][0][0])+(trivector_multiplier[p][1]*sfc_trivector[threadIdx.x][1][0]);
			midpos[1]=sfc_vertex[sfc_facet[threadIdx.x][2]][1]+(trivector_multiplier[p][0]*sfc_trivector[threadIdx.x][0][1])+(trivector_multiplier[p][1]*sfc_trivector[threadIdx.x][1][1]);
			midpos[2]=sfc_vertex[sfc_facet[threadIdx.x][2]][2]+(trivector_multiplier[p][0]*sfc_trivector[threadIdx.x][0][2])+(trivector_multiplier[p][1]*sfc_trivector[threadIdx.x][1][2]);

			trans_co[1]=(-1.0f*sunvector[1]*midpos[0])+(sunvector[0]*midpos[1]);
			shadow_radius=sqrt((crater_radius*crater_radius)-(trans_co[1]*trans_co[1]));
			trans_co[0]=(sunvector[0]*midpos[0])+(sunvector[1]*midpos[1])-shadow_radius;

			if (midpos[2]<=(a*trans_co[0]) && sfcshadow[threadIdx.x]!=0.0f) { sfcshadow[threadIdx.x]=sfcshadow[threadIdx.x]-(1.0f/(float)(m*m)); }
			if (sfcshadow[threadIdx.x]<0.0f) { sfcshadow[threadIdx.x]=0.0f; }

		}
		__syncthreads();
	
		RS_D[blockIdx.x].sfc_shadow[step][threadIdx.x]=sfcshadow[threadIdx.x];

	}

	__syncthreads();

};


// Code to set which CUDA device to use


extern "C"
void set_cuda_device() {

	cudaSetDevice(cuda_device);

};


// Code to initialise GPU constant memory values


extern "C"
void initialise_GPU_constants(GPU_Surface S_H) {

	float host_sfc_vertex[SFCNUMVTX][3];
	int host_sfc_facet[SFCNUMFACE][3];
	float host_sfc_midpos[SFCNUMFACE][3];
	float host_sfc_trivector[SFCNUMFACE][2][3];
	float host_sfc_normcalc[SFCNUMFACE][3][3];
	float host_sfc_normal[SFCNUMFACE][3];

	for (int i=0;i!=SFCNUMVTX;++i) {
		host_sfc_vertex[i][0]=S_H.SV[i].pos[0];
		host_sfc_vertex[i][1]=S_H.SV[i].pos[1];
		host_sfc_vertex[i][2]=S_H.SV[i].pos[2];
	}

	for (int j=0;j!=SFCNUMFACE;++j) {
		host_sfc_facet[j][0]=S_H.SF[j].sfcvertices[0];
		host_sfc_facet[j][1]=S_H.SF[j].sfcvertices[1];
		host_sfc_facet[j][2]=S_H.SF[j].sfcvertices[2];
		host_sfc_midpos[j][0]=S_H.SF[j].midpos[0];
		host_sfc_midpos[j][1]=S_H.SF[j].midpos[1];
		host_sfc_midpos[j][2]=S_H.SF[j].midpos[2];
		host_sfc_trivector[j][0][0]=S_H.SF[j].trivector[0][0];
		host_sfc_trivector[j][0][1]=S_H.SF[j].trivector[0][1];
		host_sfc_trivector[j][0][2]=S_H.SF[j].trivector[0][2];
		host_sfc_trivector[j][1][0]=S_H.SF[j].trivector[1][0];
		host_sfc_trivector[j][1][1]=S_H.SF[j].trivector[1][1];
		host_sfc_trivector[j][1][2]=S_H.SF[j].trivector[1][2];
		host_sfc_normcalc[j][0][0]=S_H.SF[j].normcalc[0][0];
		host_sfc_normcalc[j][0][1]=S_H.SF[j].normcalc[0][1];
		host_sfc_normcalc[j][0][2]=S_H.SF[j].normcalc[0][2];
		host_sfc_normcalc[j][1][0]=S_H.SF[j].normcalc[1][0];
		host_sfc_normcalc[j][1][1]=S_H.SF[j].normcalc[1][1];
		host_sfc_normcalc[j][1][2]=S_H.SF[j].normcalc[1][2];
		host_sfc_normcalc[j][2][0]=S_H.SF[j].normcalc[2][0];
		host_sfc_normcalc[j][2][1]=S_H.SF[j].normcalc[2][1];
		host_sfc_normcalc[j][2][2]=S_H.SF[j].normcalc[2][2];
		host_sfc_normal[j][0]=S_H.SF[j].normal[0];
		host_sfc_normal[j][1]=S_H.SF[j].normal[1];
		host_sfc_normal[j][2]=S_H.SF[j].normal[2];
	}

	cudaMemcpyToSymbol(sfc_vertex,host_sfc_vertex,SFCNUMVTX*3*sizeof(float));
	cudaMemcpyToSymbol(sfc_facet,host_sfc_facet,SFCNUMFACE*3*sizeof(int));
	cudaMemcpyToSymbol(sfc_midpos,host_sfc_midpos,SFCNUMFACE*3*sizeof(float));
	cudaMemcpyToSymbol(sfc_trivector,host_sfc_trivector,SFCNUMFACE*2*3*sizeof(float));
	cudaMemcpyToSymbol(sfc_normcalc,host_sfc_normcalc,SFCNUMFACE*3*3*sizeof(float));
	cudaMemcpyToSymbol(sfc_normal,host_sfc_normal,SFCNUMFACE*3*sizeof(float));

	float host_trivector_multiplier[m*m][2];
	int k, l, p;

	for (k=1;k!=m+1;++k) {
		for (l=1;l!=k+1;++l) {
			p=((k-1)*(k-1))+(2*(l-1));
			host_trivector_multiplier[p][0]=(float)k-(1.0f/3.0f);
			host_trivector_multiplier[p][1]=(float)l-(2.0f/3.0f);
		}
	}

	for (k=2;k!=m+1;++k) {
		for (l=1;l!=k;++l) {
			p=((k-1)*(k-1))+(2*(l-1))+1;
			host_trivector_multiplier[p][0]=(float)k-(2.0f/3.0f);
			host_trivector_multiplier[p][1]=(float)l-(1.0f/3.0f);
		}
	}

	cudaMemcpyToSymbol(trivector_multiplier,host_trivector_multiplier,m*m*2*sizeof(float));

};


// Code to initialise GPU roughness shadow arrays


extern "C"
void initialise_GPU_shadow_arrays() {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Shadow);
	cudaMalloc((void **) &RS_D, size);

};


// Code to free up GPU roughness shadow arrays


extern "C"
void freeup_GPU_shadow_arrays() {

	cudaFree(RS_D);

};


// Code to obtain GPU roughness shadow data


extern "C" 
void obtain_GPU_roughness_shadows(Roughness_Shadow *RS_H) {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Shadow);
	cudaMemcpy(RS_H,RS_D,size,cudaMemcpyDeviceToHost);

};


// Surface roughness shadowing GPU control code


extern "C"
void GPU_sfcpartialshadow(Roughness_Shadow *RS_H) {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Shadow);
	cudaMemcpy(RS_D,RS_H,size,cudaMemcpyHostToDevice);

	for (int step=0;step!=NTAU;++step) {
		sfcpartialshadow <<<ASTNUMFACE,SFCNUMFACE>>> (RS_D,step);
		cudaDeviceSynchronize();
	}

};
