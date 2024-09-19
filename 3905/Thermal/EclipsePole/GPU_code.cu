#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// In-Built Input Parameters

const int ASTNUMFACE = 6400;											// Number of shape model facets
const int SFCNUMVTX = 61;												// Number of surface roughness model vertices
const int SFCNUMFACE = 100;												// Number of surface roughness model facets
const int NTAU = 650;													// Number of model time steps (multiples of 50)
const int NX = 56;														// Number of model depth steps
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


// Surface thermal modelling values stored in GPU constant memory


__constant__ float PI;
__constant__ float rot_period;


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
	float sfc_temp[NX+1][SFCNUMFACE];

};


// Pointers to GPU device memory


Facet_Properties *FP_D;
Roughness_Shadow *RS_D;
Roughness_Illumination *RI_D;
Roughness_Globalillumination *RGI_D;
Roughness_Globalheating *RGH_D;
Roughness_Temperatures *RT_D;
float *viewfactors_D;


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
			s=sqrtf((sunvector[0]*sunvector[0])+(sunvector[1]*sunvector[1]));
			sunvector[0]=sunvector[0]/s;
			sunvector[1]=sunvector[1]/s;
			sunvector[2]=sunvector[2]/s;
		}
		__syncthreads();

		a=fabsf(sunvector[2])/sqrtf((sunvector[0]*sunvector[0])+(sunvector[1]*sunvector[1]));

		for (p=0;p!=(m*m);++p) {
		
			midpos[0]=sfc_vertex[sfc_facet[threadIdx.x][2]][0]+(trivector_multiplier[p][0]*sfc_trivector[threadIdx.x][0][0])+(trivector_multiplier[p][1]*sfc_trivector[threadIdx.x][1][0]);
			midpos[1]=sfc_vertex[sfc_facet[threadIdx.x][2]][1]+(trivector_multiplier[p][0]*sfc_trivector[threadIdx.x][0][1])+(trivector_multiplier[p][1]*sfc_trivector[threadIdx.x][1][1]);
			midpos[2]=sfc_vertex[sfc_facet[threadIdx.x][2]][2]+(trivector_multiplier[p][0]*sfc_trivector[threadIdx.x][0][2])+(trivector_multiplier[p][1]*sfc_trivector[threadIdx.x][1][2]);

			trans_co[1]=(-1.0f*sunvector[1]*midpos[0])+(sunvector[0]*midpos[1]);
			shadow_radius=sqrtf((crater_radius*crater_radius)-(trans_co[1]*trans_co[1]));
			trans_co[0]=(sunvector[0]*midpos[0])+(sunvector[1]*midpos[1])-shadow_radius;

			if (midpos[2]<=(a*trans_co[0]) && sfcshadow[threadIdx.x]!=0.0f) { sfcshadow[threadIdx.x]=sfcshadow[threadIdx.x]-(1.0f/(float)(m*m)); }
			if (sfcshadow[threadIdx.x]<0.0f) { sfcshadow[threadIdx.x]=0.0f; }

		}
		__syncthreads();
	
		RS_D[blockIdx.x].sfc_shadow[step][threadIdx.x]=sfcshadow[threadIdx.x];

	}

	__syncthreads();

};


// Determine roughness facet illumination from direct and single scattered sunlight (not including thermal radiation!)


__global__ void roughness_single_scattering(Facet_Properties *FP, Roughness_Shadow *RS, Roughness_Illumination *RI, Roughness_Globalillumination *RGI, float *viewfactors, float solar, int step) {

	__shared__ float albedo;
	__shared__ float illumination[SFCNUMFACE];

	int i;
	float shadow, illumangle, globalillumination, totalillumination;

	if (threadIdx.x==0) {
		albedo=FP[blockIdx.x].albedo;
	}
	__syncthreads();

	shadow=RS[blockIdx.x].sfc_shadow[step][threadIdx.x];
	illumangle=RS[blockIdx.x].sfc_illumangle[step][threadIdx.x];
	globalillumination=RGI[blockIdx.x].sfc_globalillumination[step][threadIdx.x];
	__syncthreads();
	
	illumination[threadIdx.x]=(solar*shadow*illumangle)+globalillumination;
	totalillumination=illumination[threadIdx.x];

	for (i=0;i!=SFCNUMFACE;++i) {
		totalillumination=totalillumination+(viewfactors[(i*SFCNUMFACE)+threadIdx.x]*albedo*illumination[i]);
	}

	if (totalillumination<=0.0f) { totalillumination=0.0f; }
	__syncthreads();

	RI[blockIdx.x].sfc_illumination[step][threadIdx.x]=totalillumination;

};


// Determine roughness facet illumination from direct and multiple scattered sunlight (not including thermal radiation!)


__global__ void roughness_multiple_scattering(Facet_Properties *FP, Roughness_Shadow *RS, Roughness_Illumination *RI, Roughness_Globalillumination *RGI, float *viewfactors, float solar, float scattering_accuracy, int step) {

	__shared__ bool iterate;
	__shared__ int iteratetest;
	__shared__ float albedo;
	__shared__ float powerold[SFCNUMFACE];

	int i;
	float sum;
	float shadow, illumangle, illumination, globalillumination, solarpower, powernew;
	
	if (threadIdx.x==0) {
		albedo=FP[blockIdx.x].albedo;
	}
	__syncthreads();
	
	shadow=RS[blockIdx.x].sfc_shadow[step][threadIdx.x];
	illumangle=RS[blockIdx.x].sfc_illumangle[step][threadIdx.x];
	illumination=RI[blockIdx.x].sfc_illumination[step][threadIdx.x];
	globalillumination=RGI[blockIdx.x].sfc_globalillumination[step][threadIdx.x];

	solarpower=(solar*illumangle*shadow)+globalillumination;
	powerold[threadIdx.x]=albedo*illumination;
	
	if (threadIdx.x==0) { iterate=true; }
	__syncthreads();

	while (iterate==true) {

		sum=0.0f;

		for (i=0;i!=SFCNUMFACE;++i) {
			sum=sum+(viewfactors[(i*SFCNUMFACE)+threadIdx.x]*powerold[i]);
		}

		powernew=albedo*(solarpower+sum);

		if (threadIdx.x==0) { iteratetest=0; }
		__syncthreads();

		if ( fabsf(powernew-powerold[threadIdx.x]) > scattering_accuracy ) { iteratetest=iteratetest+1; }
		powerold[threadIdx.x]=powernew;
		__syncthreads();
		
		if (threadIdx.x==0 && iteratetest==0) { iterate=false; }
		__syncthreads();

	}

	illumination=powernew/albedo;
	RI[blockIdx.x].sfc_illumination[step][threadIdx.x]=illumination;

};


// Zero roughness temperatures


__global__ void roughness_zero_temperatures(Roughness_Temperatures *RT) {

	int step, i;

	for (step=0;step!=NTAU;++step) {
		RT[blockIdx.x].sfc_stemp0[step][threadIdx.x]=0.0f;
		RT[blockIdx.x].sfc_stemp1[step][threadIdx.x]=0.0f;
	}
	for (i=0;i!=(NX+1);++i) {
		RT[blockIdx.x].sfc_temp[i][threadIdx.x]=0.0f;
	}

};


// Initialise roughness temperatures to begin thermal modelling


__global__ void roughness_initialise_temperatures(Facet_Properties *FP, Roughness_Illumination *RI, Roughness_Globalheating *RGH, Roughness_Temperatures *RT, float *viewfactors, float EPSSIGMA) {

	__shared__ float albedo;
	__shared__ float stemp0[SFCNUMFACE];

	int step, i;
	float temp, selfheatingtemp;

	if (threadIdx.x==0) {
		albedo=FP[blockIdx.x].albedo;
	}
	__syncthreads();

	// First initialisation step

	temp=0.0f;
	for (step=0;step!=NTAU;++step) {
		temp=temp+((1.0f/NTAU)*(((1.0f-albedo)*RI[blockIdx.x].sfc_illumination[step][threadIdx.x])+RGH[blockIdx.x].sfc_globalheating[step][threadIdx.x]));
	}
	stemp0[threadIdx.x]=powf(temp/(EPSSIGMA),0.25f);
	__syncthreads();

	// Second initialisation step

	selfheatingtemp=0.0f;
	for (i=0;i!=SFCNUMFACE;++i) {
		selfheatingtemp=selfheatingtemp+(EPSSIGMA*powf(stemp0[i],4.0f)*viewfactors[(i*SFCNUMFACE)+threadIdx.x]);
	}

	temp=0.0f;
	for (step=0;step!=NTAU;++step) {
		temp=temp+((1.0f/NTAU)*(((1.0f-albedo)*RI[blockIdx.x].sfc_illumination[step][threadIdx.x])+RGH[blockIdx.x].sfc_globalheating[step][threadIdx.x]+selfheatingtemp));
	}
	temp=powf(temp/(EPSSIGMA),0.25f);
	
	for (i=0;i!=(NX+1);++i) {
		RT[blockIdx.x].sfc_temp[i][threadIdx.x]=temp;
	}
	for (step=0;step!=NTAU;++step) {
		RT[blockIdx.x].sfc_stemp0[step][threadIdx.x]=temp;
		RT[blockIdx.x].sfc_stemp1[step][threadIdx.x]=temp;
	}

};


// Finite difference


__global__ void finite_difference(Facet_Properties *FP, Roughness_Temperatures *RT) {

	__shared__ bool converged;

	int i;
	float DZ, DTAU, MULTI0, temp[3];

	if (threadIdx.x==0) {
		converged=FP[blockIdx.x].sfc_converged;
	}
	__syncthreads();

	if (converged==false) {

		DZ=8.0f/NX;
		DTAU=(2.0f*PI)/NTAU;
		MULTI0=DTAU/(DZ*DZ);

		temp[0]=RT[blockIdx.x].sfc_temp[0][threadIdx.x];
		temp[1]=RT[blockIdx.x].sfc_temp[1][threadIdx.x];
		temp[2]=RT[blockIdx.x].sfc_temp[2][threadIdx.x];

		for (i=1;i!=NX;++i) {
			RT[blockIdx.x].sfc_temp[i][threadIdx.x]=temp[1]+(MULTI0*(temp[2]-(2.0f*temp[1])+temp[0]));
			temp[0]=temp[1];
			temp[1]=temp[2];
			if (i<=NX) { temp[2]=RT[blockIdx.x].sfc_temp[i+2][threadIdx.x]; }
		}

		RT[blockIdx.x].sfc_temp[NX][threadIdx.x]=RT[blockIdx.x].sfc_temp[NX-1][threadIdx.x];

	}

	__syncthreads();
	
};


// Boundary condition


__global__ void boundary_condition(Facet_Properties *FP, Roughness_Illumination *RI, Roughness_Globalheating *RGH, Roughness_Temperatures *RT, float *viewfactors, float EPSSIGMA, float TACC, int step) {

	__shared__ bool converged;
	__shared__ float albedo;
	__shared__ float thermal_inertia;
	__shared__ float stemp[SFCNUMFACE];

	int i;
	float DZ, MULTI1;
	float MULTI2;
	float SBC, SBCD;
	float DT, TR[2], temp;
	float incidentflux;
	float illumination;

	if (threadIdx.x==0) {
		converged=FP[blockIdx.x].sfc_converged;
	}
	__syncthreads();

	if (converged==false) {

		// Load facet properties

		if (threadIdx.x==0) {
			albedo=FP[blockIdx.x].albedo;
			thermal_inertia=FP[blockIdx.x].thermal_inertia;
		}
		__syncthreads();

		// Load surface temperature and illumination data

		TR[0]=RT[blockIdx.x].sfc_temp[0][threadIdx.x];
		temp=RT[blockIdx.x].sfc_temp[1][threadIdx.x];
		stemp[threadIdx.x]=RT[blockIdx.x].sfc_stemp1[step][threadIdx.x];
		RT[blockIdx.x].sfc_stemp0[step][threadIdx.x]=stemp[threadIdx.x];
		stemp[threadIdx.x]=EPSSIGMA*powf(stemp[threadIdx.x],4.0f);
		illumination=RI[blockIdx.x].sfc_illumination[step][threadIdx.x];
		incidentflux=RGH[blockIdx.x].sfc_globalheating[step][threadIdx.x];
		__syncthreads();

		// Determine incident flux by summing over neighbouring facets

		for (i=0;i!=SFCNUMFACE;++i) {
			incidentflux=incidentflux+(stemp[i]*viewfactors[(i*SFCNUMFACE)+threadIdx.x]);
		}

		// Apply surface boundary condition
		
		DT=1.0f;
		DZ=8.0f/NX;
		MULTI1=thermal_inertia*sqrtf(2.0f*PI/rot_period)/DZ;
		MULTI2=((1.0f-albedo)*illumination)+incidentflux;

		while (DT >= TACC) {
			SBC=MULTI2+(MULTI1*(temp-TR[0]))-(EPSSIGMA*TR[0]*TR[0]*TR[0]*TR[0]);
			SBCD=MULTI1+(4.0f*EPSSIGMA*TR[0]*TR[0]*TR[0]);
			TR[1]=TR[0]+(SBC/SBCD);
			DT=fabsf(TR[1]-TR[0]);
			TR[0]=TR[1];
		}
		__syncthreads();

		// Write updated surface temperature data

		RT[blockIdx.x].sfc_temp[0][threadIdx.x]=TR[1];
		RT[blockIdx.x].sfc_stemp1[step][threadIdx.x]=TR[1];
			
	}

	__syncthreads();

};


// Instantaneous equilibrium thermal model


__global__ void roughness_inst_equilibrium(Facet_Properties *FP, Roughness_Illumination *RI, Roughness_Globalheating *RGH, Roughness_Temperatures *RT, float *viewfactors, float EPSSIGMA, int step) {

	__shared__ bool converged;
	__shared__ float albedo;
	__shared__ float stemp[SFCNUMFACE];
	
	int i;
	float incidentflux;
	float illumination;

	if (threadIdx.x==0) {
		converged=FP[blockIdx.x].sfc_converged;
	}
	__syncthreads();

	if (converged==false) {

		// Load facet properties

		if (threadIdx.x==0) {
			albedo=FP[blockIdx.x].albedo;
		}
		__syncthreads();

		// Load surface temperature and illumination data

		stemp[threadIdx.x]=RT[blockIdx.x].sfc_stemp1[step][threadIdx.x];
		RT[blockIdx.x].sfc_stemp0[step][threadIdx.x]=stemp[threadIdx.x];
		stemp[threadIdx.x]=EPSSIGMA*powf(stemp[threadIdx.x],4.0f);
		illumination=RI[blockIdx.x].sfc_illumination[step][threadIdx.x];
		incidentflux=RGH[blockIdx.x].sfc_globalheating[step][threadIdx.x];
		__syncthreads();

		// Determine incident flux by summing over neighbouring facets

		for (i=0;i!=SFCNUMFACE;++i) {
			incidentflux=incidentflux+(stemp[i]*viewfactors[(i*SFCNUMFACE)+threadIdx.x]);
		}
		__syncthreads();

		// Instantaneous radiative equilibrium

		stemp[threadIdx.x]=powf(((((1.0f-albedo)*illumination)+incidentflux)/(EPSSIGMA)),0.25f);

		// Write updated surface temperature data

		RT[blockIdx.x].sfc_stemp1[step][threadIdx.x]=stemp[threadIdx.x];

	}

	__syncthreads();

};


// Code to test that asteroid surface roughness temperatures have converged


__global__ void converged(Facet_Properties *FP, Roughness_Temperatures *RT, float TACC) {

	__shared__ bool converged;

	int step;
	float stemp0;
	float stemp1;
	float DT;

	if (threadIdx.x==0) {
		converged=FP[blockIdx.x].sfc_converged;
	}
	__syncthreads();

	if (converged==false) {

		__syncthreads();
		if (threadIdx.x==0) {
			converged=true;
		}
		__syncthreads();

		for (step=0;step!=NTAU;++step) {
			
			stemp0=RT[blockIdx.x].sfc_stemp0[step][threadIdx.x];
			stemp1=RT[blockIdx.x].sfc_stemp1[step][threadIdx.x];

			DT=fabsf(stemp0-stemp1);
			if (DT >= TACC) { converged=false; }

		}

	}
	__syncthreads();

	if (threadIdx.x==0) {
		FP[blockIdx.x].sfc_converged=converged;
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
void initialise_GPU_constants(float rot_period_H, float viewfactors_H[SFCNUMFACE][SFCNUMFACE], GPU_Surface S_H) {

	float PI_H = 4.0f*atan(1.0f);

	cudaMemcpyToSymbol(PI,&PI_H,sizeof(float));
	cudaMemcpyToSymbol(rot_period,&rot_period_H,sizeof(float));

	cudaMalloc((void **) &viewfactors_D,SFCNUMFACE*SFCNUMFACE*sizeof(float));
	cudaMemcpy(viewfactors_D,viewfactors_H,SFCNUMFACE*SFCNUMFACE*sizeof(float),cudaMemcpyHostToDevice);

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


// Code to initialise GPU facet properties arrays


extern "C"
void initialise_GPU_facet_properties_arrays() {

	size_t size;
	size=ASTNUMFACE*sizeof(Facet_Properties);
	cudaMalloc((void **) &FP_D, size);

};


// Code to initialise GPU roughness shadow arrays


extern "C"
void initialise_GPU_shadow_arrays() {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Shadow);
	cudaMalloc((void **) &RS_D, size);

};


// Code to initialise GPU roughness illumination arrays


extern "C"
void initialise_GPU_illumination_arrays() {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Illumination);
	cudaMalloc((void **) &RI_D, size);

};


// Code to initialise GPU roughness globalillumination arrays


extern "C"
void initialise_GPU_globalillumination_arrays() {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Illumination);
	cudaMalloc((void **) &RGI_D, size);

};


// Code to initialise GPU roughness globalheating arrays


extern "C"
void initialise_GPU_globalheating_arrays() {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Globalheating);
	cudaMalloc((void **) &RGH_D, size);

};


// Code to initialise GPU roughness temperature arrays


extern "C"
void initialise_GPU_temperature_arrays() {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Temperatures);
	cudaMalloc((void **) &RT_D, size);

};


// Code to free up GPU roughness shadow arrays


extern "C"
void freeup_GPU_shadow_arrays() {

	cudaFree(RS_D);

};


// Code to free up GPU roughness globalillumination arrays


extern "C"
void freeup_GPU_globalillumination_arrays() {

	cudaFree(RGI_D);

};


// Code to free up GPU roughness globalheating arrays


extern "C"
void freeup_GPU_globalheating_arrays() {

	cudaFree(RGH_D);

};


// Code to free up GPU roughness temperature arrays


extern "C"
void freeup_GPU_temperature_arrays() {

	cudaFree(RT_D);

};


// Code to load GPU facet properties data


extern "C"
void set_GPU_facet_properties(Facet_Properties *FP_H) {

	size_t size;
	size=ASTNUMFACE*sizeof(Facet_Properties);
	cudaMemcpy(FP_D,FP_H,size,cudaMemcpyHostToDevice);

};


// Code to load GPU roughness global-illumination data


extern "C"
void load_GPU_globalillumination(Roughness_Globalillumination *RGI_H) {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Globalillumination);
	cudaMemcpy(RGI_D,RGI_H,size,cudaMemcpyHostToDevice);

};


// Code to load GPU roughness global-selfheating data


extern "C"
void load_GPU_globalheating(Roughness_Globalheating *RGH_H) {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Globalheating);
	cudaMemcpy(RGH_D,RGH_H,size,cudaMemcpyHostToDevice);

};


// Code to obtain GPU facet properties data


extern "C"
void obtain_GPU_facet_properties(Facet_Properties *FP_H) {

	size_t size;
	size=ASTNUMFACE*sizeof(Facet_Properties);
	cudaMemcpy(FP_H,FP_D,size,cudaMemcpyDeviceToHost);

};


// Code to obtain GPU roughness shadow data


extern "C" 
void obtain_GPU_roughness_shadows(Roughness_Shadow *RS_H) {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Shadow);
	cudaMemcpy(RS_H,RS_D,size,cudaMemcpyDeviceToHost);

};


// Code to obtain GPU roughness temperature data


extern "C" 
void obtain_GPU_roughness_temperatures(Roughness_Temperatures *RT_H) {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Temperatures);
	cudaMemcpy(RT_H,RT_D,size,cudaMemcpyDeviceToHost);

};


// Code to obtain GPU roughness illumination data


extern "C" 
void obtain_GPU_roughness_illumination(Roughness_Illumination *RI_H) {

	size_t size;
	size=ASTNUMFACE*sizeof(Roughness_Illumination);
	cudaMemcpy(RI_H,RI_D,size,cudaMemcpyDeviceToHost);

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


// Surface roughness single scattering GPU control code


extern "C"
void GPU_single_scattering(float solar) {

	for (int step=0;step!=NTAU;++step) {	
		roughness_single_scattering <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RS_D,RI_D,RGI_D,viewfactors_D,solar,step);
		cudaDeviceSynchronize();
	}

};


// Surface roughness multiple scattering GPU control code


extern "C"
void GPU_multiple_scattering(float solar, float scattering_accuracy) {

	for (int step=0;step!=NTAU;++step) {	
		roughness_multiple_scattering <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RS_D,RI_D,RGI_D,viewfactors_D,solar,scattering_accuracy,step);
		cudaDeviceSynchronize();
	}

};


// Surface roughness initialise temperatures GPU control code


extern "C"
void GPU_zero_temperatures() {

	roughness_zero_temperatures <<<ASTNUMFACE,SFCNUMFACE>>> (RT_D);
	cudaDeviceSynchronize();

};


// Surface roughness initialise temperatures GPU control code


extern "C"
void GPU_initialise_temperatures(float EPS, float SIGMA) {

	float EPSSIGMA;
	EPSSIGMA=EPS*SIGMA;

	roughness_initialise_temperatures <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RI_D,RGH_D,RT_D,viewfactors_D,EPSSIGMA);
	cudaDeviceSynchronize();

};


// Surface roughness finite difference GPU control code


extern "C"
void GPU_inertia_difference(float EPS, float SIGMA, float TACC) {

	float EPSSIGMA;
	EPSSIGMA=EPS*SIGMA;

	float TACC2;
	TACC2=TACC/10.0f;

	for (int step=0;step!=NTAU;++step) {
		finite_difference <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RT_D);
		cudaDeviceSynchronize();
		boundary_condition <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RI_D,RGH_D,RT_D,viewfactors_D,EPSSIGMA,TACC2,step);
		cudaDeviceSynchronize();
	}

};


// Surface roughness instantaneous equilibrium GPU control code


extern "C"
void GPU_inst_equilibrium(float EPS, float SIGMA) {

	float EPSSIGMA;
	EPSSIGMA=EPS*SIGMA;

	for (int step=0;step!=NTAU;++step) {
		roughness_inst_equilibrium <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RI_D,RGH_D,RT_D,viewfactors_D,EPSSIGMA,step);
		cudaDeviceSynchronize();
	}

};


// Temperatures converged GPU control code


extern "C"
void GPU_converged(float TACC) {

	converged <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RT_D,TACC);
	cudaDeviceSynchronize();

};