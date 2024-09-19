#include "cuda_runtime.h"
#include "device_launch_parameters.h"

const int ASTNUMFACE = 512;												// Number of shape model facets
const int SFCNUMVTX = 61;												// Number of surface roughness model vertices
const int SFCNUMFACE = 100;												// Number of surface roughness model facets
const int NUMPRIREV = 9;												
const int NTAU = NUMPRIREV*650;											// Number of model time steps (multiples of 50)
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


//struct Facet_Properties {

//	bool sfc_converged;
//	float albedo;
//	float thermal_inertia;

//};

bool *sfc_converged_D;
float *albedo_D;
float *thermal_inertia_D;
//float *rotperiod_D;
float *numrevs_D;

// Structure Roughness_Shadow containing roughness facet shadow and illumination angle information


//struct Roughness_Shadow {

//	int shadow[NTAU];
//	float sfc_shadow[NTAU][SFCNUMFACE];
//	float sfc_illumvector[NTAU][3];
//	float sfc_illumangle[NTAU][SFCNUMFACE];

//};

int *shadow_D;
float *sfc_shadow_D;
float *sfc_illumvector_D;
float *sfc_illumangle_D;


// Structure Roughness_Illumination containing roughness facet illumination flux information


//struct Roughness_Illumination {

//	float sfc_illumination[NTAU][SFCNUMFACE];

//};

float *sfc_illumination_D;


// Structure Roughness_Globalillumination containing roughness global-illumination flux information


//struct Roughness_Globalillumination {

//	float sfc_globalillumination[NTAU][SFCNUMFACE];

//};

float *sfc_globalillumination_D;


// Structure Roughness_Globalheating containing roughness facet global-selfheating flux information


//struct Roughness_Globalheating {

//	float sfc_globalheating[NTAU][SFCNUMFACE];

//};

float *sfc_globalheating_D;


// Structure Roughness_Temperatures containing roughness facet temperature information


//struct Roughness_Temperatures {

//	float sfc_stemp0[NTAU][SFCNUMFACE];
//	float sfc_stemp1[NTAU][SFCNUMFACE];
//	float sfc_temp[NX+1][SFCNUMFACE];

//};

float *sfc_stemp0_D;
float *sfc_stemp1_D;
float *sfc_temp_D;


// Pointers to GPU device memory


//Facet_Properties *FP_D;
//Roughness_Shadow *RS_D;
//Roughness_Illumination *RI_D;
//Roughness_Globalillumination *RGI_D;
//Roughness_Globalheating *RGH_D;
//Roughness_Temperatures *RT_D;
float *viewfactors_D;


// Code to determine surface partial shadowing


//__global__ void sfcpartialshadow(Roughness_Shadow *RS_D, int step) {
__global__ void sfcpartialshadow(int *shadow_D, float *sfc_shadow_D, float *sfc_illumvector_D, float *sfc_illumangle_D, int step) {

	__shared__ int shadow;
	__shared__ float sunvector[3];
	__shared__ float sfcshadow[SFCNUMFACE];

	int p;
	float midpos[3];
	float trans_co[2];
	float a, s;
	float shadow_radius;
	float crater_radius = 10.0f;

	if (threadIdx.x == 0) {
		//	shadow=RS_D[blockIdx.x].shadow[step];
		shadow = shadow_D[(blockIdx.x*NTAU) + step];
	}
	__syncthreads();

	if (shadow == 1) {

		if (threadIdx.x<4) {
			//	sunvector[threadIdx.x]=RS_D[blockIdx.x].sfc_illumvector[step][threadIdx.x];
			sunvector[threadIdx.x] = sfc_illumvector_D[(blockIdx.x*NTAU * 3) + (step * 3) + threadIdx.x];
		}
		//	sfcshadow[threadIdx.x]=RS_D[blockIdx.x].sfc_shadow[step][threadIdx.x];
		sfcshadow[threadIdx.x] = sfc_shadow_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
		__syncthreads();

		if (threadIdx.x == 0) {
			s = sqrtf((sunvector[0] * sunvector[0]) + (sunvector[1] * sunvector[1]));
			sunvector[0] = sunvector[0] / s;
			sunvector[1] = sunvector[1] / s;
			sunvector[2] = sunvector[2] / s;
		}
		__syncthreads();

		a = fabsf(sunvector[2]) / sqrtf((sunvector[0] * sunvector[0]) + (sunvector[1] * sunvector[1]));

		for (p = 0; p != (m*m); ++p) {

			midpos[0] = sfc_vertex[sfc_facet[threadIdx.x][2]][0] + (trivector_multiplier[p][0] * sfc_trivector[threadIdx.x][0][0]) + (trivector_multiplier[p][1] * sfc_trivector[threadIdx.x][1][0]);
			midpos[1] = sfc_vertex[sfc_facet[threadIdx.x][2]][1] + (trivector_multiplier[p][0] * sfc_trivector[threadIdx.x][0][1]) + (trivector_multiplier[p][1] * sfc_trivector[threadIdx.x][1][1]);
			midpos[2] = sfc_vertex[sfc_facet[threadIdx.x][2]][2] + (trivector_multiplier[p][0] * sfc_trivector[threadIdx.x][0][2]) + (trivector_multiplier[p][1] * sfc_trivector[threadIdx.x][1][2]);

			trans_co[1] = (-1.0f*sunvector[1] * midpos[0]) + (sunvector[0] * midpos[1]);
			shadow_radius = sqrtf((crater_radius*crater_radius) - (trans_co[1] * trans_co[1]));
			trans_co[0] = (sunvector[0] * midpos[0]) + (sunvector[1] * midpos[1]) - shadow_radius;

			if (midpos[2] <= (a*trans_co[0]) && sfcshadow[threadIdx.x] != 0.0f) { sfcshadow[threadIdx.x] = sfcshadow[threadIdx.x] - (1.0f / (float)(m*m)); }
			if (sfcshadow[threadIdx.x]<0.0f) { sfcshadow[threadIdx.x] = 0.0f; }

		}
		__syncthreads();

		//	RS_D[blockIdx.x].sfc_shadow[step][threadIdx.x]=sfcshadow[threadIdx.x];
		sfc_shadow_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = sfcshadow[threadIdx.x];

	}

	__syncthreads();

};


// Determine roughness facet illumination from direct and single scattered sunlight (not including thermal radiation!)


//__global__ void roughness_single_scattering(Facet_Properties *FP, Roughness_Shadow *RS, Roughness_Illumination *RI, Roughness_Globalillumination *RGI, float *viewfactors, float solar, int step) {
__global__ void roughness_single_scattering(float *albedo_D, float *sfc_shadow_D, float *sfc_illumangle_D, float *sfc_globalillumination_D, float *sfc_illumination_D, float *viewfactors, float solar, int step) {

	__shared__ float albedo;
	__shared__ float illumination[SFCNUMFACE];

	int i;
	float shadow, illumangle, globalillumination, totalillumination;

	if (threadIdx.x == 0) {
		//	albedo=FP[blockIdx.x].albedo;
		albedo = albedo_D[blockIdx.x];
	}
	__syncthreads();

	//	shadow=RS[blockIdx.x].sfc_shadow[step][threadIdx.x];
	//	illumangle=RS[blockIdx.x].sfc_illumangle[step][threadIdx.x];
	//	globalillumination=RGI[blockIdx.x].sfc_globalillumination[step][threadIdx.x];
	shadow = sfc_shadow_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
	illumangle = sfc_illumangle_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
	globalillumination = sfc_globalillumination_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
	__syncthreads();

	illumination[threadIdx.x] = (solar*shadow*illumangle) + globalillumination;
	totalillumination = illumination[threadIdx.x];

	for (i = 0; i != SFCNUMFACE; ++i) {
		totalillumination = totalillumination + (viewfactors[(i*SFCNUMFACE) + threadIdx.x] * albedo*illumination[i]);
	}

	if (totalillumination <= 0.0f) { totalillumination = 0.0f; }
	__syncthreads();

	//	RI[blockIdx.x].sfc_illumination[step][threadIdx.x]=totalillumination;
	sfc_illumination_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = totalillumination;

};


// Determine roughness facet illumination from direct and multiple scattered sunlight (not including thermal radiation!)


//__global__ void roughness_multiple_scattering(Facet_Properties *FP, Roughness_Shadow *RS, Roughness_Illumination *RI, Roughness_Globalillumination *RGI, float *viewfactors, float solar, float scattering_accuracy, int step) {
__global__ void roughness_multiple_scattering(float *albedo_D, float *sfc_shadow_D, float *sfc_illumangle_D, float *sfc_illumination_D, float *sfc_globalillumination_D, float *viewfactors, float solar, float scattering_accuracy, int step) {

	__shared__ bool iterate;
	__shared__ int iteratetest;
	__shared__ float albedo;
	__shared__ float powerold[SFCNUMFACE];

	int i;
	float sum;
	float shadow, illumangle, illumination, globalillumination, solarpower, powernew;

	if (threadIdx.x == 0) {
		//	albedo=FP[blockIdx.x].albedo;
		albedo = albedo_D[blockIdx.x];
	}
	__syncthreads();

	//	shadow=RS[blockIdx.x].sfc_shadow[step][threadIdx.x];
	//	illumangle=RS[blockIdx.x].sfc_illumangle[step][threadIdx.x];
	//	illumination=RI[blockIdx.x].sfc_illumination[step][threadIdx.x];
	//	globalillumination=RGI[blockIdx.x].sfc_globalillumination[step][threadIdx.x];
	shadow = sfc_shadow_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
	illumangle = sfc_illumangle_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
	illumination = sfc_illumination_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
	globalillumination = sfc_globalillumination_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];

	solarpower = (solar*illumangle*shadow) + globalillumination;
	powerold[threadIdx.x] = albedo*illumination;

	if (threadIdx.x == 0) { iterate = true; }
	__syncthreads();

	while (iterate == true) {

		sum = 0.0f;

		for (i = 0; i != SFCNUMFACE; ++i) {
			sum = sum + (viewfactors[(i*SFCNUMFACE) + threadIdx.x] * powerold[i]);
		}

		powernew = albedo*(solarpower + sum);

		if (threadIdx.x == 0) { iteratetest = 0; }
		__syncthreads();

		if (fabsf(powernew - powerold[threadIdx.x]) > scattering_accuracy) { iteratetest = iteratetest + 1; }
		powerold[threadIdx.x] = powernew;
		__syncthreads();

		if (threadIdx.x == 0 && iteratetest == 0) { iterate = false; }
		__syncthreads();

	}

	illumination = powernew / albedo;
	//	RI[blockIdx.x].sfc_illumination[step][threadIdx.x]=illumination;
	sfc_illumination_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = illumination;

};


// Zero roughness temperatures


//__global__ void roughness_zero_temperatures(Roughness_Temperatures *RT) {
__global__ void roughness_zero_temperatures(float *sfc_stemp0_D, float *sfc_stemp1_D, float *sfc_temp_D) {

	int step, i;

	//	for (step=0;step!=NTAU;++step) {
	//		RT[blockIdx.x].sfc_stemp0[step][threadIdx.x]=0.0f;
	//		RT[blockIdx.x].sfc_stemp1[step][threadIdx.x]=0.0f;
	//	}
	//	for (i=0;i!=(NX+1);++i) {
	//		RT[blockIdx.x].sfc_temp[i][threadIdx.x]=0.0f;
	//	}
	for (step = 0; step != NTAU; ++step) {
		sfc_stemp0_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = 0.0f;
		sfc_stemp1_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = 0.0f;
	}
	for (i = 0; i != (NX + 1); ++i) {
		sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + (i*SFCNUMFACE) + threadIdx.x] = 0.0f;
	}

};


// Initialise roughness temperatures to begin thermal modelling


//__global__ void roughness_initialise_temperatures(Facet_Properties *FP, Roughness_Illumination *RI, Roughness_Globalheating *RGH, Roughness_Temperatures *RT, float *viewfactors, float EPSSIGMA) {
__global__ void roughness_initialise_temperatures(float *albedo_D, float *sfc_illumination_D, float *sfc_globalheating_D, float *sfc_temp_D, float *sfc_stemp0_D, float *sfc_stemp1_D, float *viewfactors, float EPSSIGMA, float IR_ABS) {

	__shared__ float albedo;
	__shared__ float stemp0[SFCNUMFACE];

	int step, i;
	float temp, selfheatingtemp;

	if (threadIdx.x == 0) {
		//	albedo=FP[blockIdx.x].albedo;
		albedo = albedo_D[blockIdx.x];
	}
	__syncthreads();

	// First initialisation step

	temp = 0.0f;
	for (step = 0; step != NTAU; ++step) {
		//	temp=temp+((1.0f/NTAU)*(((1.0f-albedo)*RI[blockIdx.x].sfc_illumination[step][threadIdx.x])+RGH[blockIdx.x].sfc_globalheating[step][threadIdx.x]));
		temp = temp + ((1.0f / NTAU)*(((1.0f - albedo)*sfc_illumination_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x]) + sfc_globalheating_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x]));
	}
	stemp0[threadIdx.x] = powf(temp / (EPSSIGMA), 0.25f);
	__syncthreads();

	// Second initialisation step

	selfheatingtemp = 0.0f;
	for (i = 0; i != SFCNUMFACE; ++i) {
		selfheatingtemp = selfheatingtemp + (IR_ABS*EPSSIGMA*powf(stemp0[i], 4.0f)*viewfactors[(i*SFCNUMFACE) + threadIdx.x]);
	}

	temp = 0.0f;
	for (step = 0; step != NTAU; ++step) {
		//	temp=temp+((1.0f/NTAU)*(((1.0f-albedo)*RI[blockIdx.x].sfc_illumination[step][threadIdx.x])+RGH[blockIdx.x].sfc_globalheating[step][threadIdx.x]+selfheatingtemp));
		temp = temp + ((1.0f / NTAU)*(((1.0f - albedo)*sfc_illumination_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x]) + sfc_globalheating_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] + selfheatingtemp));
	}
	temp = powf(temp / (EPSSIGMA), 0.25f);

	for (i = 0; i != (NX + 1); ++i) {
		//	RT[blockIdx.x].sfc_temp[i][threadIdx.x]=temp;
		sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + (i*SFCNUMFACE) + threadIdx.x] = temp;
	}
	for (step = 0; step != NTAU; ++step) {
		//	RT[blockIdx.x].sfc_stemp0[step][threadIdx.x]=temp;
		//	RT[blockIdx.x].sfc_stemp1[step][threadIdx.x]=temp;
		sfc_stemp0_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = temp;
		sfc_stemp1_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = temp;
	}

};


// Finite difference


//__global__ void finite_difference(Facet_Properties *FP, Roughness_Temperatures *RT) {
__global__ void finite_difference(bool *sfc_converged_D, float *numrevs_D, float *sfc_temp_D) {

	__shared__ bool converged;
	__shared__ float numrevs;

	int i;
	float DZ, DTAU, MULTI0, temp[3];

	if (threadIdx.x == 0) {
		//	converged=FP[blockIdx.x].sfc_converged;
		converged = sfc_converged_D[blockIdx.x];
	}
	__syncthreads();

	if (converged == false) {

		if (threadIdx.x == 0) {
			numrevs = numrevs_D[blockIdx.x];
		}
		__syncthreads();

		DZ = 8.0f / (float)NX;
		DTAU = (2.0f*PI*numrevs) / (float)NTAU;
		MULTI0 = DTAU / (DZ*DZ);

		//	temp[0]=RT[blockIdx.x].sfc_temp[0][threadIdx.x];
		//	temp[1]=RT[blockIdx.x].sfc_temp[1][threadIdx.x];
		//	temp[2]=RT[blockIdx.x].sfc_temp[2][threadIdx.x];
		temp[0] = sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + (0 * SFCNUMFACE) + threadIdx.x];
		temp[1] = sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + (1 * SFCNUMFACE) + threadIdx.x];
		temp[2] = sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + (2 * SFCNUMFACE) + threadIdx.x];

		for (i = 1; i != NX; ++i) {
			//	RT[blockIdx.x].sfc_temp[i][threadIdx.x]=temp[1]+(MULTI0*(temp[2]-(2.0f*temp[1])+temp[0]));
			sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + (i*SFCNUMFACE) + threadIdx.x] = temp[1] + (MULTI0*(temp[2] - (2.0f*temp[1]) + temp[0]));
			temp[0] = temp[1];
			temp[1] = temp[2];
			//	if (i<=NX) { temp[2]=RT[blockIdx.x].sfc_temp[i+2][threadIdx.x]; }
			if (i <= NX) { temp[2] = sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + ((i + 2)*SFCNUMFACE) + threadIdx.x]; }
		}

		//	RT[blockIdx.x].sfc_temp[NX][threadIdx.x]=RT[blockIdx.x].sfc_temp[NX-1][threadIdx.x];
		sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + (NX*SFCNUMFACE) + threadIdx.x] = sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + ((NX - 1)*SFCNUMFACE) + threadIdx.x];

	}

	__syncthreads();

};


// Boundary condition


//__global__ void boundary_condition(Facet_Properties *FP, Roughness_Illumination *RI, Roughness_Globalheating *RGH, Roughness_Temperatures *RT, float *viewfactors, float EPSSIGMA, float TACC, int step) {
__global__ void boundary_condition(bool *sfc_converged_D, float *albedo_D, float *thermal_inertia_D, float *numrevs_D, float *sfc_temp_D, float *sfc_stemp0_D, float *sfc_stemp1_D, float *sfc_illumination_D, float *sfc_globalheating_D, float *viewfactors, float EPSSIGMA, float IR_ABS, float TACC, int step) {

	__shared__ bool converged;
	__shared__ float albedo;
	__shared__ float thermal_inertia;
//	__shared__ float rotperiod;
	__shared__ float numrevs;
	__shared__ float stemp[SFCNUMFACE];

	int i;
	float DZ, MULTI1;
	float MULTI2;
	float SBC, SBCD;
	float DT, TR[2], temp;
	float incidentflux;
	float illumination;

	if (threadIdx.x == 0) {
		//	converged=FP[blockIdx.x].sfc_converged;
		converged = sfc_converged_D[blockIdx.x];
	}
	__syncthreads();

	if (converged == false) {

		// Load facet properties

		if (threadIdx.x == 0) {
			//	albedo=FP[blockIdx.x].albedo;
			//	thermal_inertia=FP[blockIdx.x].thermal_inertia;
			albedo = albedo_D[blockIdx.x];
			thermal_inertia = thermal_inertia_D[blockIdx.x];
			// rotperiod = rotperiod_D[blockIdx.x];
			numrevs = numrevs_D[blockIdx.x];
		}
		__syncthreads();

		// Load surface temperature and illumination data

		//	TR[0]=RT[blockIdx.x].sfc_temp[0][threadIdx.x];
		//	temp=RT[blockIdx.x].sfc_temp[1][threadIdx.x];
		//	stemp[threadIdx.x]=RT[blockIdx.x].sfc_stemp1[step][threadIdx.x];
		//	RT[blockIdx.x].sfc_stemp0[step][threadIdx.x]=stemp[threadIdx.x];
		//	stemp[threadIdx.x]=EPSSIGMA*powf(stemp[threadIdx.x],4.0f);
		//	illumination=RI[blockIdx.x].sfc_illumination[step][threadIdx.x];
		//	incidentflux=RGH[blockIdx.x].sfc_globalheating[step][threadIdx.x];
		TR[0] = sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + (0 * SFCNUMFACE) + threadIdx.x];
		temp = sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + (1 * SFCNUMFACE) + threadIdx.x];
		stemp[threadIdx.x] = sfc_stemp1_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
		sfc_stemp0_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = stemp[threadIdx.x];
		stemp[threadIdx.x] = IR_ABS*EPSSIGMA*powf(stemp[threadIdx.x], 4.0f);
		illumination = sfc_illumination_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
		incidentflux = sfc_globalheating_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
		__syncthreads();

		// Determine incident flux by summing over neighbouring facets

		for (i = 0; i != SFCNUMFACE; ++i) {
			incidentflux = incidentflux + (stemp[i] * viewfactors[(i*SFCNUMFACE) + threadIdx.x]);
		}

		// Apply surface boundary condition

		DT = 1.0f;
		DZ = 8.0f / (float)NX;
//		MULTI1 = thermal_inertia*sqrtf(2.0f*PI / rot_period) / DZ;
		MULTI1 = thermal_inertia*sqrtf(2.0f*PI*numrevs / rot_period) / DZ;
		MULTI2 = ((1.0f - albedo)*illumination) + incidentflux;

		while (DT >= TACC) {
			SBC = MULTI2 + (MULTI1*(temp - TR[0])) - (EPSSIGMA*TR[0] * TR[0] * TR[0] * TR[0]);
			SBCD = MULTI1 + (4.0f*EPSSIGMA*TR[0] * TR[0] * TR[0]);
			TR[1] = TR[0] + (SBC / SBCD);
			DT = fabsf(TR[1] - TR[0]);
			TR[0] = TR[1];
		}
		__syncthreads();

		// Write updated surface temperature data

		//	RT[blockIdx.x].sfc_temp[0][threadIdx.x]=TR[1];
		//	RT[blockIdx.x].sfc_stemp1[step][threadIdx.x]=TR[1];
		sfc_temp_D[(blockIdx.x*(NX + 1)*SFCNUMFACE) + (0 * SFCNUMFACE) + threadIdx.x] = TR[1];
		sfc_stemp1_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = TR[1];

	}

	__syncthreads();

};


// Instantaneous equilibrium thermal model


//__global__ void roughness_inst_equilibrium(Facet_Properties *FP, Roughness_Illumination *RI, Roughness_Globalheating *RGH, Roughness_Temperatures *RT, float *viewfactors, float EPSSIGMA, int step) {
__global__ void roughness_inst_equilibrium(bool *sfc_converged_D, float *albedo_D, float *sfc_stemp0_D, float *sfc_stemp1_D, float *sfc_illumination_D, float *sfc_globalheating_D, float *viewfactors, float EPSSIGMA, float IR_ABS, int step) {

	__shared__ bool converged;
	__shared__ float albedo;
	__shared__ float stemp[SFCNUMFACE];

	int i;
	float incidentflux;
	float illumination;

	if (threadIdx.x == 0) {
		//	converged=FP[blockIdx.x].sfc_converged;
		converged = sfc_converged_D[blockIdx.x];
	}
	__syncthreads();

	if (converged == false) {

		// Load facet properties

		if (threadIdx.x == 0) {
			//	albedo=FP[blockIdx.x].albedo;
			albedo = albedo_D[blockIdx.x];
		}
		__syncthreads();

		// Load surface temperature and illumination data

		//	stemp[threadIdx.x]=RT[blockIdx.x].sfc_stemp1[step][threadIdx.x];
		//	RT[blockIdx.x].sfc_stemp0[step][threadIdx.x]=stemp[threadIdx.x];
		//	stemp[threadIdx.x]=EPSSIGMA*powf(stemp[threadIdx.x],4.0f);
		//	illumination=RI[blockIdx.x].sfc_illumination[step][threadIdx.x];
		//	incidentflux=RGH[blockIdx.x].sfc_globalheating[step][threadIdx.x];
		stemp[threadIdx.x] = sfc_stemp1_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
		sfc_stemp0_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = stemp[threadIdx.x];
		stemp[threadIdx.x] = IR_ABS*EPSSIGMA*powf(stemp[threadIdx.x], 4.0f);
		illumination = sfc_illumination_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
		incidentflux = sfc_globalheating_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
		__syncthreads();

		// Determine incident flux by summing over neighbouring facets

		for (i = 0; i != SFCNUMFACE; ++i) {
			incidentflux = incidentflux + (stemp[i] * viewfactors[(i*SFCNUMFACE) + threadIdx.x]);
		}
		__syncthreads();

		// Instantaneous radiative equilibrium

		stemp[threadIdx.x] = powf(((((1.0f - albedo)*illumination) + incidentflux) / (EPSSIGMA)), 0.25f);

		// Write updated surface temperature data

		//	RT[blockIdx.x].sfc_stemp1[step][threadIdx.x]=stemp[threadIdx.x];
		sfc_stemp1_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x] = stemp[threadIdx.x];

	}

	__syncthreads();

};


// Code to test that asteroid surface roughness temperatures have converged


//__global__ void converged(Facet_Properties *FP, Roughness_Temperatures *RT, float TACC) {
__global__ void converged(bool *sfc_converged_D, float *sfc_stemp0_D, float *sfc_stemp1_D, float TACC) {

	__shared__ bool converged;

	int step;
	float stemp0;
	float stemp1;
	float DT;

	if (threadIdx.x == 0) {
		//	converged=FP[blockIdx.x].sfc_converged;
		converged = sfc_converged_D[blockIdx.x];
	}
	__syncthreads();

	if (converged == false) {

		__syncthreads();
		if (threadIdx.x == 0) {
			converged = true;
		}
		__syncthreads();

		for (step = 0; step != NTAU; ++step) {

			//	stemp0=RT[blockIdx.x].sfc_stemp0[step][threadIdx.x];
			//	stemp1=RT[blockIdx.x].sfc_stemp1[step][threadIdx.x];
			stemp0 = sfc_stemp0_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];
			stemp1 = sfc_stemp1_D[(blockIdx.x*NTAU*SFCNUMFACE) + (step*SFCNUMFACE) + threadIdx.x];

			DT = fabsf(stemp0 - stemp1);
			if (DT >= TACC) { converged = false; }

		}

	}
	__syncthreads();

	if (threadIdx.x == 0) {
		//	FP[blockIdx.x].sfc_converged=converged;
		sfc_converged_D[blockIdx.x] = converged;
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

	cudaMemcpyToSymbol(PI, &PI_H, sizeof(float));
	cudaMemcpyToSymbol(rot_period, &rot_period_H, sizeof(float));

	cudaMalloc((void **)&viewfactors_D, SFCNUMFACE*SFCNUMFACE * sizeof(float));
	cudaMemcpy(viewfactors_D, viewfactors_H, SFCNUMFACE*SFCNUMFACE * sizeof(float), cudaMemcpyHostToDevice);

	float host_sfc_vertex[SFCNUMVTX][3];
	int host_sfc_facet[SFCNUMFACE][3];
	float host_sfc_midpos[SFCNUMFACE][3];
	float host_sfc_trivector[SFCNUMFACE][2][3];
	float host_sfc_normcalc[SFCNUMFACE][3][3];
	float host_sfc_normal[SFCNUMFACE][3];

	for (int i = 0; i != SFCNUMVTX; ++i) {
		host_sfc_vertex[i][0] = S_H.SV[i].pos[0];
		host_sfc_vertex[i][1] = S_H.SV[i].pos[1];
		host_sfc_vertex[i][2] = S_H.SV[i].pos[2];
	}

	for (int j = 0; j != SFCNUMFACE; ++j) {
		host_sfc_facet[j][0] = S_H.SF[j].sfcvertices[0];
		host_sfc_facet[j][1] = S_H.SF[j].sfcvertices[1];
		host_sfc_facet[j][2] = S_H.SF[j].sfcvertices[2];
		host_sfc_midpos[j][0] = S_H.SF[j].midpos[0];
		host_sfc_midpos[j][1] = S_H.SF[j].midpos[1];
		host_sfc_midpos[j][2] = S_H.SF[j].midpos[2];
		host_sfc_trivector[j][0][0] = S_H.SF[j].trivector[0][0];
		host_sfc_trivector[j][0][1] = S_H.SF[j].trivector[0][1];
		host_sfc_trivector[j][0][2] = S_H.SF[j].trivector[0][2];
		host_sfc_trivector[j][1][0] = S_H.SF[j].trivector[1][0];
		host_sfc_trivector[j][1][1] = S_H.SF[j].trivector[1][1];
		host_sfc_trivector[j][1][2] = S_H.SF[j].trivector[1][2];
		host_sfc_normcalc[j][0][0] = S_H.SF[j].normcalc[0][0];
		host_sfc_normcalc[j][0][1] = S_H.SF[j].normcalc[0][1];
		host_sfc_normcalc[j][0][2] = S_H.SF[j].normcalc[0][2];
		host_sfc_normcalc[j][1][0] = S_H.SF[j].normcalc[1][0];
		host_sfc_normcalc[j][1][1] = S_H.SF[j].normcalc[1][1];
		host_sfc_normcalc[j][1][2] = S_H.SF[j].normcalc[1][2];
		host_sfc_normcalc[j][2][0] = S_H.SF[j].normcalc[2][0];
		host_sfc_normcalc[j][2][1] = S_H.SF[j].normcalc[2][1];
		host_sfc_normcalc[j][2][2] = S_H.SF[j].normcalc[2][2];
		host_sfc_normal[j][0] = S_H.SF[j].normal[0];
		host_sfc_normal[j][1] = S_H.SF[j].normal[1];
		host_sfc_normal[j][2] = S_H.SF[j].normal[2];
	}

	cudaMemcpyToSymbol(sfc_vertex, host_sfc_vertex, SFCNUMVTX * 3 * sizeof(float));
	cudaMemcpyToSymbol(sfc_facet, host_sfc_facet, SFCNUMFACE * 3 * sizeof(int));
	cudaMemcpyToSymbol(sfc_midpos, host_sfc_midpos, SFCNUMFACE * 3 * sizeof(float));
	cudaMemcpyToSymbol(sfc_trivector, host_sfc_trivector, SFCNUMFACE * 2 * 3 * sizeof(float));
	cudaMemcpyToSymbol(sfc_normcalc, host_sfc_normcalc, SFCNUMFACE * 3 * 3 * sizeof(float));
	cudaMemcpyToSymbol(sfc_normal, host_sfc_normal, SFCNUMFACE * 3 * sizeof(float));

	float host_trivector_multiplier[m*m][2];
	int k, l, p;

	for (k = 1; k != m + 1; ++k) {
		for (l = 1; l != k + 1; ++l) {
			p = ((k - 1)*(k - 1)) + (2 * (l - 1));
			host_trivector_multiplier[p][0] = (float)k - (1.0f / 3.0f);
			host_trivector_multiplier[p][1] = (float)l - (2.0f / 3.0f);
		}
	}

	for (k = 2; k != m + 1; ++k) {
		for (l = 1; l != k; ++l) {
			p = ((k - 1)*(k - 1)) + (2 * (l - 1)) + 1;
			host_trivector_multiplier[p][0] = (float)k - (2.0f / 3.0f);
			host_trivector_multiplier[p][1] = (float)l - (1.0f / 3.0f);
		}
	}

	cudaMemcpyToSymbol(trivector_multiplier, host_trivector_multiplier, m*m * 2 * sizeof(float));

};


// Code to initialise GPU facet properties arrays


extern "C"
void initialise_GPU_facet_properties_arrays() {

	//size_t size;
	//size=ASTNUMFACE*sizeof(Facet_Properties);
	//cudaMalloc((void **) &FP_D, size);

	size_t size;
	size = ASTNUMFACE * sizeof(bool);
	cudaMalloc((void **)&sfc_converged_D, size);
	size = ASTNUMFACE * sizeof(float);
	cudaMalloc((void **)&albedo_D, size);
	cudaMalloc((void **)&thermal_inertia_D, size);
//	cudaMalloc((void **)&rotperiod_D, size);
	cudaMalloc((void **)&numrevs_D, size);

};


// Code to initialise GPU roughness shadow arrays


extern "C"
void initialise_GPU_shadow_arrays() {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Shadow);
	//	cudaMalloc((void **) &RS_D, size);

	size_t size;
	size = ASTNUMFACE*NTAU * sizeof(int);
	cudaMalloc((void **)&shadow_D, size);
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMalloc((void **)&sfc_shadow_D, size);
	size = ASTNUMFACE*NTAU * 3 * sizeof(float);
	cudaMalloc((void **)&sfc_illumvector_D, size);
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMalloc((void **)&sfc_illumangle_D, size);

};


// Code to initialise GPU roughness illumination arrays


extern "C"
void initialise_GPU_illumination_arrays() {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Illumination);
	//	cudaMalloc((void **) &RI_D, size);

	size_t size;
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMalloc((void **)&sfc_illumination_D, size);

};


// Code to initialise GPU roughness globalillumination arrays


extern "C"
void initialise_GPU_globalillumination_arrays() {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Illumination);
	//	cudaMalloc((void **) &RGI_D, size);

	size_t size;
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMalloc((void **)&sfc_globalillumination_D, size);

};


// Code to initialise GPU roughness globalheating arrays


extern "C"
void initialise_GPU_globalheating_arrays() {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Globalheating);
	//	cudaMalloc((void **) &RGH_D, size);

	size_t size;
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMalloc((void **)&sfc_globalheating_D, size);

};


// Code to initialise GPU roughness temperature arrays


extern "C"
void initialise_GPU_temperature_arrays() {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Temperatures);
	//	cudaMalloc((void **) &RT_D, size);

	size_t size;
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMalloc((void **)&sfc_stemp0_D, size);
	cudaMalloc((void **)&sfc_stemp1_D, size);
	size = ASTNUMFACE*(NX + 1)*SFCNUMFACE * sizeof(float);
	cudaMalloc((void **)&sfc_temp_D, size);

};


// Code to free up GPU roughness shadow arrays


extern "C"
void freeup_GPU_shadow_arrays() {

	//	cudaFree(RS_D);

	cudaFree(shadow_D);
	cudaFree(sfc_shadow_D);
	cudaFree(sfc_illumvector_D);
	cudaFree(sfc_illumangle_D);

};


// Code to free up GPU roughness globalillumination arrays


extern "C"
void freeup_GPU_globalillumination_arrays() {

	//	cudaFree(RGI_D);

	cudaFree(sfc_globalillumination_D);

};


// Code to free up GPU roughness globalheating arrays


extern "C"
void freeup_GPU_globalheating_arrays() {

	//	cudaFree(RGH_D);

	cudaFree(sfc_globalheating_D);

};


// Code to free up GPU roughness temperature arrays


extern "C"
void freeup_GPU_temperature_arrays() {

	//	cudaFree(RT_D);

	cudaFree(sfc_stemp0_D);
	cudaFree(sfc_stemp1_D);
	cudaFree(sfc_temp_D);

};


// Code to load GPU facet properties data


//extern "C"
//void set_GPU_facet_properties(Facet_Properties *FP_H) {
extern "C"
void set_GPU_facet_properties(bool *sfc_converged_H, float *albedo_H, float *thermal_inertia_H, float *numrevs_H) {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Facet_Properties);
	//	cudaMemcpy(FP_D,FP_H,size,cudaMemcpyHostToDevice);

	size_t size;
	size = ASTNUMFACE * sizeof(bool);
	cudaMemcpy(sfc_converged_D, sfc_converged_H, size, cudaMemcpyHostToDevice);
	size = ASTNUMFACE * sizeof(float);
	cudaMemcpy(albedo_D, albedo_H, size, cudaMemcpyHostToDevice);
	cudaMemcpy(thermal_inertia_D, thermal_inertia_H, size, cudaMemcpyHostToDevice);
//	cudaMemcpy(rotperiod_D, rotperiod_H, size, cudaMemcpyHostToDevice);
	cudaMemcpy(numrevs_D, numrevs_H, size, cudaMemcpyHostToDevice);

};


// Code to load GPU roughness global-illumination data


//extern "C"
//void load_GPU_globalillumination(Roughness_Globalillumination *RGI_H) {
extern "C"
void load_GPU_globalillumination(float *sfc_globalillumination_H) {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Globalillumination);
	//	cudaMemcpy(RGI_D,RGI_H,size,cudaMemcpyHostToDevice);

	size_t size;
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMemcpy(sfc_globalillumination_D, sfc_globalillumination_H, size, cudaMemcpyHostToDevice);

};


// Code to load GPU roughness global-selfheating data


//extern "C"
//void load_GPU_globalheating(Roughness_Globalheating *RGH_H) {
extern "C"
void load_GPU_globalheating(float *sfc_globalheating_H) {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Globalheating);
	//	cudaMemcpy(RGH_D,RGH_H,size,cudaMemcpyHostToDevice);

	size_t size;
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMemcpy(sfc_globalheating_D, sfc_globalheating_H, size, cudaMemcpyHostToDevice);

};


// Code to obtain GPU facet properties data


//extern "C"
//void obtain_GPU_facet_properties(Facet_Properties *FP_H) {
extern "C"
void obtain_GPU_facet_properties(bool *sfc_converged_H, float *albedo_H, float *thermal_inertia_H, float *numrevs_H) {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Facet_Properties);
	//	cudaMemcpy(FP_H,FP_D,size,cudaMemcpyDeviceToHost);

	size_t size;
	size = ASTNUMFACE * sizeof(bool);
	cudaMemcpy(sfc_converged_H, sfc_converged_D, size, cudaMemcpyDeviceToHost);
	size = ASTNUMFACE * sizeof(float);
	cudaMemcpy(albedo_H, albedo_D, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(thermal_inertia_H, thermal_inertia_D, size, cudaMemcpyDeviceToHost);
	// cudaMemcpy(rotperiod_H, rotperiod_D, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(numrevs_H, numrevs_D, size, cudaMemcpyDeviceToHost);

};


// Code to obtain GPU roughness shadow data


//extern "C" 
//void obtain_GPU_roughness_shadows(Roughness_Shadow *RS_H) {
extern "C"
void obtain_GPU_roughness_shadows(int *shadow_H, float *sfc_shadow_H, float *sfc_illumvector_H, float *sfc_illumangle_H) {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Shadow);
	//	cudaMemcpy(RS_H,RS_D,size,cudaMemcpyDeviceToHost);

	size_t size;
	size = ASTNUMFACE*NTAU * sizeof(int);
	cudaMemcpy(shadow_H, shadow_D, size, cudaMemcpyDeviceToHost);
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMemcpy(sfc_shadow_H, sfc_shadow_D, size, cudaMemcpyDeviceToHost);
	size = ASTNUMFACE*NTAU * 3 * sizeof(float);
	cudaMemcpy(sfc_illumvector_H, sfc_illumvector_D, size, cudaMemcpyDeviceToHost);
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMemcpy(sfc_illumangle_H, sfc_illumangle_D, size, cudaMemcpyDeviceToHost);

};


// Code to obtain GPU roughness temperature data


//extern "C" 
//void obtain_GPU_roughness_temperatures(Roughness_Temperatures *RT_H) {
extern "C"
void obtain_GPU_roughness_temperatures(float *sfc_stemp0_H, float *sfc_stemp1_H, float *sfc_temp_H) {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Temperatures);
	//	cudaMemcpy(RT_H,RT_D,size,cudaMemcpyDeviceToHost);

	size_t size;
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMemcpy(sfc_stemp0_H, sfc_stemp0_D, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(sfc_stemp1_H, sfc_stemp1_D, size, cudaMemcpyDeviceToHost);
	size = ASTNUMFACE*(NX + 1)*SFCNUMFACE * sizeof(float);
	cudaMemcpy(sfc_temp_H, sfc_temp_D, size, cudaMemcpyDeviceToHost);

};


// Code to obtain GPU roughness illumination data


//extern "C" 
//void obtain_GPU_roughness_illumination(Roughness_Illumination *RI_H) {
extern "C"
void obtain_GPU_roughness_illumination(float *sfc_globalillumination_H) {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Illumination);
	//	cudaMemcpy(RI_H,RI_D,size,cudaMemcpyDeviceToHost);

	size_t size;
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMemcpy(sfc_globalillumination_H, sfc_globalillumination_D, size, cudaMemcpyDeviceToHost);

};


// Surface roughness shadowing GPU control code


//extern "C"
//void GPU_sfcpartialshadow(Roughness_Shadow *RS_H) {
extern "C"
void GPU_sfcpartialshadow(int *shadow_H, float *sfc_shadow_H, float *sfc_illumvector_H, float *sfc_illumangle_H) {

	//	size_t size;
	//	size=ASTNUMFACE*sizeof(Roughness_Shadow);
	//	cudaMemcpy(RS_D,RS_H,size,cudaMemcpyHostToDevice);

	//	for (int step=0;step!=NTAU;++step) {
	//		sfcpartialshadow <<<ASTNUMFACE,SFCNUMFACE>>> (RS_D,step);
	//		cudaDeviceSynchronize();
	//	}

	size_t size;
	size = ASTNUMFACE*NTAU * sizeof(int);
	cudaMemcpy(shadow_D, shadow_H, size, cudaMemcpyHostToDevice);
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMemcpy(sfc_shadow_D, sfc_shadow_H, size, cudaMemcpyHostToDevice);
	size = ASTNUMFACE*NTAU * 3 * sizeof(float);
	cudaMemcpy(sfc_illumvector_D, sfc_illumvector_H, size, cudaMemcpyHostToDevice);
	size = ASTNUMFACE*NTAU*SFCNUMFACE * sizeof(float);
	cudaMemcpy(sfc_illumangle_D, sfc_illumangle_H, size, cudaMemcpyHostToDevice);

	for (int step = 0; step != NTAU; ++step) {
		sfcpartialshadow << <ASTNUMFACE, SFCNUMFACE >> > (shadow_D, sfc_shadow_D, sfc_illumvector_D, sfc_illumangle_D, step);
		cudaDeviceSynchronize();
	}

};


// Surface roughness single scattering GPU control code


extern "C"
void GPU_single_scattering(float solar) {

	//	for (int step=0;step!=NTAU;++step) {	
	//		roughness_single_scattering <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RS_D,RI_D,RGI_D,viewfactors_D,solar,step);
	//		cudaDeviceSynchronize();
	//	}

	for (int step = 0; step != NTAU; ++step) {
		roughness_single_scattering << <ASTNUMFACE, SFCNUMFACE >> > (albedo_D, sfc_shadow_D, sfc_illumangle_D, sfc_globalillumination_D, sfc_illumination_D, viewfactors_D, solar, step);
		cudaDeviceSynchronize();
	}

};


// Surface roughness multiple scattering GPU control code


extern "C"
void GPU_multiple_scattering(float solar, float scattering_accuracy) {

	//	for (int step=0;step!=NTAU;++step) {	
	//		roughness_multiple_scattering <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RS_D,RI_D,RGI_D,viewfactors_D,solar,scattering_accuracy,step);
	//		cudaDeviceSynchronize();
	//	}

	for (int step = 0; step != NTAU; ++step) {
		roughness_multiple_scattering << <ASTNUMFACE, SFCNUMFACE >> > (albedo_D, sfc_shadow_D, sfc_illumangle_D, sfc_illumination_D, sfc_globalillumination_D, viewfactors_D, solar, scattering_accuracy, step);
		cudaDeviceSynchronize();
	}

};


// Surface roughness initialise temperatures GPU control code


extern "C"
void GPU_zero_temperatures() {

	//	roughness_zero_temperatures <<<ASTNUMFACE,SFCNUMFACE>>> (RT_D);
	//	cudaDeviceSynchronize();

	roughness_zero_temperatures << <ASTNUMFACE, SFCNUMFACE >> > (sfc_stemp0_D, sfc_stemp1_D, sfc_temp_D);
	cudaDeviceSynchronize();

};


// Surface roughness initialise temperatures GPU control code


extern "C"
void GPU_initialise_temperatures(float EPS, float SIGMA, float IR_ABS) {

	float EPSSIGMA;
	EPSSIGMA = EPS*SIGMA;

	//	roughness_initialise_temperatures <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RI_D,RGH_D,RT_D,viewfactors_D,EPSSIGMA);
	//	cudaDeviceSynchronize();

	roughness_initialise_temperatures << <ASTNUMFACE, SFCNUMFACE >> > (albedo_D, sfc_illumination_D, sfc_globalheating_D, sfc_temp_D, sfc_stemp0_D, sfc_stemp1_D, viewfactors_D, EPSSIGMA, IR_ABS);
	cudaDeviceSynchronize();

};


// Surface roughness finite difference GPU control code


extern "C"
void GPU_inertia_difference(float EPS, float SIGMA, float IR_ABS, float TACC) {

	float EPSSIGMA;
	EPSSIGMA = EPS*SIGMA;

	float TACC2;
	TACC2 = TACC / 10.0f;

	for (int step = 0; step != NTAU; ++step) {
		//	finite_difference <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RT_D);
		finite_difference << <ASTNUMFACE, SFCNUMFACE >> > (sfc_converged_D, numrevs_D, sfc_temp_D);
		cudaDeviceSynchronize();
		//	boundary_condition <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RI_D,RGH_D,RT_D,viewfactors_D,EPSSIGMA,TACC2,step);
		boundary_condition << <ASTNUMFACE, SFCNUMFACE >> > (sfc_converged_D, albedo_D, thermal_inertia_D, numrevs_D, sfc_temp_D, sfc_stemp0_D, sfc_stemp1_D, sfc_illumination_D, sfc_globalheating_D, viewfactors_D, EPSSIGMA, IR_ABS, TACC2, step);
		cudaDeviceSynchronize();
	}

};


// Surface roughness instantaneous equilibrium GPU control code


extern "C"
void GPU_inst_equilibrium(float EPS, float SIGMA, float IR_ABS) {

	float EPSSIGMA;
	EPSSIGMA = EPS*SIGMA;

	//	for (int step=0;step!=NTAU;++step) {
	//		roughness_inst_equilibrium <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RI_D,RGH_D,RT_D,viewfactors_D,EPSSIGMA,step);
	//		cudaDeviceSynchronize();
	//	}

	for (int step = 0; step != NTAU; ++step) {
		roughness_inst_equilibrium << <ASTNUMFACE, SFCNUMFACE >> > (sfc_converged_D, albedo_D, sfc_stemp0_D, sfc_stemp1_D, sfc_illumination_D, sfc_globalheating_D, viewfactors_D, EPSSIGMA, IR_ABS, step);
		cudaDeviceSynchronize();
	}

};


// Temperatures converged GPU control code


extern "C"
void GPU_converged(float TACC) {

	//	converged <<<ASTNUMFACE,SFCNUMFACE>>> (FP_D,RT_D,TACC);
	//	cudaDeviceSynchronize();

	converged << <ASTNUMFACE, SFCNUMFACE >> > (sfc_converged_D, sfc_stemp0_D, sfc_stemp1_D, TACC);
	cudaDeviceSynchronize();

};