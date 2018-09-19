#include <cufft.h>
#include <stdio.h>
#include <mex.h>


class lpRgpu{
	//global parameters
	int Nt,Nx,Nq,Ntau,ni;
	int Ntheta,Nrho;
	int Ntheta_R2C;
	int add;
	//gpu memoy
	float* dftx;
	float* dftauq;
	float* dfl;
	float2* dflc;	
	float2* dfZ;
	int* dst;int* dstadj;int* dnvals;int* didnvals;
	
	int* didthetatx;int* didrhotx;
	int* didthetatauq;int* didrhotauq;
	float* ddthetatx;float* ddrhotx;
	float* ddthetatauq;float* ddrhotauq;
	float* demul;float* dcosmul;
	int* dreorids;int* dreoridsadj;
	float* dJ;
	//fft handles
	cufftHandle plan_forward;
	cufftHandle plan_inverse;
	cufftHandle plan_f_forward;
	cufftHandle plan_f_inverse;

	cudaError_t err;
	int MBS31,MBS32,MBS33;
public:
	//void callErr(const char* err);
	lpRgpu(int* N,float* fZ,int* st,int* stadj,int* idthetatx,int* idrhotx,int* idthetatauq,int* idrhotauq,float* dthetatx,
		float* drhotx,float* dthetatauq,float* drhotauq,float* emul,float* cosmul,float* J,int* reorids,int* reoridsadj);
	~lpRgpu();
	void fwd(float* out,float* in);
	void adj(float* out,float* in);
	void printCurrentGPUMemory(const char* str);
	void fftlp(float* out,float* in);
	void fftlpadj(float* out,float* in);
	void convtx(float* out,float* in);
	void convtauq(float* out,float* in);
	void convtauqadj(float* out,float* in);
	void convtxadj(float* out,float* in);
	void getSizes(size_t *N);

};
