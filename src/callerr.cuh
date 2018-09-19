#include<mex.h>
void callErr(const char* str)
{
	mexErrMsgTxt(str);
	mexErrMsgTxt("Reset gpu");
	cudaDeviceReset();
}
