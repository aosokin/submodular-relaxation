#include "mex.h"
#include "objfunc.h"

#include <string.h>

#define DBL_EPSILON (2.2204460492503131e-16)

#ifdef __OS_WIN__
	//Windows specific names
	#define __CURRENT_LMBMU_NAME__  LMBMU
#else
	//Linux specific names
	#define __CURRENT_LMBMU_NAME__  lmbmu_
#endif

extern "C" { 
void  __CURRENT_LMBMU_NAME__(int *n, int *na, int *mcu, int *mc, int *nw,
           double *x, double *f, double *rpar, int *ipar,
           int *iout, float *time, float *rtim, double *w);
}

static void printDoubleVector(double *x, int n)
{
  int i;
  printf("[ ");
  for(i = 0; i < n; i++)
  {
    if(x[i] >= 0.01 || x[i] < DBL_EPSILON)
      printf("%.2f ", x[i]);
    else
      printf("%.2e ", x[i]);
  }
  printf("]");
}

static void printIntVector(int *x, int n)
{
  int i;
  printf("[ ");
  for(i = 0; i < n; i++)
    printf("%d ", x[i]);
  printf("]");
}

void checkScalar(const mxArray *p, const char *name)
{
  char errorMsg[1000];
  int error = 0;
  
  errorMsg[0] = '\0';
  strcat(errorMsg, name);
  
  if(mxIsEmpty(p)){
      return;
  }
  
  if(!mxIsDouble(p) || !(mxGetNumberOfElements(p) == 1))
  {
    strcat(errorMsg, " must be a scalar.");
    error = 1;
  }
  
  if(error)
    mexErrMsgTxt(errorMsg);
}

void checkIntScalar(const mxArray *p, const char *name)
{
  char errorMsg[1000];
  int error = 0;
  
  errorMsg[0] = '\0';
  strcat(errorMsg, name);
  
  if(mxIsEmpty(p)){
      return;
  }
  
  if(!mxIsDouble(p) || !(mxGetNumberOfElements(p) == 1))
  {
    strcat(errorMsg, " must be a scalar.");
    error = 1;
  }
  double d = (double)mxGetScalar(p);
  int i = (int)mxGetScalar(p);
  if(d != i)
  {
    strcat(errorMsg, " must be an integer.");
    error = 1;
  }
  
  if(error)
    mexErrMsgTxt(errorMsg);
}

void checkArray(const mxArray *p, const char *name, int n)
{
  char errorMsg[1000];
  int error = 0;
  
  errorMsg[0] = '\0';
  strcat(errorMsg, name);
  
  if(!mxIsDouble(p))
  {
    strcat(errorMsg, " must be a scalar array.");
    error = 1;
  }
  
  int ne = mxGetNumberOfElements(p);
  if(ne != n)
  {
    sprintf(errorMsg, "The length of %s is %d != %d.", name, ne, n);
    error = 1;
  }
  
  if(error)
    mexErrMsgTxt(errorMsg);
}

void checkArgs(int nrhs, const mxArray *prhs[])
{
  if(nrhs < 3)
    mexErrMsgTxt("lmbm_driver requires at least four arguments.");
  
  if(!mxIsChar(prhs[0]))
    mexErrMsgTxt("The objective function name (argument 1) must be a string.");
  checkIntScalar(prhs[2], "n (argument 3)");
  int n = (int)mxGetScalar(prhs[2]);
  checkArray(prhs[1], "x0 (argument 2)", n);
  checkIntScalar(prhs[3], "print (argument 4)");
  if(nrhs >= 5)
    checkScalar(prhs[4], "maxtime (argument 5)");
  if(nrhs >= 6)
    checkIntScalar(prhs[5], "na (argument 6)");
  if(nrhs >= 7)
    checkIntScalar(prhs[6], "mcu (argument 7)");
  if(nrhs >= 8)
    checkIntScalar(prhs[7], "mc (argument 8)");
  
  if(nrhs >= 9)
    checkArray(prhs[8], "rpar (argument 9)", 8);
  if(nrhs >= 10)
    checkArray(prhs[9], "ipar (argument 10)", 7);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char *func_name;
  int i;
  
  checkArgs(nrhs, prhs);
  
  func_name = mxArrayToString(prhs[0]);
  set_obj_func(func_name);
  double *x = mxGetPr(prhs[1]);
  
  int n = (int)mxGetScalar(prhs[2]);
  int print = (int)mxGetScalar(prhs[3]);
  float time = 300.0f;
  int na = 2;
  int mcu = 7;
  int nw;
  int mc = 7;
  double f;
  int ipar[7];
  double rpar[8];
  
  for(i = 0; i < 7; i++)
    ipar[i] = 0;
  ipar[4] = -1;
  ipar[5] = 0;
  ipar[6] = 0;
  ipar[1] = 5000000;
  ipar[2] = 5000000;
  
  for(i = 0; i < 8; i++)
    rpar[i] = 0.0;
  rpar[3] = 1e-5;
  rpar[5] = 0.5;
  
  if(nrhs >= 10 && !mxIsEmpty(prhs[9]))
  {
    double *ipar_ptr = mxGetPr(prhs[9]);
    for(i = 0; i < 7; i++)
        if ( !mxIsNaN(ipar_ptr[i]) )
            ipar[i] = ipar_ptr[i];
// ipar[4] have to be -1 or MatLab will collapse 
    ipar[4] = -1;
  }
  if(nrhs >= 9 && !mxIsEmpty(prhs[8]))
  {
    double *rpar_ptr = mxGetPr(prhs[8]);
    for(i = 0; i < 8; i++)
        if ( !mxIsNaN(rpar_ptr[i]) )
            rpar[i] = rpar_ptr[i];
  }
  
  if(nrhs >= 8 && !mxIsEmpty(prhs[7]))
    mc = (int)mxGetScalar(prhs[7]);
  if(nrhs >= 7 && !mxIsEmpty(prhs[6]) )
    mcu = (int)mxGetScalar(prhs[6]);
  if(nrhs >= 6 && !mxIsEmpty(prhs[5]))
    na = (int)mxGetScalar(prhs[5]);
  if(nrhs >= 5 && !mxIsEmpty(prhs[4]) )
    time = (int)mxGetScalar(prhs[4]);
  
  nw = 1 + 9*n + 2*n*na + 3*na + 2*n*(mcu+1) + 
       3*(mcu+2)*(mcu+1)/2 + 9*(mcu+1);
  
  int iout[3];
  float rtim[2];
  double *w = (double*) calloc( nw, sizeof(double) );
  
  /* compute the initial objective function value */
  int n_ = n;
  __CURRENT_OBJFUNC_NAME__(&n_, x, &f, NULL);
  
  if(print == 1)
  {
    printf("---------------\n");
    printf("| Parameters: |\n");
    printf("---------------\n");
    printf("%-8s %s\n", "f:", func_name);
    printf("%-8s %d\n", "n:", n);
    printf("%-8s %d\n", "print:", print);
    printf("%-8s %f\n", "maxtime:", time);
    printf("%-8s %d\n", "na:", na);
    printf("%-8s %d\n", "mcu:", mcu);
    printf("%-8s %d\n", "mc:", mc);
    printf("rpar: ");
    printDoubleVector(rpar, 8);
    printf("\n");
    printf("ipar: ");
    printIntVector(ipar, 7);
    printf("\n");
  }
  
  __CURRENT_LMBMU_NAME__(&n, &na, &mcu, &mc, &nw, x, &f, rpar, ipar,
        iout, &time, rtim, w);


  
  if(print == 1)
  {
    printf("-----------\n");
    printf("| Output: |\n");
    printf("-----------\n");
    printf("%-16s %d\n", "Termination:", iout[2]);
    printf("%-16s %d\n", "N. iter.:", iout[0]);
    printf("%-16s %d\n", "N. func. eval.:", iout[1]);
    printf("%-16s %f\n", "Final value:", f);
    printf("%-16s %f\n", "Execution time:", rtim[0]);
  }
  
  if(nlhs >= 1)
  {
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    memcpy(mxGetPr(plhs[0]), x, n*sizeof(double));
  }
  if(nlhs >= 2)
  {
    plhs[1] = mxCreateDoubleScalar(0.0);
    *mxGetPr(plhs[1]) = f;
  }
  if(nlhs >= 3)
  {
    plhs[2] = mxCreateDoubleScalar(0.0);
    *mxGetPr(plhs[2]) = iout[0];
  }
  if(nlhs >= 4)
  {
    plhs[3] = mxCreateDoubleScalar(0.0);
    *mxGetPr(plhs[3]) = iout[1];
  }
  if(nlhs >= 5)
  {
    plhs[4] = mxCreateDoubleScalar(0.0);
    *mxGetPr(plhs[4]) = iout[2];
  }
  if(nlhs >= 6)
  {
    plhs[5] = mxCreateDoubleScalar(0.0);
    *mxGetPr(plhs[5]) = rtim[0];
  }
  
  free( w );
}
