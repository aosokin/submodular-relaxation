
#include "mex.h"
#include "objfunc.h"

#include <string.h>

static char *objfunc_name = NULL;

//void objfunc_(int *n, double *x, double *f, double *g)
void OBJFUNC(int *n, double *x, double *f, double *g)
{
  int i;
  mxArray *lhs[2];
  mxArray *rhs[1];
  
  lhs[0] = mxCreateDoubleScalar(0.0);
  lhs[1] = mxCreateDoubleMatrix(*n, 1, mxREAL);
  rhs[0] = mxCreateDoubleMatrix(*n, 1, mxREAL);
  double *arr = mxGetPr(rhs[0]);
  for(i = 0; i < *n; i++)
    arr[i] = x[i];
  
  mexCallMATLAB(2, lhs, 1, rhs, objfunc_name);
  
  *f = (double)(mxGetScalar(lhs[0]));
  arr = mxGetPr(lhs[1]);
  if(g != NULL)
  {
    for(i = 0; i < *n; i++)
      g[i] = arr[i];
  }
}

void set_obj_func(const char *name)
{
  free(objfunc_name);
  objfunc_name = (char*)malloc((strlen(name)+1) * sizeof(char));
  strcpy(objfunc_name, name);
}
