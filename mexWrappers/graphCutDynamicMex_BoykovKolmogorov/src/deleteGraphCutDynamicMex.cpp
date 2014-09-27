#include "graphCutMemory.h"
#include "graphCutMex.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
	if ( nrhs != 1 ) {
		mexErrMsgIdAndTxt("deleteGraphCutDynamicMex:inputArguments","Wrong number of input arguments, expected 1");
    }
	

	 // get graph handle
	GraphType *g = NULL;
    g = getGraphHandle(prhs[0]);

	//free memory
	delete g;
	g = NULL;
}

