#include "graphCutMemory.h"
#include "graphCutMex.h"
#include "mex.h"

#include <limits>
#include <cmath>




void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
	if ( nrhs != 2 ) {
		mexErrMsgIdAndTxt("updateUnaryGraphCutDynamicMex:parameters","Wrong number of input parameter, expected 2");
    }
	
	// set up pointers for input/ output parameters
	const mxArray* graphHandleInPtr = prhs[0]; //graphHandle
	const mxArray* updateInPtr = prhs[1]; // the update array 
	mxArray **energyOutPtr = (nlhs > 0) ? &plhs[0] : NULL; //energy
	mxArray **labelsOutPtr = (nlhs > 1) ? &plhs[1] : NULL; //labeling
	
	 // get graph handle
	GraphType *g = NULL;
    g = getGraphHandle( graphHandleInPtr );

	// get the cnahges
	if (mxGetNumberOfDimensions( updateInPtr ) != 2)	{
			mexErrMsgIdAndTxt("updateUnaryGraphCutDynamicMex:updateUnaryWrongDimension","updateUnary is not 2-dimensional");
	}
	int numChanges = mxGetM( updateInPtr );
	if (mxGetN( updateInPtr ) != 3){
		mexErrMsgIdAndTxt("updateUnaryGraphCutDynamicMex:updateUnaryWrongDimension","updateUnary is not of size #changes x 3");
	}
	if (mxGetClassID( updateInPtr ) != MATLAB_ENERGYTERM_TYPE ) {
		mexErrMsgIdAndTxt("updateUnaryGraphCutDynamicMex:updateUnaryWrongType", "updateUnary is of wrong type");
	}
	EnergyTermType* changes = (EnergyTermType*)mxGetData( updateInPtr );

	int numNodes = g -> get_node_num();
	
	//start editing graph
	for(int i = 0; i < numChanges; ++i)
		if(!isInteger(changes[i]) || changes[i] < 1 || changes[i] > numNodes){
			mexErrMsgIdAndTxt("updateUnaryGraphCutDynamicMex:updateUnaryWrongNodeId", "updateUnary has one nodeId incorrect");
		}
		else
		{
			GraphType::node_id j = (GraphType::node_id)round(changes[i] - 1);
			g -> add_tweights(j, changes[i + numChanges], changes[i + 2 * numChanges]);
			g -> mark_node(j);
		}

	if (energyOutPtr == NULL) return;
	
	*energyOutPtr = mxCreateNumericMatrix(1, 1, MATLAB_ENERGY_TYPE, mxREAL);
	*(EnergyType*)mxGetData(*energyOutPtr) = (EnergyType)(g -> maxflow(true));

	if( labelsOutPtr != NULL )	{
		*labelsOutPtr = mxCreateNumericMatrix(numNodes, 1, MATLAB_LABEL_TYPE, mxREAL);
		LabelType* segment = (LabelType*)mxGetData(*labelsOutPtr);
		for(int i = 0; i < numNodes; i++)
			segment[i] = g -> what_segment(i);
	}
}


