#include "graphCutMemory.h"
#include "graphCutMex.h"
#include "mex.h"

#include <limits>
#include <cmath>


void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
	if ( nrhs != 2 ) {
		mexErrMsgIdAndTxt("graphCutDynamicMex:parameters", "Wrong number of input input arguments, expected 2");
    }
	if (nlhs > 3) {
		mexErrMsgIdAndTxt("graphCutDynamicMex:parameters", "Too many output arguments, expected 1 - 3");
	}

	// set up pointers for input/ output parameters
	const mxArray* unaryInPtr = prhs[0]; //unary terms
	const mxArray* pairwiseInPtr = prhs[1]; //pairwise terms
	mxArray **energyOutPtr = (nlhs > 0) ? &plhs[0] : NULL; //energy
	mxArray **labelsOutPtr = (nlhs > 1) ? &plhs[1] : NULL; //labeling
	mxArray **graphHandleOutPtr = (nlhs > 2) ? &plhs[2] : NULL; //graphHandle

	int numNodes = 0;
	int numEdges = 0;
	EnergyTermType* termW = NULL;
	EnergyTermType* edges = NULL;
	
	// get unary potentials
	if ( mxGetClassID( unaryInPtr ) != MATLAB_ENERGYTERM_TYPE ) {
		mexErrMsgIdAndTxt("graphCutDynamicMex:unaryPotentials", "unaryTerms is of wrong type, expected double");
	}
	if ( mxGetNumberOfDimensions( unaryInPtr ) != 2 )	{
		mexErrMsgIdAndTxt("graphCutDynamicMex:unaryPotentials","unaryTerms is not 2-dimensional");
	}
	numNodes = mxGetM(unaryInPtr);
	if ( mxGetN( unaryInPtr ) != 2 ) {
		mexErrMsgIdAndTxt("graphCutDynamicMex:unaryPotentials","unaryTerms is of wrong size, expected #node x 2");
	}
	
	termW = (EnergyTermType*)mxGetData(unaryInPtr);


	// get pairwise potentials
	if (mxGetClassID(pairwiseInPtr) != MATLAB_ENERGYTERM_TYPE ) {
		mexErrMsgIdAndTxt("graphCutDynamicMex:pairwisePotentials", "pairwiseTerms is of wrong type, expected double");
	}
	if (mxGetNumberOfDimensions(pairwiseInPtr) != 2)	{
		mexErrMsgIdAndTxt("graphCutDynamicMex:pairwisePotentials","pairwiseTerms is not 2-dimensional");
	}
	numEdges = mxGetM(pairwiseInPtr);
	if (mxGetN(pairwiseInPtr) != 4){
		mexErrMsgIdAndTxt("graphCutDynamicMex:pairwisePotentials","pairwiseTerms is of wrong size, expected #edges x 4");
	}
	edges = (EnergyTermType*)mxGetData(pairwiseInPtr);
	
	
	// start computing
	
	//prepare graph
	GraphType *g = new GraphType( numNodes, numEdges); 
	
	for(int i = 0; i < numNodes; ++i)
	{
		g -> add_node(); 
		g -> add_tweights( i, termW[i], termW[numNodes + i]); 
	}
	
	for(int i = 0; i < numEdges; ++i)
		if(edges[i] < 1 || edges[i] > numNodes || edges[numEdges + i] < 1 || edges[numEdges + i] > numNodes || edges[i] == edges[numEdges + i] || !isInteger(edges[i]) || !isInteger(edges[numEdges + i])){
			mexErrMsgIdAndTxt("graphCutDynamicMex:pairwisePotentialsWrongIndices", "Some edge has invalid vertex numbers");
		}
		else
			if(edges[2 * numEdges + i] + edges[3 * numEdges + i] < 0){
				mexErrMsgIdAndTxt("graphCutDynamicMex:pairwisePotentialsNonsubmodular", "Some edge is non-submodular");
			}
			else
			{
				if (edges[2 * numEdges + i] >= 0 && edges[3 * numEdges + i] >= 0)
					g -> add_edge((GraphType::node_id)round(edges[i] - 1), (GraphType::node_id)round(edges[numEdges + i] - 1), edges[2 * numEdges + i], edges[3 * numEdges + i]);
				else
					if (edges[2 * numEdges + i] <= 0 && edges[3 * numEdges + i] >= 0)
					{
						g -> add_edge((GraphType::node_id)round(edges[i] - 1), (GraphType::node_id)round(edges[numEdges + i] - 1), 0, edges[3 * numEdges + i] + edges[2 * numEdges + i]);
						g -> add_tweights((GraphType::node_id)round(edges[i] - 1), 0, edges[2 * numEdges + i]); 
						g -> add_tweights((GraphType::node_id)round(edges[numEdges + i] - 1), 0 , -edges[2 * numEdges + i]); 
					}
					else
						if (edges[2 * numEdges + i] >= 0 && edges[3 * numEdges + i] <= 0)
						{
							g -> add_edge((GraphType::node_id)round(edges[i] - 1), (GraphType::node_id)round(edges[numEdges + i] - 1), edges[3 * numEdges + i] + edges[2 * numEdges + i], 0);
							g -> add_tweights((GraphType::node_id)round(edges[i] - 1),0 , -edges[3 * numEdges + i]); 
							g -> add_tweights((GraphType::node_id)round(edges[numEdges + i] - 1), 0, edges[3 * numEdges + i]); 
						}
						else
							mexErrMsgIdAndTxt("computeMarginalsMex:pairwisePotentialsStrangeError", "Something strange with an edge: you should never see this message");
			}

	//compute flow
	EnergyType flow = g -> maxflow();

	//output minimum value
	if (energyOutPtr != NULL){
		*energyOutPtr = mxCreateNumericMatrix(1, 1, MATLAB_ENERGY_TYPE, mxREAL);
		*(EnergyType*)mxGetData( *energyOutPtr ) = (EnergyType)flow;
	}

	
	//output minimum cut
	if ( labelsOutPtr != NULL ){

		*labelsOutPtr = mxCreateNumericMatrix(numNodes, 1, MATLAB_LABEL_TYPE, mxREAL);
		LabelType* segment = (LabelType*)mxGetData( *labelsOutPtr );
		for(int i = 0; i < numNodes; i++)
			segment[i] = g -> what_segment(i);
	}

	if ( graphHandleOutPtr != NULL ) {
			//create a container for the pointer 
			*graphHandleOutPtr = mxCreateNumericMatrix(1, 1, MATLAB_POINTER_TYPE, mxREAL);
    
			*(GraphHandle*)mxGetData( *graphHandleOutPtr ) = (GraphHandle)g;
	}
	else
		delete g;
}

