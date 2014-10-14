
#include "mex.h"

#include <vector>
using std::vector;

#include <cmath>

#define MATLAB_ASSERT(expr,errorId,msg) if (!(expr)) {mexErrMsgIdAndTxt(errorId,msg);}

double round(double a);
int isInteger(double a);

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
	MATLAB_ASSERT( nrhs >= 2  && nrhs <= 4, "icmPottsMex:inputParameters", "Wrong number of input input arguments, expected 2 - 4");
	MATLAB_ASSERT( nlhs <= 2, "icmPottsMex:outputParameters", "Too many output arguments, expected 0 - 2");

	// set up pointers for input/ output parameters
	const mxArray* unaryInPtr = prhs[0]; //unary terms
	const mxArray* pairwiseInPtr = prhs[1]; //pairwise terms
	const mxArray* initLabelsInPtr = (nrhs > 2) ? prhs[2] : NULL; // the initial labeling
	const mxArray* maxNumIterInPtr = (nrhs > 3) ? prhs[3] : NULL; // the maximum number of iterations

	mxArray **energyOutPtr = (nlhs > 0) ? &plhs[0] : NULL; //energy
	mxArray **labelsOutPtr = (nlhs > 1) ? &plhs[1] : NULL; //labeling
	
	int numNodes = 0;
	int numLabels = 0;
	int numEdges = 0;
	int maxIter = 10;
	
	// get unary potentials
	MATLAB_ASSERT(  mxGetClassID( unaryInPtr ) == mxDOUBLE_CLASS, "icmPottsMex:unaryPotentialsWrongType", "unaryTerms is of wrong type, expected double");
	MATLAB_ASSERT(  mxGetNumberOfDimensions( unaryInPtr ) == 2, "icmPottsMex:unaryPotentialsWrongDimensionality", "unaryTerms is not 2-dimensional");
	MATLAB_ASSERT(  mxGetPi(unaryInPtr) == NULL, "icmPottsMex:unaryPotentialsComplex",  "Unary potentials should not be complex");
	numNodes = mxGetN(unaryInPtr);
	numLabels = mxGetM(unaryInPtr);
	MATLAB_ASSERT(numNodes >= 1, "icmPottsMex:unaryPotentialsWrongNumNodes", "The number of nodes is not positive");
	MATLAB_ASSERT(numLabels >= 1, "icmPottsMex:unaryPotentialsWrongNumLabels", "The number of labels is not positive");

	double* dataCost = (double*)mxGetData(unaryInPtr);
	

	//get pairwise potentials
	MATLAB_ASSERT(mxGetClassID(pairwiseInPtr) == mxDOUBLE_CLASS, "icmPottsMex:pairwisePotentialsWrongType", "Expected mxDOUBLE_CLASS for neighbours array");
	MATLAB_ASSERT(mxIsSparse(pairwiseInPtr), "icmPottsMex:pairwisePotentialsNotSparse", "Expected sparse array for neighbours");
	MATLAB_ASSERT(mxGetN(pairwiseInPtr) == numNodes && mxGetM(pairwiseInPtr) == numNodes, "icmPottsMex:pairwisePotentialsWrongSize", "Neighbours array must be NumNodes x NumNodes in size");
	MATLAB_ASSERT(mxGetPi(pairwiseInPtr) == NULL, "icmPottsMex:pairwisePotentialsComplex",  "Pairwise potentials should not be complex");

	mwIndex colNum = (mwIndex)mxGetN(pairwiseInPtr);
	const mwIndex* ir = mxGetIr(pairwiseInPtr);
	const mwIndex* jc = mxGetJc(pairwiseInPtr);
	double*        pr = mxGetPr(pairwiseInPtr);

	//check pairwise terms
	for (mwIndex c = 0; c < colNum; ++c) {
		mwIndex rowStart = jc[c]; 
		mwIndex rowEnd   = jc[c+1]; 
		for (mwIndex ri = rowStart; ri < rowEnd; ++ri)  {
			mwIndex r = ir[ri];

			double dw = pr[ri];
			if( r < c) {
				++numEdges;
				// MATLAB_ASSERT(dw >= 0, "icmPottsMex:pairwisePotentialsNegativeValues", "Potts potentials should be positive");
			}
		}
	}

	// get the initial labeling
	vector<int> curLabeling(numNodes);
	if ( initLabelsInPtr != NULL ) {
		if ( mxIsEmpty(initLabelsInPtr) ) {
			initLabelsInPtr = NULL;
		}
		else {
			MATLAB_ASSERT(  mxGetClassID( initLabelsInPtr ) == mxDOUBLE_CLASS, "icmPottsMex:initLabelsWrongType", "initLabels is of wrong type, expected double");
			MATLAB_ASSERT(  mxGetNumberOfDimensions( initLabelsInPtr ) == 2, "icmPottsMex:initLabelsWrongDimensionality", "initLabels is not 2-dimensional");
			MATLAB_ASSERT(  mxGetPi(initLabelsInPtr) == NULL, "icmPottsMex:initLabelsComplex",  "initLabels should not be complex");
			MATLAB_ASSERT(  mxGetM(initLabelsInPtr) == numNodes && mxGetN(initLabelsInPtr) == 1, "icmPottsMex:initLabelsComplex",  "initLabels is of wrong size, expected numNodes x 1");
			
			double* labelsPtr = (double*)mxGetData(initLabelsInPtr);

			for( int iNode = 0; iNode < numNodes; ++iNode ) {
				double tmp = labelsPtr[iNode];
				MATLAB_ASSERT(  tmp >= 1 && tmp <= numLabels && isInteger(tmp), "icmPottsMex:initLabelsWrongLabel",  "initLabels contains wrong label, should be in 1,...,numLabels");
				curLabeling[iNode] = round(tmp) - 1;
			}
		}
	} 
	
	if( initLabelsInPtr == NULL ) {
		for( int iNode = 0; iNode < numNodes; ++iNode ) {
			curLabeling[iNode] = rand() % numLabels;
		}
	}

	// get the maximum number of iterations
	if ( maxNumIterInPtr != NULL) {
		if ( !mxIsEmpty(maxNumIterInPtr) ) {
			MATLAB_ASSERT(  mxGetClassID( maxNumIterInPtr ) == mxDOUBLE_CLASS, "icmPottsMex:maxNumIterWrongType", "maxNumIter is of wrong type, expected double");
			MATLAB_ASSERT(  mxGetPi(maxNumIterInPtr) == NULL, "icmPottsMex:maxNumIterComplex",  "maxNumIter should not be complex");
			MATLAB_ASSERT(  mxGetNumberOfElements( maxNumIterInPtr ) == 1, "icmPottsMex:maxNumIterWrongSize", "maxNumIter is not scalar");

			maxIter = (int)(*((double*) mxGetData( maxNumIterInPtr )));
			MATLAB_ASSERT( maxIter >= 0 , "icmPottsMex:maxNumIterNegativeValue", "maxNumIter is negative");
		}
	}

		
	// start computing

	double energy = 0.0;	
	vector< vector<int> > graphIds(numNodes);
	vector< vector<double> > graphWeights(numNodes);

	//construct and compute pairwise terms
	for (mwIndex c = 0; c < colNum; ++c) {
		mwIndex rowStart = jc[c]; 
		mwIndex rowEnd   = jc[c+1]; 
		for (mwIndex ri = rowStart; ri < rowEnd; ++ri)  {
			mwIndex r = ir[ri];

			double dw = pr[ri];
			if( r < c) {
				graphIds[r].push_back(c);
				graphWeights[r].push_back(dw);

				graphIds[c].push_back(r);
				graphWeights[c].push_back(dw);

				if  ( curLabeling[r] != curLabeling[c]) {
					energy += dw;
				}
			}
		}
	}
	// add unary terms
	for( int iNode = 0; iNode < numNodes; ++iNode ) {
		energy += dataCost[ curLabeling[iNode] + iNode * numLabels ];
	}

	
	for ( int iIter = 0; iIter < maxIter; ++iIter ) {
		bool changed = false;

		for (int iNode = 0; iNode < numNodes; ++iNode ) {
			vector<double> neighWeights(numLabels, 0.0);
			for ( int iEdge = 0; iEdge < graphIds[iNode].size(); ++iEdge) {
				neighWeights[ curLabeling[ graphIds[iNode][iEdge] ] ] += graphWeights[iNode][iEdge];
			}

			// try different labels
			for( int iLabel = 0; iLabel < numLabels; ++iLabel) {
				double diff = dataCost[ iLabel  + iNode * numLabels ] - dataCost[ curLabeling[iNode] + iNode * numLabels ];

				diff += -neighWeights[ iLabel ] + neighWeights[ curLabeling[iNode] ];

				if ( diff < 0 ) {
					curLabeling[ iNode ] = iLabel;
					energy += diff;
					changed = true;
				}

			}
		}
		if (!changed)
			break;		
	}


	
	
	//output minimum value
	if (energyOutPtr != NULL){
		*energyOutPtr = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		*(double*)mxGetData( *energyOutPtr ) = energy;
	}

	
	//output minimum cut
	if ( labelsOutPtr != NULL ){

		*labelsOutPtr = mxCreateNumericMatrix(numNodes, 1, mxDOUBLE_CLASS, mxREAL);
		double* segment = (double*)mxGetData( *labelsOutPtr );
		for(int iNode = 0; iNode < numNodes; ++iNode)
			segment[ iNode ] = curLabeling[ iNode ] + 1;
	}
}

double round(double a)
{
	return (int)floor(a + 0.5);
}

int isInteger(double a)
{
	return (a - round(a) < 1e-6);
}
