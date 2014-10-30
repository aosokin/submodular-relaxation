#include <cstdlib>
#include <stdio.h>
#include <limits>
#include <vector>
using std::vector;

#include "mex.h"

#include "robust/expand.h"

#define MATLAB_ASSERT(expr,msg) if (!(expr)) {mexErrMsgIdAndTxt( "alphaExpRobustMex:error", msg);}

#if !defined(MX_API_VER) || MX_API_VER < 0x07030000
typedef int mwSize;
typedef int mwIndex;
#endif

double round(double a);
int isInteger(double a);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MATLAB_ASSERT(nrhs == 2 || nrhs == 4, "Wrong number of input argumets, expected 2 or 4" ); \
	
	MATLAB_ASSERT(nlhs <= 4, "Too many output arguments, expected 1 - 4");
	
	//Fix input parameter order:
	const mxArray *uInPtr = (nrhs > 0) ? prhs[0] : NULL; //unary
	const mxArray *pInPtr = (nrhs > 1) ? prhs[1] : NULL; //pairwise
	const mxArray *hoInPtr = (nrhs > 2) ? prhs[2] : NULL; // high-order indices
	const mxArray *hopInPtr = (nrhs > 3) ? prhs[3] : NULL; // high-order parameters
	
	//Fix output parameter order:
	mxArray **eOutPtr = (nlhs > 1) ? &plhs[1] : NULL; //energy
	mxArray **sOutPtr = (nlhs > 0) ? &plhs[0] : NULL; //solution
	mxArray **timePlotOutPtr = (nlhs > 3) ? &plhs[3] : NULL; //time plot 
	mxArray **energyPlotOutPtr = (nlhs > 2) ? &plhs[2] : NULL; //energy plot 

	// get unary potentials
	MATLAB_ASSERT(mxGetNumberOfDimensions(uInPtr) == 2, "Unary term array is not 2-dimensional");
	MATLAB_ASSERT(mxGetPi(uInPtr) == NULL, "Unary potentials should not be complex");
	
	mwSize numNodes = mxGetN(uInPtr);
	mwSize numLabels = mxGetM(uInPtr);

	MATLAB_ASSERT(numNodes >= 1, "The number of nodes is not positive");
	MATLAB_ASSERT(numLabels >= 1, "The number of labels is not positive");
	MATLAB_ASSERT(mxGetClassID(uInPtr) == mxDOUBLE_CLASS, "Expected mxDOUBLE_CLASS for input unary term argument");
	double* termW = (double*)mxGetData(uInPtr);

	//get pairwise potentials
	MATLAB_ASSERT(mxIsSparse(pInPtr), "Expected sparse array for neighbours");
	MATLAB_ASSERT(mxGetN(pInPtr) == numNodes && mxGetM(pInPtr) == numNodes,
	              "Neighbours array must be NumNodes x NumNodes in size");
	MATLAB_ASSERT(mxGetClassID(pInPtr) == mxDOUBLE_CLASS, "Expected mxDOUBLE_CLASS for neighbours array");
	MATLAB_ASSERT(mxGetPi(pInPtr) == NULL, "Pairwise potentials should not be complex");

	mwIndex colNum = (mwIndex)mxGetN(pInPtr);
	const mwIndex* ir = mxGetIr(pInPtr);
	const mwIndex* jc = mxGetJc(pInPtr);
	double*        pr = mxGetPr(pInPtr);

	// get high-order indices
	int numHO = 0;
	vector< vector< int > > highOrderInd;
	vector< vector< double > > highOrderParam;
	if (hoInPtr != NULL) {
		MATLAB_ASSERT(mxIsCell( hoInPtr ), "HO parameter shoulf be cell array");
		numHO = mxGetNumberOfElements( hoInPtr );

		highOrderInd.resize( numHO );
		highOrderParam.resize( numHO );

		for(int iPot = 0; iPot < numHO; ++iPot) {
			mxArray *cellPtr = mxGetCell( hoInPtr, iPot );

			MATLAB_ASSERT(mxGetClassID(cellPtr) == mxDOUBLE_CLASS, "Expected mxDOUBLE_CLASS for elements of HO parameter");
			MATLAB_ASSERT(mxGetNumberOfDimensions(cellPtr) == 2, "Elements of HO parameter should be 2-dimensional");
			MATLAB_ASSERT(mxGetPi(cellPtr) == NULL, "Elements of HO parameter should not be complex");
			
			mwSize curLength = mxGetM(cellPtr);
			MATLAB_ASSERT(mxGetN(cellPtr) == 1, "Elements of HO parameter should not be column vectors");

			double* curVector = (double*)mxGetData(cellPtr);

			highOrderInd[iPot].resize( curLength );

			for(int iElem = 0; iElem < curLength; ++iElem) {
				MATLAB_ASSERT( isInteger( curVector[ iElem ] ) && curVector[ iElem ] > 0.5 && curVector[ iElem ] < numNodes + 0.5, "Elements of HO elements should be integer numbers 1,...,numNodes"  )

				highOrderInd[iPot][iElem] = round( curVector[ iElem ] ) - 1;
			}
		}
	}

	// get HOP
	if ( hopInPtr != NULL ) {
		MATLAB_ASSERT(mxGetClassID( hopInPtr ) == mxDOUBLE_CLASS, "Expected mxDOUBLE_CLASS for HOP");

        MATLAB_ASSERT(mxGetNumberOfDimensions( hopInPtr ) == 2, "HOP should be 2-dimensional");
		MATLAB_ASSERT(mxGetPi( hopInPtr ) == NULL, "HOP should not be complex");
			
		MATLAB_ASSERT(mxGetM(hopInPtr) == numHO && mxGetN(hopInPtr) == 2, "HOP parameter should be of size numHO x 2");

		double* curVector = (double*)mxGetData(hopInPtr);

		for( int iPot = 0; iPot < numHO; ++iPot) {
			highOrderParam[ iPot ].resize( 2 );

			highOrderParam[ iPot ][0] = curVector[ iPot ];
			highOrderParam[ iPot ][1] = curVector[ iPot + numHO ];

			MATLAB_ASSERT( highOrderParam[ iPot ][0] >= 0, "gamma_max of HOP should be >= 0");
			MATLAB_ASSERT( highOrderParam[ iPot ][1] >= 0, "Q of HOP should be >= 0");
		}

	}

	// Potts MRF
	//check pairwise terms
	mwSize numEdges = 0;
	for (mwIndex c = 0; c < colNum; ++c) {
			mwIndex rowStart = jc[c]; 
			mwIndex rowEnd   = jc[c+1]; 
			for (mwIndex ri = rowStart; ri < rowEnd; ++ri)  {
				mwIndex r = ir[ri];
	
				double dw = pr[ri];
				if( r < c) numEdges++;
				MATLAB_ASSERT( dw >=0, "Some Potts edge have negative coefficient!");
			}
		}

	// create the energy
	Energy<double> *energy = new Energy<double>(numLabels, numNodes, numEdges, numHO);

	//initialize random unary costs (variables x labels)
	for(int iNode = 0; iNode < numNodes; ++iNode)
		for(int iLabel = 0; iLabel < numLabels; ++iLabel){
				energy->unaryCost[iNode * numLabels + iLabel] = termW[iLabel + iNode * numLabels];
			}

	//add edges
	int iEdge = 0;
	for (mwIndex c = 0; c < colNum; ++c) {
			mwIndex rowStart = jc[c]; 
			mwIndex rowEnd   = jc[c+1]; 
			for (mwIndex ri = rowStart; ri < rowEnd; ++ri)  {
				mwIndex r = ir[ri];
	
				double dw = pr[ri];
				if( r < c) {

					energy->pairIndex[iEdge * 2] = r;
					energy->pairIndex[iEdge * 2 + 1] = c;
					energy->pairCost[iEdge] = dw;

					++iEdge;
				}
			}
		}

	
	//initialize number of elements in each segment
	for( int iHO = 0; iHO < numHO; ++iHO){
		energy->higherElements[ iHO ] = highOrderInd[ iHO ].size();
	}
	
	//allocate energy for higher order indexes
	energy->AllocateHigherIndexes();

	//initialize indexes of pixels in each segment
	for( int iHO = 0; iHO < numHO; ++iHO)
		for (int iNode = 0; iNode < highOrderInd[ iHO ].size(); ++iNode) {
			energy->higherIndex[ iHO ][ iNode ] = highOrderInd[ iHO ][ iNode ];
		}
	
	//initialize truncation ratio Q, gamma_k and gamma_max for each clique
	for(int iHO = 0; iHO < numHO; ++iHO)
	{
		//truncation ratio 30%
		energy->higherTruncation[ iHO ] = highOrderParam[ iHO ][1];

		//gamma_k
		for(int iLabel = 0; iLabel < numLabels; ++iLabel) 
			energy->higherCost[iHO * (numLabels + 1) + iLabel] = 0;

		//gamma_max
		energy->higherCost[iHO * (numLabels + 1) + numLabels] = highOrderParam[ iHO ][0];
	}

	//initialize alpha expansion - max 50 iterations
	AExpand *expand = new AExpand(energy, 50);

	//initialize solution
	int *solution = new int[numNodes];
	memset(solution, 0, numNodes * sizeof(int));

	//solve CRF
	double energyValue = expand->minimize(solution);

	//output the best energy value
	if(eOutPtr != NULL)	{
		*eOutPtr = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		*(double*)mxGetData(*eOutPtr) = (double)energyValue;
	}

	//output the best solution
	if(sOutPtr != NULL)	{
		*sOutPtr = mxCreateNumericMatrix(numNodes, 1, mxDOUBLE_CLASS, mxREAL);
		double* segment = (double*)mxGetData(*sOutPtr);
		for(int i = 0; i < numNodes; ++i)
			segment[i] = solution[i] + 1;
	}

	//output time plot
	if(timePlotOutPtr != NULL)	{
		*timePlotOutPtr = mxCreateNumericMatrix(expand->timePlot.size(), 1, mxDOUBLE_CLASS, mxREAL);
		double* time = (double*)mxGetData(*timePlotOutPtr);
		for(int i = 0; i < expand->timePlot.size(); ++i)
			time[i] = expand->timePlot[i];
	}

	//output energy plot
	if(energyPlotOutPtr != NULL)	{
		*energyPlotOutPtr = mxCreateNumericMatrix(expand->energyPlot.size(), 1, mxDOUBLE_CLASS, mxREAL);
		double* curE = (double*)mxGetData(*energyPlotOutPtr);
		for(int i = 0; i < expand->energyPlot.size(); ++i)
			curE[i] = expand->energyPlot[i];
	}

    delete[] solution;
	delete expand;
	delete energy;
}




inline double round(double a)
{
	return (int)floor(a + 0.5);
}

inline int isInteger(double a)
{
	return (a - round(a) < 1e-6);
}

