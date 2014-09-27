#include "mex.h"
#define MATLAB_ASSERT(expr,errorId,msg) if (!(expr)) {mexErrMsgIdAndTxt(errorId,msg);}

#include <vector>
using std::vector;

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
	MATLAB_ASSERT( nrhs == 2, "viterbiPottsMex:inputParameters", "Wrong number of input input arguments, expected 2");
	MATLAB_ASSERT( nlhs <= 2, "viterbiPottsMex:outputParameters", "Too many output arguments, expected 0 - 2");

	// set up pointers for input/ output parameters
	const mxArray* unaryInPtr = prhs[0]; //unary terms
	const mxArray* pairwiseInPtr = prhs[1]; //pairwise terms
	
	mxArray **energyOutPtr = (nlhs > 0) ? &plhs[0] : NULL; //energy
	mxArray **labelsOutPtr = (nlhs > 1) ? &plhs[1] : NULL; //labeling
	
	int numNodes = 0;
	int numLabels = 0;
	
	// get unary potentials
	MATLAB_ASSERT(  mxGetClassID( unaryInPtr ) == mxDOUBLE_CLASS, "viterbiPottsMex:unaryPotentialsWrongType", "unaryTerms is of wrong type, expected double");
	MATLAB_ASSERT(  mxGetNumberOfDimensions( unaryInPtr ) == 2, "viterbiPottsMex:unaryPotentialsWrongDimensionality", "unaryTerms is not 2-dimensional");
	MATLAB_ASSERT(  mxGetPi(unaryInPtr) == NULL, "viterbiPottsMex:unaryPotentialsComplex",  "Unary potentials should not be complex");
	numNodes = mxGetM(unaryInPtr);
	numLabels = mxGetN(unaryInPtr);
	MATLAB_ASSERT(numNodes >= 1, "viterbiPottsMex:unaryPotentialsWrongNumNodes", "The number of nodes is not positive");
	MATLAB_ASSERT(numLabels >= 1, "viterbiPottsMex:unaryPotentialsWrongNumLabels", "The number of labels is not positive");

	double* dataCost = (double*)mxGetData(unaryInPtr);
	

	//get pairwise potentials
	MATLAB_ASSERT(mxGetClassID(pairwiseInPtr) == mxDOUBLE_CLASS, "viterbiPottsMex:pairwisePotentialsWrongType", "Expected mxDOUBLE_CLASS for neighbours array");
	MATLAB_ASSERT((mxGetM(pairwiseInPtr) == numNodes - 1) && mxGetN(pairwiseInPtr) == 1, "viterbiPottsMex:pairwisePotentialsWrongSize", "Pairwise array must be (numNodes-1) x 1 in size");
	MATLAB_ASSERT(mxGetPi(pairwiseInPtr) == NULL, "viterbiPottsMex:pairwisePotentialsComplex",  "Pairwise potentials should not be complex");

	double* pairwiseCost = (double*)mxGetData(pairwiseInPtr);
	
	//for(int iNode = 0; iNode < numNodes - 1; ++iNode)
	//	MATLAB_ASSERT( pairwiseCost[ iNode ] >= NULL, "viterbiPottsMex:pairwisePotentialsNegative",  "Pairwise potentials should not be negative");
    
	// start computing
	vector<int> prevPosition( (numNodes - 1) * numLabels, 0 ); //caution: prevPosition and curCost are stored label-first
	vector<double> curCost( numNodes * numLabels, 0 );

	// put the first line
	for(int iLabel = 0; iLabel < numLabels; ++iLabel)
		curCost[ iLabel + 0 * numLabels ] = dataCost[ 0 + numNodes * iLabel ];

    if (numLabels > 1) {
        for(int iNode = 1; iNode < numNodes; ++iNode) {
    		// find the position of the two smallest elements
            bool firstMinFlag = (curCost[ 0 + numLabels * (iNode - 1) ] < curCost[ 1 + numLabels * (iNode - 1) ]);
    		int minPoint = firstMinFlag ? 0 : 1;
            int secondMinPoint = firstMinFlag ? 1 : 0;
            
    		for(int iLabel = 2; iLabel < numLabels; ++iLabel) {
        		if( curCost[ iLabel + numLabels * (iNode - 1) ] < curCost[ minPoint + numLabels * (iNode - 1) ] ) {
                    secondMinPoint = minPoint;
            		minPoint = iLabel;
                } else if ( curCost[ iLabel + numLabels * (iNode - 1) ] < curCost[ secondMinPoint + numLabels * (iNode - 1) ] ) {
                    secondMinPoint = iLabel;
                }
            }

            for(int iLabel = 0; iLabel < numLabels; ++iLabel){
                //  cost for staying in the same position 
    			double tmp1 = curCost[ iLabel + numLabels * (iNode - 1) ];
                
                // cost for coming form the minor second min
                int moveMinPoint = (iLabel != minPoint) ? minPoint : secondMinPoint;
        		double tmp2 =  curCost[ moveMinPoint + numLabels * (iNode - 1)] + pairwiseCost[iNode - 1];
                
            	if( tmp1 < tmp2 ){
    				curCost[iLabel + numLabels * iNode] = tmp1;
        			prevPosition[ iLabel + numLabels * (iNode - 1)] = iLabel;
    			} else {
        			curCost[ iLabel + numLabels * iNode] = tmp2;
    				prevPosition[ iLabel + numLabels * (iNode - 1)] = moveMinPoint;
    			}
    
        		// add the current unary
            	curCost[iLabel + numLabels * iNode] += dataCost[iNode + numNodes * iLabel];
            }
		}
	} else { // numLabels == 1
        for(int iNode = 1; iNode < numNodes; ++iNode) {
            curCost[iNode] = curCost[iNode - 1] + dataCost[iNode];
            prevPosition[iNode - 1] = 0;
        }
    }

	// find best energy
	int minPoint = 0;
	double energy = curCost[minPoint + numLabels * (numNodes - 1)];
	for(int iLabel = 1; iLabel < numLabels; ++iLabel)
		if(curCost[ iLabel + numLabels * (numNodes - 1) ] < energy ){
			minPoint = iLabel;
			energy = curCost[minPoint + numLabels * (numNodes - 1)];
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
		segment[ numNodes - 1 ] = minPoint + 1;
		for(int iNode = numNodes - 2; iNode >= 0; --iNode) {
			minPoint = prevPosition[minPoint + numLabels * iNode];
			segment[ iNode ] = minPoint + 1;
		}
	}

}
