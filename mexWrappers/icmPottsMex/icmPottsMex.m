% 	icmPottsMex implements the ICM algorithm to minimize an MRF energy with pairwise Potts potentials
% 
% 	Energy function:
% 	E(x)   =   \sum_p D_p(x_p)   +   \sum_pq V_pq [x_p \neq x_q]
% 	where x_p \in {0,...,K},
% 	
% 	Usage:
% 	[energy] = icmPottsMex(unaryTerms, pairwiseTerms, initLabels, maxNumIter);
% 	[energy, newLabels] = icmPottsMex(unaryTerms, pairwiseTerms, initLabels, maxNumIter);
% 	
% 	Inputs:
% 	unaryTerms - of type double, array size [numLabels, numNodes]; 
% 	pairwiseTerms - sparse matrix of type double, size [numNodes, numNodes], only upper triangle is used
% 	initLabels - the initial labeling, double, size [numNodes, 1], default : random
% 	maxNumIter - maximum number of ICM sweeps (default: 10)
% 
% 	Outputs:
% 	energy - of type double, a single number; optimal energy value
% 	newLabels - best found labeling
% 	
% 	Anton Osokin (firstname.lastname@gmail.com),  22.05.2013
