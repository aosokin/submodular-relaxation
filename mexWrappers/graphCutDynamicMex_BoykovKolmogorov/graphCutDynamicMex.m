% 	graphCutDynamicMex - Matlab wrapper to the implementation of min-cut algorithm by Yuri Boykov and Vladimir Kolmogorov:
% 	http://pub.ist.ac.at/~vnk/software/maxflow-v3.03.src.zip
%
%   This version can automatically perform reparametrization on all submodular edges.
% 	This version supports dynamic updates of unary potentials.
% 
%	 Usage:
%	[cut] = graphCutDynamicMex(unaryTerms, pairwiseTerms);
% 	[cut, labels] = graphCutDynamicMex(unaryTerms, pairwiseTerms);
% 	[cut, labels, graphHandle] = graphCutDynamicMex(unaryTerms, pairwiseTerms);
% 
% 	if graphHandle is not requested all memory is cleaned up, otherwise function deleteGraphCutDynamicMex needs to be called
%  
%	Inputs:
% 	termWeights	-	the edges connecting the source and the sink with the regular nodes (array of type double, size : [numNodes, 2])
% 				termWeights(i, 1) is the weight of the edge connecting the source with node #i
% 				termWeights(i, 2) is the weight of the edge connecting node #i with the sink
% 				numNodes is determined from the size of termWeights.
%	edgeWeights	-	the edges connecting regular nodes with each other (array of type double, array size [numEdges, 4])
% 				edgeWeights(i, 3) connects node #edgeWeights(i, 1) to node #edgeWeights(i, 2)
% 				edgeWeights(i, 4) connects node #edgeWeights(i, 2) to node #edgeWeights(i, 1)
%				The only requirement on edge weights is submodularity: edgeWeights(i, 3) + edgeWeights(i, 4) >= 0
% 
% 	Outputs:
% 	cut           -	the minimum cut value (type double)
% 	labels		-	a vector of length numNodes, where labels(i) is 0 or 1 if node #i belongs to S (source) or T (sink) respectively.
% 	graphHandle	- a single number, for direct usage in deleteGraphCutDynamicMex and updateUnaryGraphCutDynamicMex only
%
% 	To build the code in Matlab choose reasonable compiler and run build_graphCutDymanicMex.m
% 	Run example_graphCutDymanicMex.m to test the code
%
%   See also deleteGraphCutDynamicMex, updateUnaryGraphCutDynamicMex
% 
% 	Anton Osokin (firstname.lastname@gmail.com),  19.05.2013
