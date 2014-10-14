% 	updateUnaryGraphCutDynamicMex - a part of graphCutDynamicMex:
%		Matlab wrapper to the implementation of min-cut algorithm by Yuri Boykov and Vladimir Kolmogorov:
%	 	http://pub.ist.ac.at/~vnk/software/maxflow-v3.03.src.zip
% 	
%	updateUnaryGraphCutDynamicMex updates the terminal weights and computes the new cut value
% 
%	Usage:
%	[cut] = updateUnaryGraphCutDynamicMex(graphHandle, changedVertices);
%	[cut, labels] = updateUnaryGraphCutDynamicMex(graphHandle, changedVertices);
%  
%	Inputs:
%	graphHandle - a single number given by graphCutDynamicMex
%	updateUnary - of type double, array size [numChanges, 3];  ([p, sourceLink, sinkLink]); the extra cost of the terminal links of node #p
% 
%	Outputs:
%	cut         -	the minimum cut value (type double)
%	labels		-	a vector of length numNodes, where labels(i) is 0 or 1 if node #i belongs to S (source) or T (sink) respectively.
% 
%	See also deleteGraphCutDynamicMex, graphCutDynamicMex
% 
% 	Anton Osokin (firstname.lastname@gmail.com),  19.05.2013
