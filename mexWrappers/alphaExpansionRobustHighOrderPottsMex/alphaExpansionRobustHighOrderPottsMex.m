% alphaExpansionRobustHighOrderPottsMex is a mex-wrapper for the alpha-expansion algorithm to optimize MRF energy with pairwise Potts and robust high-order Potts potentials.
% 
% The original C++ code:
% http://research.microsoft.com/en-us/um/people/pkohli/code/robust.zip
%
% If you use this code please cite the following paper:
%   Pushmeet Kohli, Lubor Ladicky, Philip H. S. Torr. 
%   Robust Higher Order Potentials for Enforcing Label Consistency
%   International Journal of Computer Vision, 82(3):302-324, 2009.
% 
% Input examples:
%   alphaExpansionRobustHighOrderPottsMex(U, P)
%   alphaExpansionRobustHighOrderPottsMexs(U, P, HO, HOC)
%   
% Output examples:
%   [S, E] = alphaExpansionRobustHighOrderPottsMex(U, P, HO, HOC)
%   [S, E, energyPlot, timePlot] = alphaExpansionRobustHighOrderPottsMex(U, P, HO, HOC)
%   
% INPUT:
% 	U		- unary terms (double[numLabels, numNodes])
% 	P		- matrix of edge coefficients (sparse double[numNodes, numNodes]); only upper triangle is used
% 	HO		- groups of edges, showing high-order potentials (cell[numHO, 1], each element - vector of indices)
% 	HOP		- parameters of Robust high-order potentials (double[numHO, 2]), each row gives \gamma_max and Q, \gamma =0; \gamma_max >= 0; Q >= 0;
% 
% OUTPUT: 
%   S		- labeling that has energy E, vector numNodes * 1 of type double (indices are in [1,...,numLabels])
%   E       - energy of labeling S
% 	energyPlot, timePlot - plots of energy and time
% 
%  Anton Osokin (firstname.lastname@gmail.com),  31.05.2013
