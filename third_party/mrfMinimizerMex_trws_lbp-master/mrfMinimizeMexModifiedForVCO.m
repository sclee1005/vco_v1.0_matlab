% Function to optimize MRF energy using TRW-S or BP algorithm (wrapper to Vladimir Kolmogorov's code).
% This version assumes that pairsise potentials can be "decomposed": V_{ij}(k,l) = P(i, j) * M(k, l).
% Here P depends only on variable indeces, M depends only on labels.
% Input examples:
% mrfMinimizeMex(U, P)
% mrfMinimizeMex(U, P, M)
% mrfMinimizeMex(U, P, M, options)
% Output examples:
% S = mrfMinimizeMex(U, P, M, options)
% [S, E] = mrfMinimizeMex(U, P, M, options)
% [S, E, LB] = mrfMinimizeMex(U, P, M, options)
% 
% INPUT:
% 	U		- unary terms (double[numLabels, numNodes])
% 	P		- matrix of edge coefficients (sparse double[numNodes, numNodes]);
% 	M		- matrix of label dependencies (double[naumLabels, numLabels]); if M is not specified, Potts is assumed
% 				if you want to set options without M call: mrfMinimizeMex(U, P, [], options)
% options	- Stucture that determines metod to be used.
% 				Fields:  
% 					method		:	method to use (string: 'trw-s' or 'bp') default: 'trw-s'
% 					maxIter		:	maximum number of iterations (double) default: 100
% 					funcEps		:	If functional change is less than funcEps then stop, TRW-S only (double) default: 1e-2
% 					verbosity	:	verbosity level: 0 - no output; 1 - final output; 2 - full output (double) default: 0
% 					printMinIter:	After printMinIter iterations start printing the lower bound (double) default: 10
% 					printIter	:	and print every printIter iterations (double) default: 5
% 
% OUTPUT: 
% 	S		- labeling that has energy E, vector numNodes * 1 of type double (indeces are in [1,...,numLabels])
% 	E		- The best found energy of type double
% 	LB		- maximum value of lower bound of type double (only for TRW-S method)
%
%	by Anton Osokin, firstname.lastname@gmail.com, Summer 2011 
% modified by syshin for VCO 2017.08


