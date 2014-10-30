#ifndef __ENERGY_H__
#define __ENERGY_H__


// Energy class
// By Pushmeet Kohli, Lubor Ladicky, Phil Torr

// public variables :

// - nlabel - number of labels
// - nvar - number of nodes/variables
// - npair - number of pairwise connections
// - nhigher - number of cliques/segments
// - unaryCost - array of unary costs [nvar * nlabel]
// - pairIndex - array of pairs of indexes of nodes for each pairwise connection [2 * npair]
// - pairCosts - array of costs for pairwise Potts model [npair]
// - higherElements - array of number of elements in each segment [nhigher]
// - higherCosts - array of higher order costs {gamma_k, gammamax} for each segment [nhigher * (nlabel + 1)]
// - higherTruncation - array of Truncation parameter Q for each segment [nhigher]
// - higherIndex - array of arrays of indexes of nodes for each segment [nhigher] x [higherElements[]]


// public functions :

// Energy(int nLabel, int nVar, int nPair, int nHigher)
// - constructor for the class initializing number of labels, nodes, pairwise connections and segments

// void AllocateHigherIndexes()
// - allocates memory for list of indexes within each clique
// - this function has to be called after number of nodes in each segment (higherElements) is initialized


template<typename termType> class Energy
{
	public:
		termType *unaryCost, *pairCost, *higherCost, *higherTruncation;
		int nlabel, nvar, npair, nhigher;
		int *pairIndex, **higherIndex, *higherElements;
		
		Energy(int nLabel, int nVar, int nPair, int nHigher)
		{
			nlabel = nLabel;
			nvar = nVar;
			npair = nPair;
			nhigher = nHigher;

			unaryCost = new termType[nvar * nlabel];
			pairIndex = new int[npair * 2];
			pairCost = new termType[npair];

			higherCost = new termType[nhigher * (nlabel + 1)];
			higherElements = new int[nhigher];
			higherTruncation = new double[nhigher];
			higherIndex = new int *[nhigher];
			memset(higherIndex, 0, nhigher * sizeof(int *));
		}

		~Energy()
		{
			int i;
			if(higherIndex != NULL) for(i = 0; i < nhigher; i++) if(higherIndex[i] != NULL) delete[] higherIndex[i];
	
			if(pairCost != NULL) delete[] pairCost;
			if(unaryCost != NULL) delete[] unaryCost;
			if(pairIndex != NULL) delete[] pairIndex;

			if(higherTruncation != NULL) delete[] higherTruncation;
			if(higherCost != NULL) delete[] higherCost;
			if(higherIndex != NULL) delete[] higherIndex;
			if(higherElements != NULL) delete[] higherElements;
		}

		void AllocateHigherIndexes()
		{
			int i;
			for(i = 0;i < nhigher; i++) higherIndex[i] = new int[higherElements[i]];
		}
};

#endif
