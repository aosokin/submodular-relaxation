#include "expand.h"
#include <stdio.h>
#include <stdlib.h>

int main()
{
	int i, j;

	//initialize random seed
	srand(40000);

	//inititialize energy with 3 labels, 25 variables (grid 5x5), 40 pairwise connections, 5 segments
	int gridsize = 5, segments = 5;
	Energy<double> *energy = new Energy<double>(3, gridsize * gridsize, 2 * gridsize * (gridsize - 1), segments);

	//initialize random unary costs (variables x labels)
	for(i = 0; i < energy->nvar; i++) for(j = 0; j < energy->nlabel; j++)
	{
		energy->unaryCost[i * energy->nlabel + j] = (rand() % 100) * 0.1;
	}

	//initialize grid topology - horizontal connections
	for(i = 0; i < gridsize; i++) for(j = 0; j < gridsize - 1; j++)
	{
		energy->pairIndex[(i * (gridsize - 1) + j) * 2] = i * gridsize + j;
		energy->pairIndex[(i * (gridsize - 1) + j) * 2 + 1] = i * gridsize + j + 1;
	}

	//initialize grid topology - vertical connections
	int offset = gridsize * (gridsize - 1) * 2;
	for(i = 0; i < gridsize; i++) for(j = 0; j < gridsize - 1; j++)
	{
		energy->pairIndex[offset + (i * (gridsize - 1) + j) * 2] = j * gridsize + i;
		energy->pairIndex[offset + (i * (gridsize - 1) + j) * 2 + 1] = (j + 1) * gridsize + i;
	}

	//initialize random pariwise costs
	for(i = 0; i < energy->npair; i++)
	{
		energy->pairCost[i] =  (rand() % 20) * 0.1;
	}

	//initialize number of elements in each segment
	energy->higherElements[0] = 5;
	energy->higherElements[1] = 3;
	energy->higherElements[2] = 8;
	energy->higherElements[3] = 4;
	energy->higherElements[4] = 5;

	//allocate energy for higher order indexes
	energy->AllocateHigherIndexes();

	//initialize indexes of pixels in each segment
	energy->higherIndex[0][0] = 0;
	energy->higherIndex[0][1] = 1;
	energy->higherIndex[0][2] = 2;
	energy->higherIndex[0][3] = 5;
	energy->higherIndex[0][4] = 6;

	energy->higherIndex[1][0] = 3;
	energy->higherIndex[1][1] = 4;
	energy->higherIndex[1][2] = 8;

	energy->higherIndex[2][0] = 7;
	energy->higherIndex[2][1] = 9;
	energy->higherIndex[2][2] = 12;
	energy->higherIndex[2][3] = 13;
	energy->higherIndex[2][4] = 14;
	energy->higherIndex[2][5] = 18;
	energy->higherIndex[2][6] = 19;
	energy->higherIndex[2][7] = 24;

	energy->higherIndex[3][0] = 10;
	energy->higherIndex[3][1] = 11;
	energy->higherIndex[3][2] = 15;
	energy->higherIndex[3][3] = 20;

	energy->higherIndex[4][0] = 16;
	energy->higherIndex[4][1] = 17;
	energy->higherIndex[4][2] = 21;
	energy->higherIndex[4][3] = 22;
	energy->higherIndex[4][4] = 23;

	//initialize truncation ratio Q, gamma_k and gamma_max for each clique
	for(i = 0; i < energy->nhigher; i++)
	{
		//truncation ratio 30%
		energy->higherTruncation[i] = 0.3;

		//gamma_k
		for(j = 0; j < energy->nlabel; j++) energy->higherCost[i * (energy->nlabel + 1) + j] = (rand() % 100) * 0.1;

		//gamma_max
		energy->higherCost[i * (energy->nlabel + 1) + energy->nlabel] = 100 * 0.1;
	}

	//initialize alpha expansion - max 10 iterations
	AExpand *expand = new AExpand(energy, 10);

	//initialize solution
	int *solution = new int[gridsize * gridsize];
	memset(solution, 0, gridsize * gridsize * sizeof(int));

	//solve CRF
	expand->minimize(solution);

	//print solution
	for(i = 0; i < gridsize; i++)
	{
		for(j = 0; j < gridsize; j++) printf("%d ", solution[i * gridsize + j]);
		printf("\n");
	}

	//free memory
	delete[] solution;
	delete expand;
	delete energy;

	return 0;
}
