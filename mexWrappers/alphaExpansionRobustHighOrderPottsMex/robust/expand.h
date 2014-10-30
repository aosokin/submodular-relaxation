#include <math.h>
#include <time.h>
#include <stdio.h>
#include <sys/timeb.h>
#include "graph.h"
#include "energy.h"

//AOSOKIN
#include <time.h>

// Alpha expansion class
// By Pushmeet Kohli, Lubor Ladicky, Phil Torr

// public functions :

// AExpand(Energy<termType> *e, int maxIter)
// - constructor for the class for solving energy e with maxIter number of iterations 

// void minimize(int *solution)
// - solves energy and saves into the solution array


typedef double termType;
typedef Graph<termType, termType, termType> Grapht;

class AExpand
{
	public:

		//time measurements: Anton
		vector<double> timePlot;
		vector<double> energyPlot;

		AExpand(Energy<termType> *e, int maxIter)
		{
			maxiter = maxIter;
			energy = e;
			nvar = energy->nvar;
			npair = energy->npair;
			nhigher = energy->nhigher;
			nlabel = energy->nlabel;
		}

		termType minimize(int *solution)
		{
			int label_buf_num;
			int step;
			termType E_old;

			//init time measurements: Anton
			timePlot.resize(maxiter * nlabel + 1, -1);
			energyPlot.resize(maxiter * nlabel + 1, -1);
			clock_t tStart = clock();

			label_map = solution;

			g = new Grapht(nvar + 2 * nhigher, npair);
			nodes = new Grapht::node_id[nvar + 2 * nhigher];

			E = compute_energy();

			//update time measurements: Anton
			timePlot[0] = (double)(clock() - tStart) / CLOCKS_PER_SEC;
			energyPlot[0] = E;

			label_buf_num = nlabel;
			step = 0;

			int iter, label;
			for(iter = 0; (iter < maxiter) && (label_buf_num > 0); iter++)
			{
				for(label = 0; label < nlabel; label++)
				{
					E_old = E;
					expand(label);
					g->reset();
					step++;
	        
					E = compute_energy();

					//update time measurements: Anton
					timePlot[1 + label + iter * nlabel] = (double)(clock() - tStart) / CLOCKS_PER_SEC;
					energyPlot[1 + label + iter * nlabel] = E;


					if(E_old == E) label_buf_num--; 
					else label_buf_num = nlabel - 1;
				}
			}
			delete g;
			delete[] nodes;

			return E;
		}

	private :
		int nvar, npair, nhigher, nlabel;
		Grapht *g;
		Grapht::node_id *nodes;
		termType E;

		Energy<termType> *energy;
		int *label_map;
		int maxiter, i, j;

		void expand(unsigned char label)
		{
			termType newE, constE = 0;
			unsigned char label_bar, *is_active;
			
			is_active = new unsigned char[nvar];
			
			for(i = 0; i < nvar; i++)
			{
				label_bar = label_map[i];
				if(label_bar == label)
				{
					is_active[i] = true;
					constE += energy -> unaryCost[i * nlabel + label];
				}
				else
				{
					is_active[i] = false;
					nodes[i] = g -> add_node();
					g->add_tweights(nodes[i], energy->unaryCost[i * nlabel + label], energy->unaryCost[i * nlabel + label_bar]);
				}
			}

			int from, to;
			termType weight;

			for(i = 0; i < npair; i++)
			{
				from = energy->pairIndex[2*i];	
				to = energy->pairIndex[2*i+1];
				weight = energy->pairCost[i];

				if(is_active[from] && is_active[to]) continue;
				else if((is_active[from]) && (!is_active[to])) g->add_tweights(nodes[to], 0, weight);
				else if((!is_active[from]) && (is_active[to])) g->add_tweights(nodes[from], 0, weight);
				else
				{
					if(label_map[from] == label_map[to]) g -> add_edge(nodes[from], nodes[to], weight, weight);
					else
					{
						g->add_tweights(nodes[from], 0, weight);
						g->add_edge(nodes[from], nodes[to], 0, weight);
					}
				}
			}
			termType lambda_a, lambda_b, lambda_m, gamma_b;
			int maxLabel, number_old;

			for(i = 0;i < nhigher; i++)
			{
				maxLabel = getMaxLabel(i);
	            
				lambda_m = energy->higherCost[i * (nlabel + 1) + nlabel];
				lambda_a = energy->higherCost[i * (nlabel + 1) + label];

				nodes[2*i+nvar] = g->add_node();
				g->add_tweights(nodes[2 * i + nvar],0,lambda_m - lambda_a);
				for(j = 0; j < energy->higherElements[i]; j++)
				{
					if (!is_active[energy->higherIndex[i][j]])
						g->add_edge(nodes[2 * i + nvar],nodes[energy->higherIndex[i][j]], 0, (lambda_m - lambda_a) / energy->higherTruncation[i]);
				}
				if((maxLabel == -1) || (maxLabel == label))
				{
					number_old = 0;
					lambda_b = energy->higherCost[i * (nlabel + 1) + nlabel];
				}
				else
				{
					number_old = cardinality(i, maxLabel);
					lambda_b = energy->higherCost[i * (nlabel + 1) + maxLabel] + (energy->higherElements[i] - number_old)
							*(energy->higherCost[i * (nlabel + 1) + nlabel] - energy->higherCost[i * (nlabel + 1) + maxLabel])*(1 / energy->higherTruncation[i]);

					gamma_b = energy->higherCost[i * (nlabel + 1) + maxLabel];

					nodes[2*i+nvar+1] = g->add_node();
					g->add_tweights(nodes[2 * i + nvar + 1],lambda_m - lambda_b,0);
					for(j = 0; j < energy->higherElements[i]; j++)
						if (label_map[energy->higherIndex[i][j]] == maxLabel)
							g->add_edge(nodes[2 * i + nvar + 1],nodes[energy->higherIndex[i][j]], (lambda_m - gamma_b) / energy->higherTruncation[i],0);
				}
				constE -= lambda_m - (lambda_a + lambda_b);
			}

			g -> maxflow();
			for(i = 0; i<nvar; i++) if((!is_active[i]) && (g->what_segment(nodes[i]) == Grapht::SINK)) label_map[i] = label;
			
			newE = compute_energy();
			delete[] is_active;
		}


		int getMaxLabel(int i)
		{
			int j;
			int *num_labels = new int[nlabel];

			for(j = 0;j < nlabel; j++)
				num_labels[j] = 0;

			for(j = 0;j < energy->higherElements[i]; j++)
				num_labels[label_map[energy->higherIndex[i][j]]]++;

			int number=0, maxLabel;

			for(j = 0;j < nlabel; j++)
			{
				if(number <= num_labels[j])
				{
					number = num_labels[j];
					maxLabel = j;
				}
			}

			delete[] num_labels;
			if(number > (energy->higherElements[i] - energy->higherTruncation[i])) return maxLabel;
			else return -1;
		}


		int cardinality(int i, int label)
		{
			int count_label = 0, j;

			for(j = 0;j<energy->higherElements[i]; j++)
				if(label_map[energy->higherIndex[i][j]] == label) 
					count_label++;
			
			return count_label;
		}


		termType compute_energy()
		{
			termType tempE = 0;
			int i, j;
			
			for(i = 0; i < nvar; i++)
				tempE += energy->unaryCost[i * nlabel + label_map[i]];
			
			for(i = 0; i < npair; i++)
			{
				if(label_map[energy->pairIndex[2 * i]] != label_map[energy->pairIndex[2 * i + 1]])
					tempE += energy->pairCost[i];
			}
			
			int *num_labels = new int[nlabel];

			for(i = 0; i < nhigher; i++)
			{
				for(j = 0; j < nlabel; j++) num_labels[j] = 0;

				for(j = 0; j < energy->higherElements[i]; j++) num_labels[label_map[energy->higherIndex[i][j]]]++;

				termType cost, minCost = energy->higherCost[(nlabel + 1) * i + nlabel];

				for(j = 0;j < nlabel; j++)
				{
					cost = energy->higherCost[(nlabel + 1) * i + j] + (energy->higherElements[i] - num_labels[j])
							* (energy->higherCost[(nlabel + 1) * i + nlabel]-energy->higherCost[(nlabel + 1) * i + j]) * (1 / energy->higherTruncation[i]);
					if (minCost >= cost) minCost = cost;
				}
				
				tempE += minCost;
			}
			delete[] num_labels;

			return tempE;
		}
};
