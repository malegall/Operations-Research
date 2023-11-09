#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "TP1Functions.h"
#include<ilcplex/cplex.h>


int main()
{
	int rval = 0;

	//TP #
	int TPnumber = 6;

	//Structures
	graph G;
	tree T;

	int i,n;

	char* instancename = (char*)malloc(256*sizeof(char));
	strcpy(instancename,"TP4instance.csv");

	//strcpy(instancename,"test.csv");


	switch(TPnumber)
	{
		//TP PRIM: minimum weight spanning trees
		case 1:
			{
				//We read the instance
				readTP1Instance(&G,instancename);

				n = G.n;
				//We initialize the tree structure
				initialize_tree(&T,n);
				//We run Prim from all starting nodes
				for( i = 0 ; i < n ; i++)
				{
					//Prim
					prim(&G,i,&T);
					//We display the tree in the terminal
					display_tree(&T);
				}
				//Freeings
				free_graph(&G);
				free_tree(&T);
				break;
			}
			//TP maximum bipartite matchings
		case 2:
			{
				//We read the instance

				//We initialize//malloc the matching structure

				//We find an initial matching with a greedy algorithm

				//We find a maximum matching from this initial solution

				//Freeings
				//	graph
				//	tree
				break;
			}

			//TP DFS Dijsktra
		case 4:
			{
				readTP1Instance(&G,instancename);

				check_if_eulerian(&G);
				n = G.n;
				//We initialize the tree structure
				initialize_tree(&T,n);

				int s = 0;
				int t = G.n;

				//DFS on G starting from s;
				DFS(&G,s,&T);
				display_tree(&T);
				free_tree(&T);



				//Dijkstra on G between s and t;
				initialize_tree(&T,n);
				dijkstra(&G,s,t,&T);
				display_tree(&T);


				//Freeings
				free_graph(&G);
				free_tree(&T);

				break;
			}
			//TP Eulerian tour
		case 5:
			{
				readTP1Instance(&G,instancename);

				check_if_eulerian(&G);

				//Freeings
				free_graph(&G);

				break;
			}

			//TP Eulerianize 
		case 6:
			{
				readTP1Instance(&G,instancename);

				check_if_eulerian(&G);

				make_graph_eulerian(&G);

				//Freeings
				free_graph(&G);

				break;
			}

			//MISC
		case -1:
			{
				n = 100;
				generate_eulerian_graph(&G,n);
				break;
			}
		default:
			{
				fprintf(stderr,"Wrong TP #. Stopping.\n");
				exit(1);
			}
	}	

	return rval;
}

