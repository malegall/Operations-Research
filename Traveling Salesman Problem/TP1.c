#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "TP1Functions.h"
#include <ilcplex/cplex.h>


int main()
{
	int rval = 0;

	//TP #
	int TPnumber = 6;

	//Structures
	graph G;
	tree T;

	int i,j,n;

	char* instancename = (char*)malloc(256*sizeof(char));
	strcpy(instancename,"K10.csv");

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
				int t = 1;

				//DFS on G starting from s;
				//DFS(&G,s,&T);
				//display_tree(&T);
				//free_tree(&T);



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

				//Freeings
				free_graph(&G);

				break;
			}

			//TP Eulerianize 
		case 6:
			{
				readTP1Instance(&G,instancename);
				
				tree T;
				initialize_tree(&T,G.n);
				
				// On trouve un arbre couvrant de poids minimal
				prim(&G,0,&T);
				display_tree(&T);

				// On sélectionne les sommets de degré impair
				node* oddnodes = (node*)malloc(sizeof(node)*T.n);
				int odd_nodes_size = odd_nodes(&G,&T,oddnodes);
				printf("%d\n",odd_nodes_size);
				for(i=0; i<odd_nodes_size; i++){
					printf("Noeud %d\n", oddnodes[i].id);
				}

				
				// On crée le graphe complet des sommets impairs
				tree* dijkstraT = (tree*)malloc(sizeof(tree)*odd_nodes_size*(odd_nodes_size-1)/2);
				for(i=0; i<odd_nodes_size*(odd_nodes_size-1)/2; i++){
					initialize_tree(&dijkstraT[i],G.n);
				}
				// On cherche le poids à mettre pour chaque arrête
				int k = 0;
				for (i = 0; i < odd_nodes_size; i++){
					for (j = i + 1; j < odd_nodes_size; j++){
						printf("Node %d and %d\n", oddnodes[i].id, oddnodes[j].id);
						dijkstra(&G, oddnodes[i].id, oddnodes[j].id, &dijkstraT[k++]);
					}
				}
				
				// On crée le graphe complet des sommets de degré impair
				graph Kn;
				create_complete_graph(&Kn,oddnodes,odd_nodes_size,dijkstraT);
				display_graph(&Kn);
		
				// On trouve le matching parfait de poids min dans ce graphe complet
				edge**M = (edge**)malloc(sizeof(edge*)*Kn.n/2);
				for(i=0; i<Kn.n/2; i++){
					M[i] = (edge*)malloc(sizeof(edge));
				}
				int matching_size = min_weight_perfect_matching(&Kn,M);
				
				// On unit le matching avec l'arbre couvrant de G
				graph U;
				union_matching_tree(&U,&T,M,matching_size);
				display_tree(&T);
				display_graph(&U);
				
				// On trouve un tour eulérien
				tree fT;
				initialize_tree(&fT,U.m+1);
				fleury(&U,&fT,0);
				display_tree(&fT);
				
				// On construit le cycle hamiltonien
				tree HamiltonC;
				initialize_tree(&HamiltonC,U.n+1);
				hamiltonian_cycle(&U,&fT,&HamiltonC);
				display_tree(&HamiltonC);

				// On écrit le cycle dans le fichier .csv
				char* filename = (char*)malloc(256*sizeof(char));
				strcpy(filename,"sauvegarde.csv");
				write_hamiltonian_cycle(&HamiltonC,U.n,filename);

				//Freeings
				free(oddnodes);		
				free_graph(&G);
				free_graph(&Kn);
				free_graph(&U);
				for(i=0; i<odd_nodes_size*(odd_nodes_size-1)/2; i++){
					free_tree(&dijkstraT[i]);
				}
				free_tree(&T);
				free_tree(&fT);
				free_tree(&HamiltonC);

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

