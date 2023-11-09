#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define LARGE_NUMBER 1e9

typedef struct edge edge;


typedef struct node 
{
	int id;

	int nb_neighbors;
	edge** neighbors;

	int cur;

} node;

typedef struct edge 
{
	int id;
	int i;
	int j;
	double cost;
} edge;



typedef struct graph
{
	int n;
	int m;

	edge* edges;
	node* nodes;

	//maximum edge cost
	double max_edge_cost;


} graph;

typedef struct tree 
{
	//The size of the tree is n-1 edges
	int n;

	//array of pointers to the edges in the tree
	edge** edges;

	//Tree total cost
	double cost;
} tree;





int readTP1Instance(graph* Gptr, char* instanceFileName);
int free_graph(graph*Gptr);
int display_graph(graph*Gptr);
int dijkstra(graph*Gptr,int root_node);
int giseven(graph*Gptr);
