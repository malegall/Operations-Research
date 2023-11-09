#ifndef __Gardien__
#define __Gardien__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include<ilcplex/cplex.h>
#include <assert.h>
#define LARGE_NUMBER 1e9



typedef struct edge 
{
	int id;
	int i;
	int j;
	double cost;

	//Indicates if the edge belongs to the spanning tree
	int is_in_tree;

	//Indicates if the edge has been removed or not
	int displayed;

} edge;

typedef struct node 
{
	int id;

	//Indicates if the node already belongs to the spanning tree
	int is_in_tree;

	//number of incident edges
	int degree;
	//Array of pointers to edges: each one points to an incident edge
	edge** edges;

	//temp index: indicates the current number of neighbors stored
	int tmp_degree;



	//DFS 
	//father indicates the pointer on the edge that goes to the father in the DFS tree
	edge* father;
	int visited;
	int beg;
	int end;


	//Dijkstra
	double d;
	int in_S;//0: untouched, 1: in S, 2: permanently labeled



} node;



typedef struct graph
{
	int n;
	int m;

	edge* edges;
	node* nodes;

	//maximum edge cost
	double max_edge_cost;


	//DFS parameter: current dfs "period"
	int t;


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

//IP Mproblem structure
typedef struct IP_problem
{
	//Internal structures
	CPXENVptr env;
	CPXLPptr lp;
	//number of variables
	int nv;
	//number of constraints
	int nc;

	//Output solution
	double *x;
	//Costs array
	double *cost;
	//TYpe of the variables (binary or continuous)
	char *c_type;
	//Bounds over the variables
	double *up_bound;
	double *low_bound;
	//Names of the variables
	char** var_name;

	//Right hand section of the constraints
	double *rhs;
	//Sense of the constraints
	char *sense;
	//Left hand section for data constraints
	int *rmatbeg;
	int *rmatind;
	double *rmatval;
	int nz;
	//Names of the constraints 
	char** const_name;
	//Solution status
	int solstat;
	double objval;
} IP_problem;




int readTP1Instance(graph* Gptr, char* instanceFileName);
int get_node_by_id(graph*Gptr, int id);
void display_graph(graph*Gptr);
int initialize_tree(tree* Tptr,int n);
int prim(graph* Gptr,int origin_node, tree*Tptr);
int display_tree(tree* Tptr);
int free_graph(graph*Gptr);
int free_tree(tree*Tptr);
int generate_eulerian_graph(graph*Gptr,int n);
int DFS(graph*Gptr,int s,tree*Tptr);
int select_min_d_node(graph* Gptr, node**nptr);
int check_dijkstra_sanity(graph*Gptr);
int dijkstra(graph*Gptr,int s,int t,tree*Tptr);
int check_if_eulerian(graph*Gptr);
int make_graph_eulerian(graph* Gptr);
int min_weight_perfect_matching(graph*Gptr, edge**M);
int in_list(int id, node*list, int n);
int degree_in_tree(int x, tree*Tptr);
int odd_nodes(graph*Gptr, tree*Tptr, node*oddnodes);
void create_complete_graph(graph*Kn, node*oddnodes, int n, tree*Tptr);
void union_matching_tree(graph*Gptr, tree*T, edge**M, int matching_size);
int all_in_tree(graph *Gptr);
int eulerian_prim(graph *Gptr, int origin_node, tree *Tptr);
int alone(graph*Gptr,int s);
int is_bridge(graph *Gptr, edge *e, int s, int*nbnodes);
int empty_graph(graph*Gptr);
void fleury(graph*Gptr, tree*Tptr, int s);
void remove_duplicate_nodes(graph*Gptr,tree*Tptr);
void hamiltonian_cycle(graph * G,tree * cycleE,tree * cycleH);
void write_hamiltonian_cycle(tree*Tptr, int n, char*filename);

#endif