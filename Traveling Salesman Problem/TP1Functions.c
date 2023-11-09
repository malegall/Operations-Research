#include "TP1Functions.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include<stdio.h>
#include<ilcplex/cplex.h>

int readTP1Instance(graph* Gptr, char* instanceFileName)
{
	int rval = 0;

	//We open the file
	FILE* fin = fopen(instanceFileName,"r");

	//We read the number of nodes and vertices
	int n,m;
	int e,i,j;

	rval = fscanf(fin,"%d,%d\n", &n,&m);
	Gptr->n = n;
	Gptr->m = m;
	fprintf(stderr,"n = %d\t m = %d\n",Gptr->n,Gptr->m);

	//We malloc and fill the nodes in their structure
	Gptr->nodes= (node*)malloc(sizeof(node)*n);
	node* iptr;
	for( i = 0 ;  i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		iptr->id = i;
		//The node is initially not in the tree
		iptr->is_in_tree = 0;
		//The node has initially degree zero
		iptr->degree= 0;
		iptr->tmp_degree= 0;
	}


	//We malloc, read and fill the edges in their structure
	Gptr->edges = (edge*)malloc(sizeof(edge)*m);
	double cost;
	edge* eptr;
	double max_edge_cost = -LARGE_NUMBER;
	for( e = 0 ;  e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		rval=fscanf(fin,"%d,%d,%lf\n", &i,&j,&cost);
		eptr->id = e;
		eptr->i = i;
		eptr->j = j;
		eptr->cost = cost;
		//We update the largest edge cost
		if(max_edge_cost < cost)
			max_edge_cost = cost;
		//The edge is initially not in the tree
		eptr->is_in_tree = 0;
		eptr->displayed = 1;
		//We update the degree of the edge's endpoints
		//tail node
		iptr = &(Gptr->nodes[i]);
		iptr->degree++;
		//head node
		iptr = &(Gptr->nodes[j]);
		iptr->degree++;

		//fprintf(stderr,"%d,%d,%d,%lf\n",eptr->ideptr->i,eptr->j,eptr->cost);
	}
	//We store the largest edge cost in the graph structure
	Gptr->max_edge_cost = max_edge_cost;


	//We close the file
	fclose(fin);

	//We fill in the neighbors for each node:
	//Note that we are working with undirected graphs
	//We malloc their arrays
	int degree;
	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		degree = iptr->degree;
		iptr->edges = (edge**)malloc(sizeof(edge*)*degree);
	}
	//We loop over the edges to fill the neighbors arrays
	int tmp;
	for( e = 0 ;  e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		i = eptr->i;
		j = eptr->j;

		//tail node
		iptr = &(Gptr->nodes[i]);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;

		//head node
		iptr = &(Gptr->nodes[j]);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;

		//fprintf(stderr,"%d,%d,%d,%lf\n",eptr->ideptr->i,eptr->j,eptr->cost);
	}

	return rval;
}

int get_node_by_id(graph*Gptr, int id){
	int i = 0;
	int b = 0;
	int rval = 0;
	node*nptr;
	while(i<Gptr->n && !b){
		nptr = &(Gptr->nodes[i]);
		if(nptr->id == id){
			rval = i;
			b = 1;
		}
	}
	return rval;
}

void display_graph(graph*Gptr){
	int i,e,n,m;

	n = Gptr->n;
	m = Gptr->m;
	fprintf(stderr,"The graph has %d nodes and %d edges:\n",n,m);
	fprintf(stderr,"Displaying nodes:\n");

	node*iptr;
	edge*eptr;

	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		fprintf(stderr,"\tNode %d has %d neighbors:\n",iptr->id, iptr->degree);

		for( e = 0 ; e < iptr->degree ; e++)
		{
			eptr = iptr->edges[e];
			if(eptr->displayed){
				fprintf(stderr,"\t\t(%d,%d)\t[w=%lf]\n",eptr->i,eptr->j,eptr->cost);
			}
			
		}
	}
}

int initialize_tree(tree* Tptr,int n)
{
	int rval = 0;

	Tptr->n = n;
	Tptr->edges = (edge**)malloc(sizeof(edge*) * (n-1));
	Tptr->cost = 0;


	return rval;
}

int prim(graph* Gptr, int origin_node, tree*Tptr)
{
	int rval = 0;

	int m = Gptr->m;
	int n = Gptr->n;
	int e,i,j;

	//We reinitialize the edges and nodes to not be in the tree
	for( e = 0 ;  e < m ; e++)
		Gptr->edges[e].is_in_tree = 0;
	for( i = 0 ;  i < n ; i++)
		Gptr->nodes[i].is_in_tree = 0;
	//We mark the origin node as in the tree
	Gptr->nodes[origin_node].is_in_tree = 1;

	fprintf(stderr,"Starting PRIM's algorithm from node %d\n",origin_node);



	//Prim's main body
	double mincost;
	double maxcost = Gptr->max_edge_cost;
	edge*eptr;
	node*iptr,*jptr;
	int candidate_e;
	int tree_size = 0;
	for(tree_size = 0 ; tree_size < n-1 ; tree_size++)
	{
		mincost = maxcost; 
		candidate_e = -1;
		//We look for the edge (i,j) with minimum weight such that 
		//	(i,j) not already in the tree
		//	ANDi NOT
		//	  i and j both in the tree
		//	OR  i and j both NOT in the tree 

		for( e = 0 ;  e < m ; e++)
		{
			//We fetch the pointer to the edge
			eptr = &(Gptr->edges[e]);

			//If it is already in the tree, skip it
			if(eptr->is_in_tree)
				continue;

			//We fetch the pointers to the endpoint nodes
			i = eptr->i;
			iptr = &(Gptr->nodes[i]);
			j = eptr->j;
			jptr = &(Gptr->nodes[j]);

			//If i and j are in the tree, skip it
			if(	(iptr->is_in_tree) && (jptr->is_in_tree)	)
				continue;

			//If i and j are both outside the tree skip it
			if(	!(iptr->is_in_tree) && !(jptr->is_in_tree)	)
				continue;

			//Otherwise, we consider it and check its cost
			if(eptr->cost < mincost)
			{
				mincost = eptr->cost;	
				candidate_e = eptr->id;;
			}
		}

		//If we did not find a candidate, the graph is not connected 
		if(candidate_e == -1)
		{
			fprintf(stderr,"The graph is disconnected. Aborting Prim\n");
		}
		//Otherwise, we add the least cost edge we found to the tree
		//We fetch the pointer to the edge
		eptr = &(Gptr->edges[candidate_e]);
		//We mark it as in the tree
		eptr->is_in_tree = 1;

		//We fetch the pointers to the endpoint nodes
		i = eptr->i;
		iptr = &(Gptr->nodes[i]);
		j = eptr->j;
		jptr = &(Gptr->nodes[j]);
		//And we mark both nodes as marked (one of them already is)
		iptr->is_in_tree = 1;
		jptr->is_in_tree = 1;
	}	

	//We fill the tree 
	Tptr->cost = 0;
	int count = 0;
	//We fill the tree structure
	for( e = 0 ; e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		//If the edge is not marked, we skip it
		if(!eptr->is_in_tree)
			continue;
		//Otherwise, we put it in the tree
		Tptr->edges[count++] = eptr;
		Tptr->cost += eptr->cost;
	}

	return rval;
}

int display_tree(tree* Tptr)
{
	int rval = 0;
	int n = Tptr->n;
	int e;
	edge*eptr;

	fprintf(stderr,"XXXXXXXX tree XXXXXXXX\n");
	fprintf(stderr,"SIZE:\t%d\n",Tptr->n);
	fprintf(stderr,"WEIGHT:\t%lf\n",Tptr->cost);
	for( e = 0 ; e < n-1 ; e++)
	{
		eptr = Tptr->edges[e];
		fprintf(stderr,	"Edge #%d (%d,%d)\tcost:\t%lf\n",
				eptr->id,eptr->i, eptr->j, eptr->cost);	
		
	}
	fprintf(stderr,"XXXXXXXXXXXXXXXXXXXXXX\n\n\n");

	return rval;
}


int free_graph(graph*Gptr)
{
	int rval=0;

	int n = Gptr->n;
	int i;
	node* iptr;
	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		free(iptr->edges);
	}
	free(Gptr->edges);
	free(Gptr->nodes);

	return rval;
}

int free_tree(tree*Tptr)
{
	int rval=0;

	free(Tptr->edges);

	return rval;
}


int generate_eulerian_graph(graph*Gptr,int n)
{
	int rval = 0;

	int** incidence = (int**)malloc(sizeof(int*)*n);
	int i,j;
	for( i = 0 ; i < n ; i++)
	{
		incidence[i] = (int*)malloc(sizeof(int)*n);
		for( j = 0 ; j < n ; j++)
			incidence[i][j] = 0;
		//There is a cycle linking everybody: the graph is connected
		if(i<n-1)
			incidence[i][i+1] = 1;
	}
	//We close the cycle
	incidence[0][n-1] = 1;

	//We keep the degrees
	int * degree = (int*)malloc(sizeof(int)*n);
	for( i = 0 ; i < n ; i++)
		degree[i] = 2;


	//We currently have n edges
	int m = n;
	//Number of edges to reach
	int M = 7*n;

	int e, connected;
	//We randomly add edges until we reach M
	for( e = n ; e < M ; e++)
	{
		//We draw a non connected pair
		do
		{
			i = rand() % n;
		}while(degree[i] == n-1);

		do
		{
			j = rand() % n;
			//If already connected, we pass
			if(incidence[i][j] || incidence[j][i])
				connected = 1;
			else connected = 0;
		}while((degree[j] == n-1) || (i == j) || (connected));

		degree[i]++;
		degree[j]++;
		incidence[i][j] = 1;
		m++;
	}

	//Then we add edges until all degrees are even
	int odd = 1;
	int testodd;
	do
	{
		odd = 0;
		//We look for an odd degree pair
		for( i = 0 ; i < n ; i++)
		{
			testodd = degree[i]/2;
			testodd = 2*testodd;	
			if(testodd != degree[i])
			{
				odd = 1;
				break;
			}
		}
		if(odd == 0)
			break;
		odd = 0;
		//We look for another odd degree vertex
		for( j = i+1 ; j < n ; j++)
		{
			if(incidence[i][j] || incidence[j][i])
				continue;
			testodd = degree[j]/2;
			testodd = 2*testodd;	
			if(testodd != degree[j])
			{
				odd = 1;
				break;
			}
		}
		if(odd == 0)
		{
			fprintf(stderr,"Impossible, there sould be at least two odd degree vertices\n");
			exit(0);
		}
		degree[i]++;
		degree[j]++;
		incidence[i][j] = 1;
		m++;
	}while(odd);

	//sanity check
	for( i = 0 ; i < n ; i++)
		for( j = 0 ; j < n ; j++)
			if(incidence[i][j] && incidence[j][i])
			{
				fprintf(stderr,"Symmetry\n");	
				exit(0);
			}
	int deg;
	for( i = 0 ; i < n ; i++)
	{
		deg=0;
		for( j = 0 ; j < n ; j++)
			if(incidence[i][j] || incidence[j][i])
				deg++;

		if(deg != 2*(deg/2))
		{
			fprintf(stderr,"Node %d has odd degree %d\n",i,deg);
			exit(0);
		}
	}







	//we scramble all indices;
	int*correspondance = (int*)malloc(sizeof(int)*n);
	int*taken = (int*)malloc(sizeof(int)*n);

	for( i = 0 ; i < n ; i++)
	{
		correspondance[i] = -1;
		taken[i] = 0;
	}

	for( i = 0 ; i < n ; i++)
	{
		do
		{
			j = rand() % n;

		}while(taken[j]);
		taken[j] = 1;
		correspondance[i] = j;
	}


	int cost;
	FILE* fin = fopen("TP5instance.csv","w");

	fprintf(fin,"%d,%d\n",n,m);
	for( i = 0 ; i < n ; i++)
		for( j = 0 ; j < n ; j++)
		{
			if(incidence[i][j] == 0)
				continue;

			cost = rand() % 15;
			cost++;
			fprintf(fin,"%d,%d,%d\n",correspondance[i],correspondance[j],cost);
		}


	fclose(fin);

	free(degree);
	free(correspondance);
	free(taken);

	for( i = 0 ; i < n ; i++)
		free(incidence[i]);

	free(incidence);


	return rval;
}

int DFS_explore(graph*Gptr, int s)
{
	int rval = 0;

	node*iptr = &(Gptr->nodes[s]);
	iptr->visited = 1;
	Gptr->t++;
	iptr->beg=Gptr->t;
	int degree = iptr->degree;
	int e;
	edge*eptr;
	int i,j;
	node*jptr;

	for( e  = 0 ;  e  < degree ; e++)
	{
		eptr = iptr->edges[e];
		if(eptr->displayed){
			//printf("Edge #%d ((%d,%d) [%lf])\n",eptr->id,eptr->i,eptr->j,eptr->cost);
			i = eptr->i;
			j = eptr->j;
			//We select the neighboring node
			if(i != s)
				j=i;
			jptr = &(Gptr->nodes[j]);

			//If neighbor already visited, skip it
			if(jptr->visited)
				continue;
			//Otherwise, we explore it
			jptr->visited = 1;
			jptr->father = eptr;
			DFS_explore(Gptr,j);
		}
		
	}
	iptr->visited = 2;
	Gptr->t++;
	iptr->end = Gptr->t;

	return rval;
}



int DFS(graph*Gptr,int s, tree*Tptr)
{
	int rval = 0;

	//We reinitialize the DFS parameters
	Gptr->t = 0;

	int n = Gptr->n;
	int i;
	node*iptr;
	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		iptr->father = NULL;
		iptr->visited= 0;
		iptr->beg= -1;
		iptr->end= -1;
	}

	i = s;
	int j;
	while(i != -1)
	{
		//We explore the node 
		iptr = &(Gptr->nodes[i]);
		if(iptr->visited == 0)
			DFS_explore(Gptr,i);

		//We pick the next node to explore
		i=-1;
		for( j = 0 ; j < n ; j++)
		{
			iptr = &(Gptr->nodes[j]);
			if(iptr->visited)
				continue;
			i=j;
			break;
		}
	}

	//We now find the set of edges composing the tree
	int cnt = 0;
	Tptr->cost = 0;
	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		//If we are at the root
		if(iptr->father == NULL)
			continue;

		Tptr->edges[cnt++]=iptr->father;
		Tptr->cost += iptr->father->cost;
	}
	//fprintf(stderr,"n = %d, cnt = %d\n",n,cnt);
	//exit(0);




	return rval;
}

int select_min_d_node(graph* Gptr, node**nptr)
{
	int rval = 0;

	int i;
	int n = Gptr->n;

	double min_d = LARGE_NUMBER;; 
	node*iptr;

	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		if(iptr->in_S != 1)
			continue;

		if(iptr->d > min_d)
			continue;

		min_d = iptr->d;
		*nptr = iptr;		
	}

	return rval;
}

int check_dijkstra_sanity(graph*Gptr)
{
	int rval = 0;

	int m = Gptr->m;
	int e,i,j;
	edge*eptr;
	node*iptr,*jptr;

	int sanity = 1;

	for( e = 0 ; e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		i = eptr->i;
		iptr = &(Gptr->nodes[i]);
		j = eptr->j;
		jptr = &(Gptr->nodes[j]);

		if(iptr->d + eptr->cost < jptr->d)
		{
			fprintf(stderr,
					"Dijkstra failed for edge %d between nodes %d and %d: %lf+%lf = %lf < %lf\n",
					e,i,j,iptr->d,eptr->cost,iptr->d + eptr->cost,jptr->d);
			sanity = 0;
		}
		if(jptr->d + eptr->cost < iptr->d)
		{
			fprintf(stderr,
					"Dijkstra failed for edge %d between nodes %d and %d: %lf+%lf = %lf < %lf\n",
					e,j,i,jptr->d,eptr->cost,jptr->d + eptr->cost,iptr->d);
			sanity = 0;
		}
	}
	if(sanity)
		fprintf(stderr,"Potentials are consistent. Disjktra ended successfully.\n");
	else
		fprintf(stderr,"Potentials are NOT consistent. Disjktra failed.\n");

	return rval;
}

int all_in_tree(graph *Gptr)
{
	for (int i = 0; i < Gptr->n; i++)
	{
		node *iptr = &(Gptr->nodes[i]);
		if (iptr->in_S != 2)
		{
			return 0;
		}
	}

	return 1;
}

int dijkstra(graph *Gptr, int s, int t, tree *Tptr)
{
	int rval = 0;

	int n = Gptr->n;
	int i, j;
	node *iptr, *jptr;
	for (i = 0; i < n; i++)
	{
		iptr = &(Gptr->nodes[i]);
		iptr->father = NULL;
		iptr->d = LARGE_NUMBER;
		//All nodes are outside S at first
		iptr->in_S = 0;
	}

	i = s;
	iptr = &(Gptr->nodes[i]);
	iptr->father = NULL;
	iptr->d = 0;
	iptr->in_S = 1;

	edge * eptr;
	int e;
	int m = Gptr->m;
	for (e = 0; e < m; e++)
	{
		eptr = &(Gptr->edges[e]);
		eptr->is_in_tree = 0;
	}

	//Dijkstra loop
	while (1)
	{
		//We select the node with lowest label in S
		iptr = NULL;
		select_min_d_node(Gptr, &iptr);
		if (iptr->id == t)
		{
			//iptr->id == t
			Tptr->cost = iptr->d;
			break;
		}

		/*fprintf(stderr,"Node %d (%lf) selected",iptr->id, iptr->d);
		if(iptr->father !=NULL)
		fprintf(stderr,"(father: %d)",iptr->father->id);
		fprintf(stderr,"\n");
		*/

		i = iptr->id;

		//We delete it from S
		iptr->in_S = 2;

		//We check all incident edges
		for (e = 0; e < iptr->degree; e++)
		{
			eptr = iptr->edges[e];
			//We catch the index of the neighboring node
			if (eptr->i == i)
				j = eptr->j;
			else
				j = eptr->i;
			//and its pointer
			jptr = &(Gptr->nodes[j]);
			//fprintf(stderr,"\tChecking node %d (%lf) via edge %d (%lf)\n",
			// j,jptr->d,eptr->id, eptr->cost);

			//if already permanently labeled, skip it
			if (jptr->in_S == 2)
				continue;

			//If untouched yet, add it to S
			if (jptr->in_S == 0)
				jptr->in_S = 1;

			//If improvable, do it
			if (jptr->d > iptr->d + eptr->cost)
			{
				if (jptr->father != NULL)
					jptr->father->is_in_tree = 0;
				jptr->father = eptr;
				jptr->father->is_in_tree = 1;
				jptr->d = iptr->d + eptr->cost;
			}
		}
	}

	check_dijkstra_sanity(Gptr);

	//We fill the tree 

	int count = 0;
	//We fill the tree structure
	for (e = 0; e < m; e++)
	{
		eptr = &(Gptr->edges[e]);
		//If the edge is not marked, we skip it
		if (!eptr->is_in_tree)
			continue;
		//Otherwise, we put it in the tree
		Tptr->edges[count++] = eptr;
	}

	return rval;
}

int check_if_eulerian(graph*Gptr)
{
	int rval = 0;

	int n = Gptr->n;
	int i;
	node*iptr;
	int degree;
	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		degree = iptr->degree;
		if(degree != 2*(degree/2))
		{
			fprintf(stderr,"The graph is not eulerian: node %d has odd degree %d\n",
					i, degree);
			return rval;
		}
	}
	fprintf(stderr,"The graph is Eulerian\n");

	return rval;
}

int min_weight_perfect_matching(graph*Gptr, edge**M)
{
	int rval = 0;

	IP_problem ip_prob;
	IP_problem* ip_prob_ptr = &ip_prob;
	ip_prob_ptr->env = NULL;
	ip_prob_ptr->lp = NULL;
	ip_prob_ptr->env = CPXopenCPLEX (&rval);
	if(rval) fprintf(stderr,"ERROR WHILE CALLING CPXopenCPLEX\n");
	if ( ip_prob_ptr->env == NULL ) 
	{
		char  errmsg[1024];
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (ip_prob_ptr->env, rval, errmsg);
		fprintf (stderr, "%s", errmsg);
		exit(0);	
	}

	//We create the MIP problem
	ip_prob_ptr->lp = CPXcreateprob (ip_prob_ptr->env, &rval, "TP1");
	if(rval) fprintf(stderr,"ERROR WHILE CALLING CPXcreateprob\n");

	rval = CPXsetintparam (ip_prob_ptr->env, CPX_PARAM_DATACHECK, CPX_ON); 
	rval = CPXsetintparam (ip_prob_ptr->env, CPX_PARAM_SCRIND, CPX_ON);

	//Number of nodes in the graph
	int n = Gptr->n;
	//Number of edges in the graph
	int m = Gptr->m;
	//Number of variables (=nb of edges)
	int nv = m;

	//We fill our arrays
	//Memory
	ip_prob_ptr->nv = nv;
        ip_prob_ptr->x = (double*)malloc(sizeof(double)*nv);
        ip_prob_ptr->cost = (double*)malloc(sizeof(double)*nv);
        ip_prob_ptr->c_type = (char*)malloc(sizeof(char)*nv);
        ip_prob_ptr->up_bound = (double*)malloc(sizeof(double)*nv);
        ip_prob_ptr->low_bound = (double*)malloc(sizeof(double)*nv);
	ip_prob_ptr->var_name = (char**)malloc(sizeof(char*)*nv);



	int e,i,id = 0;
	//Structures keeping the index of each variable
	int*id_x_e = (int*)malloc(sizeof(int)*m);

	edge* eptr;
	for( e = 0 ; e < m ; e++)
	{
		//We keep the id
		id_x_e[e] = id;

		//we fetch the pointer on  the current edge e
		eptr = &(Gptr->edges[e]);
		//We generate the variable attributes
		ip_prob_ptr->x[id] = 0;
		ip_prob_ptr->cost[id] = eptr->cost;
		ip_prob_ptr->c_type[id] = 'B'; //binary
		ip_prob_ptr->up_bound[id] = 1;
		ip_prob_ptr->low_bound[id] = 0;
		ip_prob_ptr->var_name[id] = (char*)malloc(sizeof(char)*1024);
	        snprintf(       ip_prob_ptr->var_name[id],
        	                1024,
                	        "x_e%d",
                        	e);
		id++;
	}


	rval = CPXnewcols( ip_prob_ptr->env, ip_prob_ptr->lp, 
			nv, 
			ip_prob_ptr->cost, 
			ip_prob_ptr->low_bound,
			ip_prob_ptr->up_bound,
			ip_prob_ptr->c_type,
			ip_prob_ptr->var_name);
	if(rval)
		fprintf(stderr,"CPXnewcols returned errcode %d\n",rval);



	//Constraints part
        ip_prob_ptr->rhs = (double*)malloc(sizeof(double));
        ip_prob_ptr->sense = (char*)malloc(sizeof(char));
        ip_prob_ptr->rmatbeg = (int*)malloc(sizeof(int));
	ip_prob_ptr->nz = m;	//maximum size for nonzeroes in a constraint


        ip_prob_ptr->rmatind = (int*)malloc(sizeof(int)*nv);
        ip_prob_ptr->rmatval = (double*)malloc(sizeof(double)*nv);
	ip_prob_ptr->const_name = (char**)malloc(sizeof(char*));
	ip_prob_ptr->const_name[0] = (char*)malloc(sizeof(char)*1024);

	//We fill what we can 
	ip_prob_ptr->rmatbeg[0] = 0;

	//We generate and add each constraint to the model
	//Perfect matching constraints
	ip_prob_ptr->rhs[0] = 1;
	ip_prob_ptr->sense[0] = 'E'; //equal
	node* ndptr;
	int nz;
	for( i = 0 ; i < n ; i++)
	{
		//Constraint name
	        snprintf(       ip_prob_ptr->const_name[0],
        	                1024,
                	        "node_matched_i%d",
                        	i);
		id=0;

		//We pick a pointer to the node #i
		ndptr = &(Gptr->nodes[i]);
		//The number of nonzeroes in this constraint 
		//is equal to the number of incident edges
		nz = ndptr->degree;
		//variables x_ij coefficients
		for( e = 0 ; e < nz ; e++)
		{
			//We fetch a pointer on the current edge
			eptr = ndptr->edges[e];
		        ip_prob_ptr->rmatind[id] = id_x_e[eptr->id];
        		ip_prob_ptr->rmatval[id] =  1;
			id++;
		}
		rval = CPXaddrows( ip_prob_ptr->env, ip_prob_ptr->lp, 
			0,//No new column
			1,//One new row
			nz,//Number of nonzero coefficients
			ip_prob_ptr->rhs, 
			ip_prob_ptr->sense, 
			ip_prob_ptr->rmatbeg, 
			ip_prob_ptr->rmatind, 
			ip_prob_ptr->rmatval,
			NULL,//No new column
			ip_prob_ptr->const_name );
		if(rval)
			fprintf(stderr,"CPXaddrows returned errcode %d\n",rval);

	}


	//We write the problem for debugging purposes, can be commented afterwards
	rval = CPXwriteprob (ip_prob_ptr->env, ip_prob_ptr->lp, "min_weight_perfect_matching.lp", NULL);
	if(rval)
		fprintf(stderr,"CPXwriteprob returned errcode %d\n",rval);

	//We solve the model
	rval = CPXmipopt (ip_prob_ptr->env, ip_prob_ptr->lp);
	if(rval)
		fprintf(stderr,"CPXmipopt returned errcode %d\n",rval);

	rval = CPXsolwrite( ip_prob_ptr->env, ip_prob_ptr->lp, "min_weight_perfect_matching.sol" );
	if(rval)
		fprintf(stderr,"CPXsolwrite returned errcode %d\n",rval);

	//We get the objective value
	rval = CPXgetobjval( ip_prob_ptr->env, ip_prob_ptr->lp, &(ip_prob_ptr->objval) );
	if(rval)
		fprintf(stderr,"CPXgetobjval returned errcode %d\n",rval);

	//We get the best solution found 
	rval = CPXgetobjval( ip_prob_ptr->env, ip_prob_ptr->lp, &(ip_prob_ptr->objval) );
	rval = CPXgetx( ip_prob_ptr->env, ip_prob_ptr->lp, ip_prob_ptr->x, 0, nv-1 );
	if(rval)
		fprintf(stderr,"CPXgetx returned errcode %d\n",rval);

	
	//We display the solution
	double tolerance = 0.0001;
	for( e = 0 ; e < m ; e++) 
	{
		id = id_x_e[e];
		if(ip_prob_ptr->x[id] > 1-tolerance){	
			eptr = &(Gptr->edges[e]);
			fprintf(stderr,"Edge #%d ((%d,%d) [%lf]) is used in the perfect matching\n",
					e,eptr->i,eptr->j,eptr->cost);
		}
	}

	//Build the matching here from the CPLEX solution
	int k = 0;
	for(e=0; e<m; e++){
		id = id_x_e[e];
		if(ip_prob_ptr->x[id] <= 1-tolerance)
			continue;
		eptr = &(Gptr->edges[e]);
		
		M[k++] = eptr;
	}

	return k;
}


int make_graph_eulerian(graph* Gptr)
{
	int rval = 0;
	tree T;

	//Determine set of nodes with odd degree
	node** oddnodes = (node**)malloc(sizeof(node*)*Gptr->n);
	int oddlength = 0;
	int i;
	for(i=0; i<Gptr->n; i++){
		if(Gptr->nodes[i].degree%2){
			oddnodes[i] = &(Gptr->nodes[i]);
			oddlength++;
		}
	}
	//compute the length of the shortest path between all of them with lengths L_ij
	int j;
	double oddcost[oddlength][oddlength];
	for(i=0; i<oddlength; i++){
		for(j=i+1; j<oddlength; j++){
			dijkstra(Gptr,oddnodes[i]->id,oddnodes[j]->id,&T);
			oddcost[i][j] = T.cost;
			oddcost[j][i] = oddcost[i][j];
		}
	}
	//construct a complete graph over the odd degree nodes, with edge length L_ij
	
	//Compute a minimum weight perfect matching in it
	edge**M = (edge**)malloc(sizeof(edge*)*Gptr->n/2);
	min_weight_perfect_matching(Gptr,M);

	//Add the edges of the matching to G


	return rval;
}

int in_list(int id, node*list, int n){
	int rval = 0;
	int k = 0;
	while(k<n && !rval){
		if(list[k].id == id){
			rval = 1;
		}
		k++;
	}

	return rval;
}

int degree_in_tree(int x, tree*Tptr){
	int i;
	int n = Tptr->n;
	int degree = 0;

	for(i=0; i<n-1; i++){
		if(Tptr->edges[i]->i == x){
			degree++;
		}
		else{
			if(Tptr->edges[i]->j == x){
				degree++;
			}
		}
	}

	return degree;
}

int odd_nodes(graph*Gptr, tree*Tptr, node*oddnodes){
	int n = Tptr->n;
	int i;
	int k = 0;

	for(i=0; i<n-1; i++){
		if(degree_in_tree(Tptr->edges[i]->i,Tptr) % 2 != 0){
			if(!in_list(Gptr->nodes[Tptr->edges[i]->i].id, oddnodes,k)){
				oddnodes[k++] = Gptr->nodes[Tptr->edges[i]->i];
			}
		}

		if(degree_in_tree(Tptr->edges[i]->j,Tptr) % 2 != 0){
			if(!in_list(Gptr->nodes[Tptr->edges[i]->j].id, oddnodes,k)){
				oddnodes[k++] = Gptr->nodes[Tptr->edges[i]->j];
			}
		}
	}

	return k;
}

void create_complete_graph(graph*Kn, node*oddnodes, int n, tree*Tptr){
	int i,j;
	int m = n*(n-1)/2;

	Kn->n = n;
	Kn->m = m;

	Kn->nodes = (node*)malloc(sizeof(node)*n);
	node* iptr;
	for(i=0; i<n; i++){
		iptr = &(Kn->nodes[i]);
		iptr->id = oddnodes[i].id;
		iptr->is_in_tree = 0;
		iptr->degree= 0;
		iptr->tmp_degree= 0;
	}

	Kn->edges = (edge*)malloc(sizeof(edge)*m);
	double cost;
	double max_edge_cost = -LARGE_NUMBER;
	edge* eptr;
	int k = 0;
	for(i=0; i<n; i++){
		for(j=i+1; j<n; j++){
			eptr = &(Kn->edges[k]);
			eptr->id = k;
			eptr->i = oddnodes[i].id;
			eptr->j = oddnodes[j].id;
			cost = Tptr[k].cost;
			eptr->cost = cost;
			eptr->is_in_tree = 0;
			eptr->displayed = 1;
			if(max_edge_cost < cost){
				max_edge_cost = cost;
			}
			//We update the degree of the edge's endpoints
			//tail node
			iptr = &(Kn->nodes[i]);
			iptr->degree++;
			//head node
			iptr = &(Kn->nodes[j]);
			iptr->degree++;
			k++;
		}
	}
	Kn->max_edge_cost = max_edge_cost;

	int degree;
	for(i=0; i<n; i++){
		iptr = &(Kn->nodes[i]);
		degree = iptr->degree;
		iptr->edges = (edge**)malloc(sizeof(edge*)*degree);
	}

	int tmp;
	k = 0;
	for(i=0; i<n; i++){
		for(j=i+1; j<n; j++){
			eptr = &(Kn->edges[k]);

			//tail node
			iptr = &(Kn->nodes[i]);
			tmp = iptr->tmp_degree;
			iptr->edges[tmp] = eptr;
			iptr->tmp_degree++;

			//head node
			iptr = &(Kn->nodes[j]);
			tmp = iptr->tmp_degree;
			iptr->edges[tmp] = eptr;
			iptr->tmp_degree++;
			k++;
		}
	}
}

void union_matching_tree(graph*Gptr, tree*T, edge**M, int matching_size){
	int e,i,j;
	int n = T->n;
	int m = (T->n)-1+matching_size;
	Gptr->n = n;
	Gptr->m = m;

	Gptr->nodes = (node*)malloc(sizeof(node)*n);
	node* iptr;
	for(i=0; i<n; i++){
		iptr = &(Gptr->nodes[i]);
		iptr->id = i;
		iptr->is_in_tree = 0;
		iptr->degree = 0;
		iptr->tmp_degree = 0;
	}
	
	Gptr->edges = (edge*)malloc(sizeof(edge)*m);
	edge* eptr;
	double cost;
	double max_edge_cost = -LARGE_NUMBER;
	// Arretes de l'arbre
	for(e=0; e<n-1; e++){
		eptr = &(Gptr->edges[e]);
		eptr->id = e;
		eptr->i = T->edges[e]->i;
		eptr->j = T->edges[e]->j;
		eptr->cost = T->edges[e]->cost;
		cost = eptr->cost;
		if(max_edge_cost < cost){
			max_edge_cost = cost;
		}
		eptr->is_in_tree = T->edges[e]->is_in_tree;
		eptr->displayed = 1;
		i = eptr->i;
		j = eptr->j;
		// Mise à jour du degré des noeuds
		iptr = &(Gptr->nodes[i]);
		iptr->degree++;
		iptr = &(Gptr->nodes[j]);
		iptr->degree++;
	}
	// Arretes du matching
	for(e=n-1; e<n-1+matching_size; e++){
		eptr = &(Gptr->edges[e]);
		eptr->id = e;
		eptr->i = M[e-(n-1)]->i;
		eptr->j = M[e-(n-1)]->j;
		eptr->cost = M[e-(n-1)]->cost;
		cost = eptr->cost;
		if(max_edge_cost < cost){
			max_edge_cost = cost;
		}
		eptr->is_in_tree = M[e-(n-1)]->is_in_tree;
		eptr->displayed = 1;
		i = eptr->i;
		j = eptr->j;
		// Mise à jour du degré des noeuds
		iptr = &(Gptr->nodes[i]);
		iptr->degree++;
		iptr = &(Gptr->nodes[j]);
		iptr->degree++;
	}
	Gptr->max_edge_cost = max_edge_cost;

	int degree;
	for(i=0; i<n; i++){
		iptr = &(Gptr->nodes[i]);
		degree = iptr->degree;
		iptr->edges = (edge**)malloc(sizeof(edge*)*degree);
	}

	int tmp;
	for(e=0; e<m; e++){
		eptr = &(Gptr->edges[e]);
		i = eptr->i;
		j = eptr->j;

		iptr = &(Gptr->nodes[i]);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;

		iptr = &(Gptr->nodes[j]);
		tmp = iptr->tmp_degree;
		iptr->edges[tmp] = eptr;
		iptr->tmp_degree++;
	}
}


int eulerian_prim(graph *Gptr, int origin_node, tree *Tptr)
{
	int rval = 1;

	int m = Gptr->m;
	int n = Gptr->n;
	int e, i, j;

	//We reinitialize the edges and nodes to not be in the tree
	for (e = 0; e < m; e++)
		Gptr->edges[e].is_in_tree = 0;
	for (i = 0; i < n; i++)
		Gptr->nodes[i].is_in_tree = 0;
	//We mark the origin node as in the tree
	Gptr->nodes[origin_node].is_in_tree = 1;

	fprintf(stderr, "Starting PRIM's algorithm from node %d\n", origin_node);

	//Prim's main body
	double mincost;
	double maxcost = Gptr->max_edge_cost;
	edge * eptr;
	node *iptr, *jptr;
	int candidate_e;
	int tree_size = 0;
	for (tree_size = 0; tree_size < n - 1; tree_size++)
	{
		mincost = maxcost;
		candidate_e = -1;
		//We look for the edge (i,j) with minimum weight such that 
		// (i,j) not already in the tree
		// ANDi NOT
		// i and j both in the tree
		// OR i and j both NOT in the tree 

		for (e = 0; e < m; e++)
		{
			//We fetch the pointer to the edge
			eptr = &(Gptr->edges[e]);

			//If it is already in the tree, skip it
			if (eptr->is_in_tree || eptr->displayed == 0)
				continue;

			//We fetch the pointers to the endpoint nodes
			i = get_node_by_id(Gptr,eptr->i);
			iptr = &(Gptr->nodes[i]);
			j = get_node_by_id(Gptr,eptr->j);
			jptr = &(Gptr->nodes[j]);
			//printf("i%d j%d\n",i,j);

			//If i and j are in the tree, skip it
			if ((iptr->is_in_tree) && (jptr->is_in_tree))
				continue;

			//If i and j are both outside the tree skip it
			if (!(iptr->is_in_tree) && !(jptr->is_in_tree))
				continue;

			//Otherwise, we consider it and check its cost
			if (eptr->cost < mincost)
			{
				mincost = eptr->cost;
				candidate_e = eptr->id;;
			}
		}

		//If we did not find a candidate, the graph is not connected 
		if (candidate_e == -1)
		{
			fprintf(stderr, "The graph is disconnected. Aborting Prim\n");
			return 0;
		}

		//Otherwise, we add the least cost edge we found to the tree
		//We fetch the pointer to the edge
		eptr = &(Gptr->edges[candidate_e]);
		//We mark it as in the tree
		eptr->is_in_tree = 1;

		//We fetch the pointers to the endpoint nodes
		i = get_node_by_id(Gptr,eptr->i);
		iptr = &(Gptr->nodes[i]);
		j = get_node_by_id(Gptr,eptr->j);
		jptr = &(Gptr->nodes[j]);
		//And we mark both nodes as marked (one of them already is)
		iptr->is_in_tree = 1;
		jptr->is_in_tree = 1;
	}

	//We fill the tree 
	Tptr->cost = 0;
	int count = 0;
	//We fill the tree structure
	for (e = 0; e < m; e++)
	{
		eptr = &(Gptr->edges[e]);
		//If the edge is not marked, we skip it
		if (!eptr->is_in_tree)
			continue;
		//Otherwise, we put it in the tree
		Tptr->edges[count++] = eptr;
		Tptr->cost += eptr->cost;
	}
	return rval;
}
// cette fonction renvoit 1 si il trouve un arbre couvrant, 0 si le graphe est déconnecté
// cela va nous permettre de savoir si une arête est un pont ou non

int alone(graph*Gptr,int s){
	int b = 1;
	int i = 0;
	node*nptr = &(Gptr->nodes[s]);
	while(b && i<nptr->degree){
		if(nptr->edges[i]->displayed){
			b = 0;
		}
		i++;
	}
	return b;
}

int is_bridge(graph *Gptr, edge *e, int s, int*nbnodes){
	int i,start,end;
	edge*eptr;
	node*iptr;

	// On supprime l'arête du graphe
	for(i=0; i<Gptr->m; i++){
		eptr = &(Gptr->edges[i]);
		if(eptr->id == e->id){
			eptr->displayed = 0;
			break;
		}
	}
	// On remet à 0 les attributs du DFS
	for(i=0; i<Gptr->n; i++){
		iptr = &(Gptr->nodes[i]);
		iptr->father = NULL;
		iptr->visited= 0;
		iptr->beg= -1;
		iptr->end= -1;
	}

	if(s == e->i){
		start = e->i;
		end = e->j;
	}
	else{
		start = e->j;
		end = e->i;
	}
	if(alone(Gptr,start)){
		start = end;
		(*nbnodes)--;
	}
	DFS_explore(Gptr,start);
	int count = 0;
	for(i=0; i<Gptr->n; i++){
		iptr = &(Gptr->nodes[i]);
		if(iptr->visited){
			count++;
		}
	}

	printf("nbnodes%d count%d\n",*nbnodes,count);
	return !(*nbnodes == count);
}

int empty_graph(graph*Gptr){
	int i = 0, rval = 1;
	edge* eptr;
	while(i<Gptr->m && rval){
		eptr = &(Gptr->edges[i]);
		if(eptr->displayed){
			rval = 0;
		}
		i++;
	}

	return rval;
}

void fleury(graph*Gptr, tree*Tptr, int s){
	int i,k,b, size = 0;
	int nbnodes = Gptr->n;
	node*nptr;
	edge*eptr;
	i = s;

	while(!empty_graph(Gptr)){
		b = 1;
		k = 0;
		nptr = &(Gptr->nodes[i]);
		while(b && k<nptr->degree){
			eptr = nptr->edges[k++];
			//printf("Edge #%d ((%d,%d) [%lf])\n",eptr->id,eptr->i,eptr->j,eptr->cost);
			if(!eptr->displayed){
				continue;
			}
			
			if(!is_bridge(Gptr,eptr,i,&nbnodes)){
				b = 0;
			}
			else{
				eptr->displayed = 1;
			}
		}
		if(b){
			printf("bug\n");
			break;
		}

		//printf("Edge #%d [%d %d]\n",eptr->id,eptr->i,eptr->j);
		Tptr->edges[size++] = eptr;
		Tptr->cost += eptr->cost;
		eptr->displayed = 0;
		if(i == eptr->i){
			i = eptr->j;
		}
		else{
			i = eptr->i;
		}
		printf("i %d\n",i);
	}
	for(i=0; i<Gptr->m; i++){
		eptr = Tptr->edges[i];
		eptr->displayed = 1;
	}
}

void hamiltonian_cycle(graph * G,tree * cycleE,tree * cycleH){


	int n = G->n;
	for(int i = 0 ;  i < n ; i++)
	{
		node * iptr = &(G->nodes[i]);
		//The node is initially not in the tree
		iptr->is_in_tree = 0;

	}
	
	
	int startnode = cycleE->edges[0]->i;
	int currnode = startnode;
	G->nodes[cycleE->edges[0]->i].is_in_tree = 1;

	int m = cycleE->n;
	int k = 0;
	int cost = 0;
	

	for(int e = 0 ; e < m - 2;e++){
		//si le noeud courant est la composante i de l'arete
		if(currnode == cycleE->edges[e]->i)
		{
			//si la composante j n'est pas dans le chemin hamiltonien
			if(!G->nodes[cycleE->edges[e]->j].is_in_tree)
			{
				//si il n'y a pas de noeud intermediaire deja dans le cycle
				if(startnode == currnode)
				{
					cycleH->edges[k] = cycleE->edges[e];
					cost = cycleE->edges[e]->cost;
				}
				//s'il y a eu un noeud intermediaire on crée une nouvelle arrete entre le depart et l'arrivée
				else
				{	
					cost += cycleE->edges[e]->cost;
					edge * ed = malloc(sizeof(edge));
					ed->cost = cost;
					ed->i = startnode;
					ed->j = cycleE->edges[e]->j;
					ed->id = cycleE->edges[e]->id;
					cycleH->edges[k] = ed;				
				}

				G->nodes[cycleE->edges[e]->j].is_in_tree = 1;
				currnode = cycleE->edges[e]->j;
				cycleH->cost += cost;
				startnode = currnode;
				cost = 0;
				k++;
			}
			else
			{
				//si le noeud courant est deja dans le cycle, on stocke le cout et on avance
				cost += cycleE->edges[e]->cost;				
				currnode = cycleE->edges[e]->j;				
			}
		}
		else{
			//meme processus mais si le noeud courant est la composante j de l'arete
			if(currnode == cycleE->edges[e]->j)
			{
				if(!G->nodes[cycleE->edges[e]->i].is_in_tree)
				{	
					if(startnode == currnode)
					{
						cycleH->edges[k] = cycleE->edges[e];	
						cost = cycleE->edges[e]->cost;
					}
					else
					{
						cost += cycleE->edges[e]->cost;
						edge * ed = malloc(sizeof(edge));
						ed->cost = cost;
						ed->i = startnode;
						ed->j = cycleE->edges[e]->i;
						ed->id = cycleE->edges[e]->id;

						cycleH->edges[k] = ed;						
					}

					G->nodes[cycleE->edges[e]->i].is_in_tree = 1;
					currnode = cycleE->edges[e]->i;
					cycleH->cost += cost;
					startnode = currnode;
					cost = 0;
					k++;
				}
				else
				{
					cost += cycleE->edges[e]->cost;
					currnode = cycleE->edges[e]->i;
				}
			}
			else
			{
				//si le noeud courant est ni dans i ou j de l'arete en cours
				printf("Erreur dans la création du cycle Hamiltonien\n");
				break;
			}
		}
	}
	
	//on complete le cycle
	cost += cycleE->edges[m-2]->cost;

	edge * ed = malloc(sizeof(edge));
	ed->cost = cost;
	ed->i = startnode;
	ed->j = G->edges[0].i;
	ed->id = cycleE->edges[m-2]->id;

	cycleH->edges[cycleH->n - 2] = ed;
	cycleH->cost += cost;
}

void write_hamiltonian_cycle(tree*Tptr, int n, char*filename){
	FILE* f = fopen(filename,"w");
	int i;
	edge*eptr;

	if(f==NULL){
		printf("Erreur d'ouverture\n");
	}
	else{
		fprintf(f,"%d %f\n",n,Tptr->cost);
		for(i=0; i<(Tptr->n)-1; i++){
			eptr = Tptr->edges[i];
			fprintf(f,"%d %d %f\n",eptr->i,eptr->j,eptr->cost);
		}
	}
	
}