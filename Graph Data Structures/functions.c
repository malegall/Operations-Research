#include "functions.h"

int readTP1Instance(graph* Gptr, char* instanceFileName)
{
	int rval = 0;

	//We open the file
	FILE* fin = fopen(instanceFileName,"r");

	//We read the number of nodes and vertices
	int n,m;
	fscanf(fin,"%d,%d\n", &n,&m);
	Gptr->n = n;
	Gptr->m = m;
	//fprintf(stderr,"n = %d\t m = %d\n",Gptr->n,Gptr->m);

	//We malloc, read and fill the edges in their structure
	Gptr->edges = (edge*)malloc(sizeof(edge)*m);
	int e,i,j;
	double cost;
	edge* eptr;
	double max_edge_cost = -LARGE_NUMBER;
	for( e = 0 ;  e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		fscanf(fin,"%d,%d,%lf\n", &i,&j,&cost);
		eptr->id = e;
		eptr->i = i;
		eptr->j = j;
		eptr->cost = cost;
		//We update the largest edge cost
		if(max_edge_cost < cost)
			max_edge_cost = cost;

		//fprintf(stderr,"%d,%d,%d,%lf\n",eptr->ideptr->i,eptr->j,eptr->cost);
	}
	//We store the largest edge cost in the graph structure
	Gptr->max_edge_cost = max_edge_cost;

	//We malloc and fill the nodes in their structure
	Gptr->nodes= (node*)malloc(sizeof(node)*n);
	node* iptr;
	for( i = 0 ;  i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		iptr->id = i;
		iptr->nb_neighbors = 0;
		iptr->neighbors = NULL;
		iptr->cur = 0;

	}


	//We loop over the edges to calculate the number of neighbors
	for( e = 0 ; e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		i = eptr->i;
		iptr = &(Gptr->nodes[i]);
		iptr->nb_neighbors++;
	}

	//We now allocate the neighbors arrays
	for( i = 0 ; i < n ; i++)
	{
		iptr = &(Gptr->nodes[i]);
		iptr->neighbors = (edge**)malloc(sizeof(edge*)*iptr->nb_neighbors);
	}


	//We filll the neighbors arrays:
	//We loop over the edges
	for( e = 0 ; e < m ; e++)
	{
		eptr = &(Gptr->edges[e]);
		i = eptr->i;
		iptr = &(Gptr->nodes[i]);
		iptr->neighbors[iptr->cur++] = eptr;
	}

	//We close the file
	fclose(fin);

	return rval;
}


int display_graph(graph*Gptr)
{
	int rval = 0; 
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
		fprintf(stderr,"\tNode %d has %d neighbors:\n",i, iptr->nb_neighbors);

		for( e = 0 ; e < iptr->nb_neighbors ; e++)
		{
			eptr = iptr->neighbors[e];
			fprintf(stderr,"\t\t(%d,%d)\t[w=%lf]\n",eptr->i,eptr->j,eptr->cost);
		}
	}	

	return rval;
}



int free_graph(graph*Gptr)
{
	int rval;

	int n,i;
	n = Gptr->n;
	for( i = 0 ; i < n ; i++)
		free(Gptr->nodes[i].neighbors);


	free(Gptr->edges);
	free(Gptr->nodes);

	return rval;
}

void init_graph(graph*Gptr, int root_node)
{
	
}

int dijkstra(graph*Gptr,int root_node)
{
	int rval = 0;



	return rval;
}

int giseven(graph*Gptr)
{
	int b = 1;
	int nb = Gptr->n;
	int i=0;

	while(b==1 && i<nb)
	{
		if(((Gptr->nodes[i]).nb_neighbors)%2)
		{
			b = 0;
		}
		i++;
	}

	return b;
}