#include "Functions.h"

void readTP3Instance(graph* Gptr, char* instanceFileName)
{
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
}


void display_graph(graph*Gptr)
{
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
		fprintf(stderr,"\tNode %d : visited = %d, father = %d, dates = [%d,%d]\n"
                ,i, iptr->visited,iptr->father,iptr->dates[0],iptr->dates[1]);
        fprintf(stderr,"\tNode %d has %d neighbors:\n",i, iptr->nb_neighbors);
		for( e = 0 ; e < iptr->nb_neighbors ; e++)
		{
			eptr = iptr->neighbors[e];
			fprintf(stderr,"\t\t(%d,%d)\t[w=%lf]\n",eptr->i,eptr->j,eptr->cost);
		}
	}	
}


void free_graph(graph*Gptr)
{
	int n,i;
	n = Gptr->n;
	for( i = 0 ; i < n ; i++)
		free(Gptr->nodes[i].neighbors);


	free(Gptr->edges);
	free(Gptr->nodes);
}

void explore(node*Node, graph*Gptr, int*t)
{
    int i;
    int id = Node->id;
    Node->visited = 1;
    Node->dates[0] = ++*t;

    for(i=0; i<Node->nb_neighbors; i++)
    {
        if(Gptr->nodes[Node->neighbors[i]->j].id != id)
        {
            if(Gptr->nodes[Node->neighbors[i]->j].visited == 0)
            {
                Gptr->nodes[Node->neighbors[i]->j].visited = 1;
                Gptr->nodes[Node->neighbors[i]->j].father = id;
                explore(&(Gptr->nodes[Node->neighbors[i]->j]),Gptr,t);
            }
        }
        else
        {
            if(Gptr->nodes[Node->neighbors[i]->i].visited == 0)
            {
                Gptr->nodes[Node->neighbors[i]->i].visited = 1;
                Gptr->nodes[Node->neighbors[i]->i].father = id;
                explore(&(Gptr->nodes[Node->neighbors[i]->i]),Gptr,t);
            }
        }
    }

    Node->dates[1] = ++*t;
}

void DFS(graph*Gptr)
{
    int n = Gptr->n;
    int i,j;
    int t = 0;
    node* jptr;
    for(j=0; j<n; j++)
    {
        jptr = &(Gptr->nodes[j]);
        jptr->visited = 0;
        jptr->father = -1;
        jptr->dates[0] = 0;
        jptr->dates[1] = 0;
    }

    for(i=0; i<n; i++)
    {
        if(Gptr->nodes[i].visited == 0)
        {
            explore(&(Gptr->nodes[i]),Gptr,&t);
        }
    }
}