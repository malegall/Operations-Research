#include "functions.h"

int main()
{
	int rval = 0;
	int root = 0;
	int i, n;

	// Structures
	graph G;
	tree T;

	readTP1Instance(&G, "TP0instance.csv");
	display_graph(&G);

	dijkstra(&G, root);

	printf("Pair : %d\n", giseven(&G));
	
	// Freeings
	free_graph(&G);

	return rval;
}
