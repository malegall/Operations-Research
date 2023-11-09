#include "TP1Functions.h"

int main()
{
	//Structures
	graph G;
    readTP1Instance(&G,"instance.csv");
    display_graph(&G);

    tree*T = (tree*)malloc(sizeof(tree));
    T->edges = (edge**)malloc(sizeof(edge*)*((G.n)-1));
    T->n = (G.n)-1;
    
    Prim(T,&G);

    display_tree(T);

	return 0;
}
