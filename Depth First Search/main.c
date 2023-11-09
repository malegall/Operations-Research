#include "Functions.h"

int main()
{
	graph G;
    readTP3Instance(&G,"TP3instance.csv");
    DFS(&G);
    display_graph(&G);
}
