# Traveling Salesman Problem (TSP) - Christofides Algorithm

## Overview of the Traveling Salesman Problem (TSP)

The Traveling Salesman has a challenge: visiting multiple cities with different distances between each city as quickly as possible. The goal is to find the shortest circuit, passing through each city exactly once.

## Overview of the Christofides Algorithm

The Christofides algorithm is a method designed to solve the Traveling Salesman Problem. Invented by Nico Christofides in 1976, it guarantees finding a path of less than 1.5 times the optimal path. The algorithm involves finding the minimum spanning tree (using Prim's algorithm), creating a complete graph with odd-degree vertices in the spanning tree, finding shortest paths between odd-degree nodes (using Dijkstra's algorithm), obtaining a minimum-weight perfect matching in the complete graph of odd-degree vertices, and finally, combining the minimum spanning tree and the perfect matching to form an Eulerian cycle and then a Hamiltonian cycle.

## Prim's Algorithm

Prim's algorithm is a greedy algorithm that returns a minimum-weight spanning tree. The code takes a pointer to a graph object containing a list of edges and nodes, the index of the starting node (origin node), and a pointer to a tree object that will be filled with the minimum-weight spanning tree found.

The algorithm starts by resetting all edges and nodes to be not in the tree. Then, the starting node is marked as in the tree, and Prim's algorithm begins. It repeatedly adds the cheapest edge connecting a tree node to an out-of-tree node until all nodes are in the tree. In each iteration, the code traverses all unmarked edges to find the cheapest edge that satisfies the mentioned conditions. If such an edge is found, it is added to the tree by marking it as in the tree and updating the status of its nodes.

At the end of the algorithm, the minimum-weight tree is filled in the tree object using the marked edges in the original graph.

## Creating the Complete Graph of Odd-Degree Vertices

To choose the edges incident to odd-degree vertices to add and create an Eulerian graph for a Hamiltonian cycle, we create a complete graph of odd-degree vertices. The edge values between odd-degree vertices correspond to the shortest paths (obtained with Dijkstra's algorithm) in the original graph. This complete graph is then used to find a minimum-weight perfect matching.

### Dijkstra's Algorithm

Dijkstra's algorithm is a shortest path algorithm used to find the path between two nodes in a weighted graph. It initializes all nodes as not connected to the shortest path tree, giving them an infinite distance to the origin. It visits the node with the smallest distance to the origin and marks it as visited by setting its distance to the origin as the value of the edge connecting them. For each unvisited neighbor of the current node, it calculates the distance passing through the current node. If this distance is shorter than the currently stored distance for the neighbor, it updates the stored distance. These steps are repeated until the desired node is reached or all nodes are visited. If waiting to visit all nodes, it results in a consistent Dijkstra tree where all distances from the origin are minimal, but distances between non-origin nodes are not necessarily minimal.

### Complete Graph Creation Function
The function for creating the complete graph takes as input the array of Dijkstra trees, where each tree gives the distance between two defined points, an allocated complete graph Kn, representing the complete graph of odd-degree vertices, and the array of odd-degree vertices needed for instantiation. It first allocates memory for the nodes and edges of the complete graph. Then, for each pair of nodes, it sets the cost of the corresponding Dijkstra tree as the edge value. The function also updates the degree of each node and stores the maximum edge cost attribute.

Therefore, the odd-degree vertices remain the same, with edge values being the minimal distance between pairs of points.

## Minimum Weight Perfect Matching in the Complete Graph of Odd-Degree Vertices

We perform a minimum weight perfect matching in the complete graph of odd-degree vertices to find the edges to add for creating an Eulerian tour. To achieve this, we use the CPLEX library, a solver for mathematical optimization problems.

The min_weight_perfect_matching function takes as input a pointer to a graph structure Gptr and a pointer to an edge structure M, which will store the matching edges. The function returns an integer value, the optimal value of the objective function.

First, the function initializes a problem structure IPproblem to store the data required by CPLEX. Then, it creates a CPLEX environment using the CPXopenCPLEX function. If the environment cannot be created, the program exits.

Next, the function creates a Mixed Integer Programming (MIP) problem using the CPXcreateprob function. It sets certain parameters using CPXsetintparam and CPXsetdblparam.

The function then generates variables and constraints necessary for the perfect matching problem using the input graph data. It creates a binary variable xe for each edge in the graph and generates constraints to ensure that each node in the graph is paired with exactly one other node. The objective function is defined to minimize the sum of weights of edges included in the matching.

Finally, the function calls CPXmipopt to solve the problem and retrieves the optimal solution using CPXgetx. The matched edges are stored in the M array, and the optimal value of the objective function is returned.

## Union of the Minimum Weight Perfect Matching and the Minimum Weight Spanning Tree

The union of the minimum weight perfect matching and the minimum weight spanning tree is achieved using the union_matching_tree function. This union creates a new graph that combines the edges of the minimum spanning tree and the minimum weight perfect matching, creating the necessary conditions to construct an Eulerian cycle.

The function takes as input a pointer to a graph structure Gptr, a tree T, a list of edges M, and the size of the M list. It creates a new graph Gptr that is the union of the tree T and the edge list M.

This function starts by initializing the nodes of Gptr by assigning them a unique identifier, marking them as not in the tree T, initializing their degree and a temporary degree to 0. Then, it creates the edges of Gptr by copying the edges of T and M, updating their cost and marking them as in the tree T or not, and updating the degrees of adjacent nodes for each edge.

The code allocates memory for edge lists adjacent to each node of Gptr and updates these lists for each node using the edges of the graph.

## Eulerian Cycle

To obtain an Eulerian cycle on the union of the minimum weight perfect matching and the minimum weight spanning tree, we use Fleury's algorithm. It is essential to determine if an edge is a bridge for the algorithm's proper functioning. For this purpose, a displayed integer was added to the edge structure to indicate whether the edge is still part of the graph or not (1 if the edge still exists, 0 otherwise).

We use a function is_bridge that takes a pointer to a graph structure Gptr, a pointer to an edge structure e, a vertex defined by an integer s, and a pointer to the number of accessible nodes represented by nbnodes. First, the function removes the edge from the graph and chooses one of the two vertices of the removed edge to be the starting vertex. If one of the two vertices is excluded from the graph, the other vertex is chosen. Then, a DFS_explore is performed from the selected starting vertex, counting the number of visited vertices after the exploration. If the number of visited vertices equals the number of accessible nodes in the graph, then the edge is not a bridge. Otherwise, the algorithm did not explore the graph entirely, indicating multiple connected components after the edge is removed, making the edge a bridge.

The Eulerian path is stored in a tree structure. The fleury function takes as input a pointer to a graph structure Gptr, a pointer to a tree structure T, and a starting vertex s. While the graph is not empty (empty if all edges are not displayed), an incident edge to the chosen vertex is selected. If this edge is not a bridge, it is added to the tree, the incident node becomes the new starting vertex, and the current node is advanced without changing the starting node. If the current node is already in the tree, the edge cost is simply added to the local cost, and the current node is advanced without changing the starting node.

If the current node is neither the starting nor the incident node, it indicates an error in the creation of the Hamiltonian cycle, or the input Eulerian cycle was incorrectly instantiated. The algorithm then displays an error message and exits the loop.

After traversing all edges except the last one of the Eulerian cycle, the algorithm adds the last edge of the Eulerian cycle to the total cost of the tree. It then creates an edge between the last visited node and the first node of the graph (storing the total cost of the edge) and adds this edge to the Hamiltonian cycle tree.

Finally, the algorithm returns the Hamiltonian cycle (the weight is always the same as that of the Eulerian cycle).

## Hamiltonian Cycle

From the Eulerian cycle, we create a Hamiltonian cycle that passes through all nodes exactly once. The algorithm involves traversing the Eulerian cycle and removing duplicates so that each node is visited only once (though it saves the cost). The algorithm takes as input the original graph G, and two trees (cycleE and cycleH) used to store the obtained Eulerian cycle (a cycle that passes through each edge of the graph exactly once) and the Hamiltonian cycle. The algorithm relies on the fact that the edges of the Eulerian cycle are sorted in the order of traversal to create a sorted Hamiltonian cycle. (Note that the starting node changes during the algorithm to take the value of the last visited node, the new starting point for the next edge of the graph).

The algorithm starts by initializing all nodes of the graph as not in the Hamiltonian cycle (isintree = 0). It then selects the first node of the Eulerian cycle (cycleE) and adds it to the cycle tree (isintree = 1).

Next, for each edge of the Eulerian cycle (except the last one), the algorithm checks if the current node is the endpoint i or j of the edge. If it is i and j is not already in the Hamiltonian cycle, the algorithm adds an edge between the starting node i and j (storing the total cost of the edge) and marks j as part of the Hamiltonian cycle. The current and starting nodes then become j. If the current node is already in the Hamiltonian cycle, the algorithm simply adds the edge cost to the local cost and advances the current node without changing the starting node.

If the current node is j and i is not already in the Hamiltonian cycle, the algorithm adds an edge between the starting node and i (storing the total cost of the edge) and marks i as part of the Hamiltonian cycle. The current and starting nodes then become i. If the current node is already in the Hamiltonian cycle, the algorithm simply adds the edge cost to the local cost and advances the current node without changing the starting node.

If the current node is neither i nor j, it indicates an error in the creation of the Hamiltonian cycle or that the input Eulerian cycle was incorrectly instantiated. The algorithm then displays an error message and exits the loop.

After traversing all edges except the last one of the Eulerian cycle, the algorithm adds the last edge of the Eulerian cycle to the total cost of the Hamiltonian cycle tree. It then creates an edge between the last visited node and the first node of the graph (storing the total cost of the edge) and adds this edge to the Hamiltonian cycle tree.

Finally, the algorithm returns the Hamiltonian cycle (the weight is always the same as that of the Eulerian cycle).

## Conclusion

In conclusion, the Traveling Salesman no longer has a problem. Now, there's an algorithm that, given a graph of cities and distances as input, provides a satisfying tour of less than 1.5 times the optimal solution, regardless of the complexity of the input graph. The realization of this comprehensive algorithm has been interesting, as it required the use of many concepts and algorithms from graph theory. Moreover, its implementation concretizes graph theory, as the Christofides algorithm has been used for various applications, such as vehicle route planning, integrated circuit design, robot path planning, and task scheduling. Finally, it's worth noting that the Christofides algorithm has also served as the basis for many other heuristic methods for the TSP, making it a valuable and versatile tool. In this sense, this work has been exciting and useful to accomplish.
