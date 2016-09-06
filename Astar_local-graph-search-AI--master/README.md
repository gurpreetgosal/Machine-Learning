# local-graph-search-AI-
Implementation of local graph search algorithm for minimizing the total cost of ordering of vertices of a graph w.r.t all vertices.

We are given vertices V1, …,Vn and possible parent sets for each vertex. Each parent set has an associated cost. Let O be an ordering (a permutation) of the vertices. We say that a parent set of a vertex Vi is consistent with an ordering O if all of the parents come before the vertex in the ordering. Let mc(Vi, O) be the minimum cost of the parent sets of vertex Vi that are consistent with ordering O. The task is to find an ordering O that minimizes the total cost: mc(V1, O) + … + mc(Vn, O).

This problem is NP-Complete and so presumably intractable. Thus, it makes sense to consider using local search algorithm to solve instances of the problem.

Two data files are provided and the code parses the graph data given in the files and implements local search. To work on a particular file modify the file name accordingly in .py file.

The data file for this example would look like:

5

1 4

{},153.466

{3},96.093

{4},97.913

{5},99.835

2 4

{},141.023

{3},122.446
{4},121.576

{5},123.398

3 6

{},169.482

{1},112.109

{2},150.906

{1,2},107.516

{4},51.681

{5},41.775

The first line of the file is the number of vertices. The next lines of the file consist of the variable name
and the number of parent sets, followed by each parent set and its cost on a separate line.
