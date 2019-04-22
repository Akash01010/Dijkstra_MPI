# Dijkstra_MPI
Parallelising Dijkstra algorithm using MPI programming


I have to parallelize Dijkstra algorithm. Dijkstra is a greedy algorithm and main overhead in the algorithm is while it calculates the minimum weight vertex to be explored next. I have parallelized this particular step of calculating minimum of weights of the vertices unexplored.

To Run the program:
Step 1: Type mpicc -o dijkstra dijkstra.c -lm
Step 2: mpiexec -n 1 ./dijkstra
(It will ask to input no of vertices and adjacency matrix and starting node. One example for testing is provided in Test_example.txt. Graph is attached for visual representation of testing example).

It will show the Shortest path from starting node to all other vertices. 
