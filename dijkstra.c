#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <mpi.h>
#ifndef MAX
#define MAX 10
#endif
int min(int a,int b) {
	return (((a)<(b))?(a):(b));
}
int max(int a,int b) {
	return (((a)>(b))?(a):(b));
}

int getMinDist(int n,int distance[],int visited[]);

void dijkstra(int **graph,int n,int startnode);

int main(int argc, char* argv[])
{
	int my_rank,n,u,**graph,size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	double start,end;
	if(my_rank==0){
		//Taking input of the graph
		graph=(int**)malloc(MAX*sizeof(int*));
		for(int i=0;i<MAX;i++){
			graph[i]=(int*)malloc(MAX*sizeof(int));
		}
		printf("Enter no. of vertices:");
		scanf("%d",&n);
		printf("\nEnter the adjacency matrix:\n");
		
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				scanf("%d",&graph[i][j]);
			}
		}
		
		printf("\nEnter the starting node:");
		scanf("%d",&u);
		//Startint clock time
		start=MPI_Wtime();
	}
	//This is to set the value of number of nodes as needed by each process for calling dijkstra in next line
	//Starting node need not be broadcasted
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//Each process 
	dijkstra(graph,n,u);
	if(my_rank==0){
		//ending clock time
		end=MPI_Wtime();
		printf("\nThe algorithm took %f seconds to run with %d processes\n", end-start,size);
	}
	MPI_Finalize();
	return 0;
}

int getMinDist(int n,int distance[],int visited[]){
	int mindistance=INT_MAX,nextnode=INT_MAX;
	for(int i=0;i<n;i++){
		if(distance[i]<mindistance&&!visited[i]){
			mindistance=distance[i];
			nextnode=i;
		}
	}
	return nextnode;
}

void dijkstra(int **graph,int n,int startnode)
{
	int my_rank,num_process;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_process);
	int **cost,*distance,*predecessor,*visited;
 	if(my_rank==0){
 		//I have dynamically allocated memory as if it would have been alloted outside then every process would have alloted that memory space and for no use.
 		//This will initialize cost matrix, distance matrix and predecessor matrix
 		cost=(int**)malloc(MAX*sizeof(int*));
		for(int i=0;i<MAX;i++){
			cost[i]=(int*)malloc(MAX*sizeof(int));
		}
 		distance=(int*)malloc(MAX*sizeof(int));
 		predecessor=(int*)malloc(MAX*sizeof(int));
 		visited=(int*)malloc(MAX*sizeof(int));
		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++)
				if(graph[i][j]==0)
					cost[i][j]=INT_MAX;
				else
					cost[i][j]=graph[i][j];
		for(int i=0;i<n;i++)
		{
			distance[i]=cost[startnode][i];
			predecessor[i]=startnode;
			visited[i]=0;
		}
		
		distance[startnode]=0;
		visited[startnode]=1;
	}
	int count=1;
	while(count<n-1)
	{
		int *batch_size,batches=ceil((float)n/num_process),*displs;
		if(my_rank==0){
			//This will be helpful in case the number of nodes is not multiple of number of processes
			batch_size=(int*)malloc(num_process*sizeof(int));
			displs=(int*)malloc(num_process*sizeof(int));
			int count=0;
			for(int i=0;i<num_process;i++){
				batch_size[i]=min((i+1)*batches,n)-min(i*batches,n);
				displs[i]=(count);
				count+=batch_size[i];
			}
		}
		//This will tell each processor how many elements  from distance and visited array to expect.
		MPI_Scatter(batch_size,1,MPI_INT,&batches,1,MPI_INT,0,MPI_COMM_WORLD);
		int dist[batches],vis[batches];
		//This will send parts of distance and visited arrays to different processors
		MPI_Scatterv(distance,batch_size,displs,MPI_INT,dist,batches,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Scatterv(visited,batch_size,displs,MPI_INT,vis,batches,MPI_INT,0,MPI_COMM_WORLD);
		//Here processors will calculated the minimum distance node to be explored next according to the subarray they have
		int mindistanceNode=getMinDist(batches,dist,vis);

		//This will set the min distance node value encountered by each process
		int mindist[2],mindistance[2],proc_rank;
		if(mindistanceNode!=INT_MAX){
			mindist[0]=dist[mindistanceNode],mindist[1]=my_rank;
		}else{
			mindist[0]=INT_MAX,mindist[1]=my_rank;
		}
		//This will store the minimum value and the processor sending the minimum value in the mindistance array.
		MPI_Reduce(mindist,mindistance,1,MPI_2INT,MPI_MINLOC,0,MPI_COMM_WORLD);
		if(my_rank==0){
			proc_rank=mindistance[1];
		}

		//This is to get the vertex which need to be explored next by each processor
		if(mindistanceNode!=INT_MAX){
			mindist[0]=dist[mindistanceNode],mindist[1]=mindistanceNode;
		}else{
			mindist[0]=INT_MAX,mindist[1]=mindistanceNode;
		}
		//This will help process 0 get which vertex need to be explored next
		MPI_Reduce(mindist,mindistance,1,MPI_2INT,MPI_MINLOC,0,MPI_COMM_WORLD);
		int nextnode;
		if(my_rank==0){
			mindistance[1]=mindistance[1]+(displs[proc_rank]);
			nextnode=mindistance[1];
			visited[nextnode]=1;
			//This will update values and see if going through the newly explored node there is shorter path to any other unexplored node.
			for(int i=0;i<n;i++){
				if(!visited[i]){
					if(cost[nextnode][i]!=INT_MAX && mindistance[0]+cost[nextnode][i]<distance[i]){
						distance[i]=mindistance[0]+cost[nextnode][i];
						predecessor[i]=nextnode;
					}
				}
			}
		}
		count++;
	}
	//This will print the results
	if(my_rank==0){
		for(int i=0;i<n;i++)
			if(i!=startnode)
			{
				printf("\nDistance of node%d=%d",i,distance[i]);
				printf("\nPath=%d",i);
				
				int j=i;
				do
				{
					j=predecessor[j];
					printf("<-%d",j);
				}while(j!=startnode);
		}
	}
}
