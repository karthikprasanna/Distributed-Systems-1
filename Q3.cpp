#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

bool is_alive(int **grid, int* top, int* bottom, int N, int M, int i, int j)
{
    if (i >= 0 && i < N)
    {
        if (j >= 0 && j < M)
        {
            return grid[i][j];
        }
        
    }
    else if (i == -1 && top)
    {
        if (j >= 0 && j < M)
        {
            return top[j];
        }
    }
    else if (i == N && bottom)
    {
        if (j >= 0 && j < M)
        {
            return bottom[j];
        }
    }
    
    return 0;
}

int num_neighbours(int **grid, int* top, int* bottom, int N, int M, int i, int j)
{
    return is_alive(grid, top, bottom, N, M, i-1, j-1)+is_alive(grid, top, bottom, N, M, i-1, j)+is_alive(grid, top, bottom, N, M, i-1, j+1)
        +is_alive(grid, top, bottom, N, M, i+1, j-1)+is_alive(grid, top, bottom, N, M, i+1, j)+is_alive(grid, top, bottom, N, M, i+1, j+1)
            +is_alive(grid, top, bottom, N, M, i, j-1)+is_alive(grid, top, bottom, N, M, i, j+1);
    
}


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int recieve_configurations;
    int N, M, T;
    int **grid, **gathered;
   
    int *flatten_send, *flatten_recieve, *buffer;
    int counts_send[size], displacements[size], count_recieve;

    if (my_rank == 0)
    {
        cin>>N>>M>>T;
        flatten_send = (int*)malloc(N * M * sizeof(int));
        grid = (int**)malloc(N * sizeof(int*));
        for (int i = 0; i < N; ++i)
        {
            grid[i] = &flatten_send[M * i];
            for (int j = 0; j < M; j++)
            {
                cin>>grid[i][j];
            }
            
        }

        int block_size = N / size;
        

        if (block_size * size != N)
        {
            block_size++;
        }

        for (int i = 0; i < size; i++) {
            counts_send[i] = block_size * M;
            displacements[i] = i * block_size * M;

            if (i >= N)
            {
                counts_send[i] = 0;
                displacements[i] = 0;
            }
        }

        if (N % size != 0)
        {
            counts_send[size - 1] = N * M - block_size * (size - 1) * M;
        }
        
        flatten_recieve = (int*)malloc(N * M * sizeof(int));
        gathered = (int**)malloc(N * sizeof(int*));
        for (int i = 0; i < N; ++i)
        {
            gathered[i] = flatten_recieve + i * M;
        }

    }
    MPI_Bcast(&counts_send, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
 
    count_recieve = counts_send[my_rank];
    buffer = (int*)malloc(count_recieve * sizeof(int));
    if (my_rank == 0)
    MPI_Scatterv(grid[0], counts_send, displacements, MPI_INT, buffer, count_recieve, MPI_INT, 0, MPI_COMM_WORLD);
    else
    MPI_Scatterv(NULL, NULL, NULL, MPI_INT, buffer, count_recieve, MPI_INT, 0, MPI_COMM_WORLD);
    
    int n = count_recieve / M;
    
    int** block = (int**)malloc(n * sizeof(int*));
    int** simulated = (int**)malloc(n * sizeof(int*));
    int *flatten_simulated = (int*)malloc(count_recieve * sizeof(int));
    int *flatten_block = (int*)malloc(count_recieve * sizeof(int));
    int *bottom = NULL;
    int *top = NULL;

    if (my_rank > 0)
    {
        top = (int*)malloc(M * sizeof(int));
    }
    
    if (my_rank < size - 1)
    {
        bottom = (int*)malloc(M * sizeof(int));
    }
    

    for (int i = 0; i < n; i++)
    {
        block[i] = flatten_block + i * M;
        simulated[i] = flatten_simulated + i * M;
        for (int j = 0; j < M; j++)
        {
            block[i][j] = buffer[i * M + j];
        }
        
    }

    if (my_rank % 2 == 0)
    {
        if (my_rank > 0)
        {
            MPI_Send(block[0], M, MPI_INT, my_rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(top, M, MPI_INT, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (my_rank < size - 1)
        {
            MPI_Send(block[n - 1], M, MPI_INT, my_rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(bottom, M, MPI_INT, my_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }
    }
    else
    {
        if (my_rank > 0)
        {
            MPI_Recv(top, M, MPI_INT, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(block[0], M, MPI_INT, my_rank - 1, 0, MPI_COMM_WORLD);
        }

        if (my_rank < size - 1)
        {
            MPI_Recv(bottom, M, MPI_INT, my_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(block[n - 1], M, MPI_INT, my_rank + 1, 0, MPI_COMM_WORLD);

        }
    }

    // cout<<endl;
    // cout<<top<<bottom<<endl;
    int count;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < M; j++)
        {
            count = num_neighbours(block, top, bottom, n, M, i, j);
            // cout<<count<<' ';
            if (is_alive(block, top, bottom, n, M, i, j))
            {
                if (count < 2 || count > 3)
                {
                    simulated[i][j] = 0;
                }
                else
                {
                    simulated[i][j] = 1;
                }
            }
            else
            {
                if (count == 3)
                {
                    simulated[i][j] = 1;
                }
                else
                {
                    simulated[i][j] = 0;
                }
                
            }
        }
    // cout<<endl;    
    }

    if (my_rank == 0)
    {
        MPI_Gatherv(simulated[0], count_recieve, MPI_INT, gathered[0], counts_send, displacements, MPI_INT, 0, MPI_COMM_WORLD);
        cout<<endl;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                cout<<gathered[i][j]<<' ';
            }
            cout<<endl;
            
        }
        free(flatten_send);
        free(flatten_recieve);
        free(flatten_block);
        free(flatten_simulated);
        free(top);
        free(bottom);
        free(block);
        free(simulated);
        free(gathered);
        free(grid);
    }
    else
    {
        MPI_Gatherv(simulated[0], count_recieve, MPI_INT,NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    free(buffer);

    MPI_Finalize();
    return EXIT_SUCCESS;    
}

/*
10 10 7
0 0 1 0 0 0 0 0 0 0
1 0 1 0 0 0 0 0 0 0
0 1 1 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0

4 4 1
0 0 0 0
0 0 0 0
0 0 1 0
0 0 0 0


*/

