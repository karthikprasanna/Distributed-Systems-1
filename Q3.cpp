#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

bool is_alive(int **grid, int *top, int *bottom, int N, int M, int i, int j)
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

int num_neighbours(int **grid, int *top, int *bottom, int N, int M, int i, int j)
{
    return is_alive(grid, top, bottom, N, M, i - 1, j - 1) + is_alive(grid, top, bottom, N, M, i - 1, j) + is_alive(grid, top, bottom, N, M, i - 1, j + 1) + is_alive(grid, top, bottom, N, M, i + 1, j - 1) + is_alive(grid, top, bottom, N, M, i + 1, j) + is_alive(grid, top, bottom, N, M, i + 1, j + 1) + is_alive(grid, top, bottom, N, M, i, j - 1) + is_alive(grid, top, bottom, N, M, i, j + 1);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int recieve_configurations;
    int N, M, T;
    int **grid;

    int *flatten_send, *buffer;
    int counts_send[size], displacements[size], count_recieve;

    if (my_rank == 0)
    {
        cin >> N >> M >> T;
        flatten_send = (int *)malloc(N * M * sizeof(int));
        grid = (int **)malloc(N * sizeof(int *));
        for (int i = 0; i < N; ++i)
        {
            grid[i] = &flatten_send[M * i];
            for (int j = 0; j < M; j++)
            {
                cin >> grid[i][j];
            }
        }

        int block_size = N / size;
        int rem = N % size;

        for (int i = 0; i < size; i++)
        {
            if (rem)
            {
                counts_send[i] = (block_size + 1) * M;
                rem--;
            }
            else
            {
                if (i >= N)
                {
                    counts_send[i] = 0;
                }
                else
                    counts_send[i] = block_size * M;
            }
            if (i == 0 || i >= N)
            {
                displacements[i] = 0;
            }
            else
                displacements[i] = displacements[i - 1] + counts_send[i - 1];
        }
    }
    MPI_Bcast(&counts_send, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T, 1, MPI_INT, 0, MPI_COMM_WORLD);

    count_recieve = counts_send[my_rank];
    int n = count_recieve / M;

    int **block = (int **)malloc(n * sizeof(int *));
    int **simulated = (int **)malloc(n * sizeof(int *));
    int *flatten_simulated = (int *)malloc(count_recieve * sizeof(int));
    int *flatten_block = (int *)malloc(count_recieve * sizeof(int));
    int *bottom = NULL;
    int *top = NULL;

    if (my_rank > 0)
    {
        top = (int *)malloc(M * sizeof(int));
    }

    if (my_rank < size - 1)
    {
        bottom = (int *)malloc(M * sizeof(int));
    }

    buffer = (int *)malloc(count_recieve * sizeof(int));
    while (T--)
    {

        if (my_rank == 0)
            MPI_Scatterv(grid[0], counts_send, displacements, MPI_INT, buffer, count_recieve, MPI_INT, 0, MPI_COMM_WORLD);
        else
            MPI_Scatterv(NULL, NULL, NULL, MPI_INT, buffer, count_recieve, MPI_INT, 0, MPI_COMM_WORLD);

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

        int count;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < M; j++)
            {
                count = num_neighbours(block, top, bottom, n, M, i, j);
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
        }
        if (my_rank == 0)
        {

            MPI_Gatherv(simulated[0], count_recieve, MPI_INT, grid[0], counts_send, displacements, MPI_INT, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Gatherv(simulated[0], count_recieve, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);
        }
    }
    if (my_rank == 0)
    {
        cout << endl;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                cout << grid[i][j] << ' ';
            }
            cout << endl;
        }
        free(flatten_send);
        free(grid);
    }
    free(flatten_block);
    free(flatten_simulated);
    free(top);
    free(bottom);
    free(block);
    free(simulated);
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
