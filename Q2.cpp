#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int N;
    int **grid;

    int *flatten_send, *buffer;
    int counts_send[size], displacements[size], count_recieve;

    if (my_rank == 0)
    {
        cin >> N;
        flatten_send = (int *)malloc(N * N * sizeof(int));
        grid = (int **)malloc(N * sizeof(int *));

        for (int i = 0; i < N; ++i)
        {
            grid[i] = &flatten_send[N * i];
            for (int j = 0; j < N; j++)
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
                counts_send[i] = (block_size + 1) * N;
                rem--;
            }
            else
            {
                if (i >= N)
                {
                    counts_send[i] = 0;
                }
                else
                    counts_send[i] = block_size * N;
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

    count_recieve = counts_send[my_rank];
    int n = count_recieve / N;

    int **block = (int **)malloc(n * sizeof(int *));
    int *flatten_block = (int *)malloc(count_recieve * sizeof(int));

    int *row = (int *)malloc(N * sizeof(int));

    int *column = (int *)malloc(N * sizeof(int));

    buffer = (int *)malloc(count_recieve * sizeof(int));

    for (int k = 0; k < N; k++)
    {

        if (my_rank == 0)
        {
            for (int i = 0; i < N; i++)
            {
                row[i] = grid[k][i];
                column[i] = grid[i][k];
            }

            MPI_Scatterv(grid[0], counts_send, displacements, MPI_INT, buffer, count_recieve, MPI_INT, 0, MPI_COMM_WORLD);
        }
        else
            MPI_Scatterv(NULL, NULL, NULL, MPI_INT, buffer, count_recieve, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Bcast(row, N, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(column, N, MPI_INT, 0, MPI_COMM_WORLD);
        for (int i = 0; i < n; i++)
        {
            block[i] = flatten_block + i * N;
            for (int j = 0; j < N; j++)
            {
                block[i][j] = buffer[i * N + j];
            }
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (column[i] != -1 && row[j] != -1)
                {
                    if (block[i][j] == -1)
                    {
                        block[i][j] = column[i] + row[j];
                    }
                    else
                        block[i][j] = min(block[i][j], column[i] + row[j]);
                }
            }
        }

        if (my_rank == 0)
        {

            MPI_Gatherv(block[0], count_recieve, MPI_INT, grid[0], counts_send, displacements, MPI_INT, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Gatherv(block[0], count_recieve, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);
        }
    }

    if (my_rank == 0)
    {
        cout << endl;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                cout << grid[i][j] << ' ';
            }
            cout << endl;
        }
        free(flatten_send);
        free(grid);
    }

    free(flatten_block);
    free(row);
    free(column);
    free(block);
    free(buffer);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

/*
2
0 25
-1 0

3
0 1 43
1 0 6
-1 -1 0

4
0 5 -1 10
-1 0 3 -1
-1 -1 0 1
-1 -1 -1 0
*/
