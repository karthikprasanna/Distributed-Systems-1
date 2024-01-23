#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

void calculateShortestDistance(vector<vector<int>> A, vector<vector<int>> A_Final, int k)
{
    int N = A.size();
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (A_Final[i][k] != INT_MAX && A_Final[k][j] != INT_MAX)
                A_Final[i][j] = min(A_Final[i][j], A_Final[i][k] + A_Final[k][j]);
            
        }
        
    }
    
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int P = size;
    int N;
    cin>>N;

    vector<vector<int>> A(N, vector<int>(N));
    vector<vector<int>> A_Final(N, vector<int>(N));

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            cin>>A[i][j];
            if (A[i][j] == -1)
            {
                A_Final[i][j] = INT_MAX;
            }
            else
            {
                A_Final[i][j] = A[i][j];
            }
        }
    }
    for (int k = 0; k < N; k++)
    {
        calculateShortestDistance(A, A_Final, k);
    }

    if (my_rank == 0)
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                if (A_Final[i][j] != INT_MAX)
                {
                    cout<<A_Final[i][j]<<" ";
                }
                else
                {
                    cout<<"-1 ";
                }                

            }
            cout<<endl;
        }   
    }
    
    MPI_Finalize();
    return 0;    
}

