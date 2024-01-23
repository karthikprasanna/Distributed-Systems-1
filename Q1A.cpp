#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;
int configurations = 0;
bool isValid(vector<string> board, int i, int j) {
    int n = board.size();
    // check coulumn
    for (int k = i + 1; k < n; ++k) {
        if (board[k][j] == 'Q') {
            return false;
        }
    }

    // check main diagonal
    for (int k = i + 1, l = j + 1; k < n and l < n; ++k, ++l) {
        if (board[k][l] == 'Q') {
            return false;
        }
    }

    // check alternative diagonal
    for (int k = i + 1, l = j - 1; k < n and l >= 0; ++k, --l) {
        if (board[k][l] == 'Q') {
            return false;
        }
    }

    return true;

}

void solver(vector<string> board, int row) {
    if (row == -1) {
        configurations++;
        return;
    }

    int n = board.size();

    for (int i = 0; i < n; ++i) {
        if (isValid(board, row, i)) {
            board[row][i] = 'Q';
            solver(board, row - 1);
            board[row][i] = '.';
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

    int recieve_configurations;
    int n;
    if (my_rank == 0)
    {
        cin>>n;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    string str = "";
    for (int i = 0; i < n; ++i) {
        str.push_back('.');
    }
    vector<string> board(n, str);
    
    int block_size = n / size;

    if (block_size * size != n)
    {
        block_size++;
    }

    for (int i = my_rank * block_size; i < (my_rank + 1) * block_size && i < n; i++)
    {
    board[n - 1][i] = 'Q';
    solver(board, n - 2);
    board[n - 1][i] = '.';
    }
    

    MPI_Reduce(&configurations, &recieve_configurations, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        cout<<recieve_configurations<<endl;
    }
    MPI_Finalize();
    return 0;    
}

