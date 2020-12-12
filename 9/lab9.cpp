#include <iostream>
#include <cmath>
#include <mpi.h>


int main(int argc, char **argv)
{
  int process_id;
  int size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  srand (process_id);


  if (process_id == 0)
  {
    const int n_bookkeepers = size - 1;
    int report[n_bookkeepers];

    for (int i = 0; i < n_bookkeepers; i++)
    {
      MPI_Status status;
      int data[2];

      int ierr = MPI_Recv(&data, 2, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

      if (ierr != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, 1);

      const int income = data[0];
      const int expenses = data[1];

      std::cout << "Received income = " << income << " expenses = " << expenses << std::endl;

      report[i] = income - expenses;
    }

    int sum = 0;
    std::cout << "Report:" << std::endl;
    for (int i = 0; i < n_bookkeepers; i++)
    {
      std::cout << report[i] << " ";
      sum += report[i];
    }
    std::cout << std::endl << "Total sum = " << sum << std::endl;
  }
  else
  {
    const int income = rand() % 100 + 1;
    const int expenses = rand() % 50 + 1;

    const int data[] = {income, expenses};
    MPI_Send(&data, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
    std::cout << "Sended income = " << income << " expenses = " << expenses << std::endl;
  }

  MPI_Finalize();
  return 0;
}