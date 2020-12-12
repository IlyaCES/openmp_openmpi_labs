#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <mpi.h>


void scatterv(const int data[], const int sendcounts[], const int displs[], 
              int recvbuffer[], const int recvsize, const int root)
{
  int process_id;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (process_id == root)
  {
    for (int i = 0; i < size; i++)
    {
      if (i != process_id)
      {
        const int from = displs[i];
        const int to = from + sendcounts[i];

        int arrayToSend[sendcounts[i]];
        for (int j = from, k = 0; j < to, k < sendcounts[i]; j++, k++)
        {
          arrayToSend[k] = data[j];
        }

        MPI_Send(&arrayToSend, sendcounts[i], MPI_INT, i, 0, MPI_COMM_WORLD);

      }
      const int from = displs[process_id];
      const int to = from + sendcounts[process_id];

      for (int j = from, k = 0; j < to, k < sendcounts[i]; j++, k++)
      {
          recvbuffer[k] = data[j];
      }
    }
  }
  else
  {
    MPI_Status status;
    int ierr = MPI_Recv(recvbuffer, sendcounts[process_id], MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

    if (ierr != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, 1);
  }
}


int main(int argc, char **argv)
{
  int process_id;
  int size;
  int data[] = {1, 2 , 3, 4 , 5, 6, 7};
  int displs[] = {0, 2, 5};
  int sendcounts[] = {2, 3, 2};
  int root = 1;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int recvbuffer[sendcounts[process_id]];

  // scatterv(data, sendcounts, displs, recvbuffer, sendcounts[process_id], root);
  // MPI_Scatterv(data, sendcounts, displs, MPI_INT,
  //              recvbuffer, sendcounts[process_id], MPI_INT, root, MPI_COMM_WORLD);

  for (int i = 0; i < 100000000; i++)
  {
    // MPI_Scatterv(data, sendcounts, displs, MPI_INT,
    //            recvbuffer, sendcounts[process_id], MPI_INT, root, MPI_COMM_WORLD);
    scatterv(data, sendcounts, displs, recvbuffer, sendcounts[process_id], root);
  }
  
  std::cout << "Process " << process_id << " received ";
  for (int k = 0; k < sendcounts[process_id]; k++)
  {
    std::cout << " " << recvbuffer[k] << " ";
  }
  std::cout << std::endl;

  MPI_Finalize();
  return 0;
}