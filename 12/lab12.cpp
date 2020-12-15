#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <mpi.h>


void print_matrix(double *matrix, const int N, const int M)
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
      std::cout << *((matrix+i*M) + j) << " ";
    }
    std::cout << std::endl;
  }
}

void multiply_matrices(double *m1, const int N1, const int M1,
                          double *m2, const int N2, const int M2,
                          double *c)
{

  for (int i = 0; i < N1; i++)
  {
    for (int j = 0; j < M2; j++)
    {
      double sum = 0;
      for (int k = 0; k < M1; k++)
        sum += (*((m1+i*M1) + k)) * (*((m2+k*M2) + j));
      *((c + i*M2) + j) = sum;
    }
  }
}


int main(int argc, char **argv)
{ 
  const int AROWS = 4;
  const int ACOLS = 4;

  double A[AROWS][ACOLS] = {
    {1, 2, 1, 0},
    {3, 4, 0, 1},
    {1, 0, 2, 0},
    {0, 1, 0, 2}
  };

  const int BROWS = 4;
  const int BCOLS = 4;

  // double B[BROWS][BCOLS] = {
  //   {1, 2},
  //   {3, 4},
  //   {1, 0},
  //   {0, 1}
  // };

  double B[BROWS][BCOLS] = {
    {1, 1, 1, 1},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}
  };

  double C[AROWS][BCOLS];

  const int LOCAL_A_ROWS = 2;
  const int LOCAL_A_COLS = 2;

  const int LOCAL_B_ROWS = 2;
  const int LOCAL_B_COLS = 2;

  const int LOCAL_C_ROWS = LOCAL_A_ROWS;
  const int LOCAL_C_COLS = LOCAL_B_COLS;

  double localA[LOCAL_A_ROWS][LOCAL_A_COLS];
  double localB[LOCAL_B_ROWS][LOCAL_B_COLS];
  double localC[LOCAL_C_ROWS][LOCAL_C_COLS];

  const int DECOMP_A_ROWS = 2;
  const int DECOMP_A_COLS = 2;
  const int DECOMP_B_ROWS = 2;
  const int DECOMP_B_COLS = 2;
  const int DECOMP_C_ROWS = 2;
  const int DECOMP_C_COLS = 2;


  int process_id;
  int size;


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  MPI_Comm comm_3D;
	int const ndim = 3;
	int dims[ndim] = {2, 2, 2};
	int periods[ndim] = {0, 0, 0};

  MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, 0, &comm_3D);


  MPI_Comm comm_xy;
  MPI_Comm comm_yz;
  MPI_Comm comm_xz;

  int xy_belongs[] = {1, 1, 0};
  int yz_belongs[] = {0, 1, 1};
  int xz_belongs[] = {1, 0, 1};

  MPI_Cart_sub(comm_3D, xy_belongs, &comm_xy);
  MPI_Cart_sub(comm_3D, yz_belongs, &comm_yz);
  MPI_Cart_sub(comm_3D, xz_belongs, &comm_xz);

  MPI_Comm comm_z;
  MPI_Comm comm_x;
  MPI_Comm comm_y;

  int z_belongs[] = {0, 0, 1};
  int x_belongs[] = {1, 0, 0};
  int y_belongs[] = {0, 1, 0};

  MPI_Cart_sub(comm_3D, z_belongs, &comm_z);
  MPI_Cart_sub(comm_3D, x_belongs, &comm_x);
  MPI_Cart_sub(comm_3D, y_belongs, &comm_y);


  const int A_BLOCK_ROWS = 2;
  const int A_BLOCK_COLS = 2;
  const int B_BLOCK_ROWS = BROWS / 2;
  const int B_BLOCK_COLS = BCOLS / 2;
  const int C_BLOCK_ROWS = A_BLOCK_ROWS;
  const int C_BLOCK_COLS = B_BLOCK_COLS;
  const int A_COLS = 4;

  MPI_Datatype blockA;
  MPI_Datatype blockA2;

  MPI_Datatype blockB;
  MPI_Datatype blockB2;

  MPI_Datatype blockC;
  MPI_Datatype blockC2;

  MPI_Type_vector(A_BLOCK_ROWS, A_BLOCK_COLS, A_COLS, MPI_DOUBLE, &blockA2);
  MPI_Type_create_resized(blockA2, 0, sizeof(double), &blockA);
  MPI_Type_commit(&blockA);

  int sizes[] = {BROWS, BCOLS};
  int subsizes[] = {B_BLOCK_ROWS, B_BLOCK_COLS};
  int starts [] = {0, 0};
  MPI_Type_create_subarray(2, sizes, subsizes,
                            starts, MPI_ORDER_C, MPI_DOUBLE, &blockB2);
  MPI_Type_create_resized(blockB2, 0, sizeof(double), &blockB);
  MPI_Type_commit(&blockB);

  MPI_Type_vector(C_BLOCK_ROWS, C_BLOCK_COLS, BCOLS, MPI_DOUBLE, &blockC2);
  MPI_Type_create_resized(blockC2, 0, sizeof(double), &blockC);
  MPI_Type_commit(&blockC);


  int disps[DECOMP_A_ROWS*DECOMP_A_COLS];
  int counts[DECOMP_A_ROWS*DECOMP_A_COLS];
  for (int ii=0; ii<DECOMP_A_ROWS; ii++) {
      for (int jj=0; jj<DECOMP_A_COLS; jj++) {
          disps[ii*DECOMP_A_COLS+jj] = ii*4*2+jj*2;
          counts [ii*DECOMP_A_COLS+jj] = 1;
      }
  }

  int dispsB[DECOMP_B_ROWS*DECOMP_B_COLS];
  int countsB[DECOMP_B_ROWS*DECOMP_B_COLS];
  for (int ii=0; ii<DECOMP_B_ROWS; ii++) {
      for (int jj=0; jj<DECOMP_B_COLS; jj++) {
          dispsB[ii*DECOMP_B_COLS+jj] = ii*BCOLS*B_BLOCK_ROWS+jj*B_BLOCK_COLS;
          countsB[ii*DECOMP_B_COLS+jj] = 1;
      }
  }

  int dispsC[DECOMP_C_ROWS*DECOMP_C_COLS];
  int countsC[DECOMP_C_ROWS*DECOMP_C_COLS];
  for (int ii=0; ii<DECOMP_C_ROWS; ii++) {
      for (int jj=0; jj<DECOMP_C_COLS; jj++) {
          dispsC[ii*DECOMP_C_COLS+jj] = ii*BCOLS*C_BLOCK_ROWS+jj*C_BLOCK_COLS;
          countsC[ii*DECOMP_C_COLS+jj] = 1;
      }
  }
  

  int coords[3];
  MPI_Cart_coords(comm_3D, process_id, 3, coords);

  if (process_id == 0)
  {
    MPI_Scatterv(&A, counts, disps, blockA, localA, 2*2, MPI_DOUBLE, 0, comm_xy);

    MPI_Scatterv(&B, countsB, dispsB, blockB, localB, B_BLOCK_ROWS*B_BLOCK_COLS, MPI_DOUBLE, 0, comm_yz);
  }
  else
  {

    if (coords[2] == 0) 
    {
      MPI_Scatterv(&A, counts, disps, blockA, &localA, 2*2, MPI_DOUBLE, 0, comm_xy);
    }
    if (coords[0] == 0)
    {

      MPI_Scatterv(&B, countsB, dispsB, blockB, localB, B_BLOCK_ROWS*B_BLOCK_COLS, MPI_DOUBLE, 0, comm_yz);
    }
  }


  int rank_in_y;
  MPI_Comm_rank(comm_y, &rank_in_y);


  MPI_Bcast(&localA, LOCAL_A_ROWS*LOCAL_A_COLS, MPI_DOUBLE, 0, comm_z);
  MPI_Bcast(&localB, B_BLOCK_ROWS*B_BLOCK_COLS, MPI_DOUBLE, 0, comm_x);

  multiply_matrices((double *)localA, LOCAL_A_ROWS, LOCAL_A_COLS,
                    (double *)localB, LOCAL_B_ROWS, LOCAL_B_COLS, 
                    (double *)localC);

  if (rank_in_y == 0) 
  {
    MPI_Reduce(MPI_IN_PLACE, &localC, LOCAL_C_ROWS*LOCAL_C_COLS, MPI_DOUBLE, MPI_SUM, 0, comm_y);
  }
  else
  {
    MPI_Reduce(&localC, NULL, LOCAL_C_ROWS*LOCAL_C_COLS, MPI_DOUBLE, MPI_SUM, 0, comm_y);
  }


  if (coords[1] == 0)
  {
    MPI_Gatherv(&localC, LOCAL_C_ROWS*LOCAL_C_COLS, MPI_DOUBLE, &C, countsC, dispsC, blockC, 0, comm_xz);
  }
  
  
  if (process_id == 0)
  {
    std::cout << "Result C = " << std::endl;
    print_matrix((double *)C, BROWS, BCOLS);
  }


  MPI_Type_free(&blockB);
  MPI_Type_free(&blockA);
  MPI_Type_free(&blockC);
  MPI_Finalize();
  return 0;
}