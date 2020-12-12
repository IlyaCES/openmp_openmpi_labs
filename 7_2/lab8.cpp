#include <iostream>
#include <cmath>


int main(int argc, char **argv)
{
  // const int n_sets = 3;
  // const int set_dim = 2;
  // const int input_sets[n_sets][set_dim] = {{1, 2}, {3, 4}, {5, 6}};

  // const int n_sets = 2;
  // const int set_dim = 3;
  // const int input_sets[n_sets][set_dim] = {{1, 2, 3}, {4, 5, 6}};

  int n_sets = 6;
  int set_dim = 7;
  int input_sets[n_sets][set_dim];
  for (int i = 0; i < n_sets; i++)
  {
    for (int j = 0; j < set_dim; j++)
    {
      input_sets[i][j] = i + j;
    }
  }

  const int result_len = static_cast<int>(pow(set_dim, n_sets));
  int result[result_len][n_sets];

  # pragma omp parallel for firstprivate(set_dim, n_sets) collapse(2)
  for (int i = 0; i < result_len; i++)
  {
    for (int j = 0; j < n_sets; j++) 
    {
      const int index = (i / (static_cast<int>(pow(set_dim, j)))) % set_dim;
      result[i][j] = input_sets[j][index];
    }
  }


  // for(int i=0; i<result_len; i++)
  // {
  //  for(int j=0; j<n_sets; j++)
  //    {
  //     std::cout<<" "<<result[i][j]<<" ";
  //    }
  //  std::cout<<"\n";
  // }
}