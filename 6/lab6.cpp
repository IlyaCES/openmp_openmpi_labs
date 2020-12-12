#include <iostream>


int f(const int array[], const int size)
{
  int result = 0;
  #pragma omp parallel for shared(array) reduction(+:result)
  for (int i = 0; i < size - 1; i += 2)
  {
    result += array[i] - array[i+1];
  }
  return result;
}


int main(int argc, char **argv)
{
  const int n = 100;
  int A[n];

  for (int i = 0; i < n; i++)
  {
    A[i] = n - i;
  }

  const int result = f(A, n);
  std::cout << "result: " << result << "\n";
}