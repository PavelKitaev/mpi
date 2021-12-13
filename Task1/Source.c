#include "mpi.h" 
#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>

long long Rubbish()
{
  const int n = 25;

  long long res = 0;
  int M1[25][25] = { 0 };
  int M2[25][25] = { 0 };
  int MR[25][25] = { 0 };

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      M1[i][j] = 100 + rand() % 10000;
    }
  }

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      M2[i][j] = 1000 + rand() % 10000;
    }
  }

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      MR[i][j] = 0;
      for (int k = 0; k < n; k++)
      {
        MR[i][j] += M1[i][k] * M2[k][j];
      }
      res += MR[i][j];
      //cout << MR[i][j] << " ";
    }
    //cout << endl;
  }
  return res;
}

double SerialAlg()
{
  const unsigned long n = 500000000;    //Количество прямоугольников
  double step;
  double PI25DT = 3.141592653589793238462643; //Табличное значение 
  double pi = 0;
  double sum = 0.0;
  double x;
  unsigned long i;

  //вычисляем шаг интегрирования и интегральную сумму
  step = 1.0 / (double)n;
  for (i = 1; i <= n; i++)
  {
    x = step * ((double)i - 0.5);
    sum += (4.0 / (1.0 + x * x));
  }
 
  pi = step * sum;

  return pi;
}

int main(int argc, char** argv)
{
  int ProcRank, ProcNum, done = 0, n = 500000000;
  unsigned long i;
  double PI25DT = 3.141592653589793238462643;
  double serial_pi, serial_time_start, serial_time_stop, serial_time;

  double parallel_time_start, parallel_time_stop, parallel_time;
  double mypi, pi, step, sum = 0.0, x;

  MPI_Init(&argc, &argv); //Инизиализация МПИ
  MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
  MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

  // главный цикл
  if (ProcRank == 0)
  {
    parallel_time_start = MPI_Wtime();
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (n > 0)
  {
    // вычисление локальных сумм
    step = 1.0 / (double)n;
    for (i = ProcRank + 1; i <= n; i += ProcNum)
    {
      x = step * ((double)i - 0.5);
      sum += (4.0 / (1.0 + x * x));
    }

    mypi = step * sum;
    MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  if (ProcRank == 0)
  {
    // вывод результата
    parallel_time_stop = MPI_Wtime();
    parallel_time = parallel_time_stop - parallel_time_start;
    printf("PARALLEL\n");
    printf("Reference PI:     %.25f\n", PI25DT);
    printf("Calculated value: %.25f\n", pi);
    printf("Error:            %.25f\n", fabs(pi - PI25DT));
    printf("Parallel time:    %f\n", parallel_time);
  }

  if (ProcRank == 0)
  {
    serial_time_start = MPI_Wtime();
    serial_pi = SerialAlg(); //Последовательная версия
    serial_time_stop = MPI_Wtime();

    serial_time = serial_time_stop - serial_time_start;
    printf("\n\n\nSERIAL\n");
    printf("Reference PI:     %.25f\n", PI25DT);
    printf("Calculated value: %.25f\n", serial_pi);
    printf("Error:            %.25f\n", fabs(serial_pi - PI25DT));
    printf("Serial_time:      %f\n", serial_time);
  }

  MPI_Finalize();
  return 0;
}