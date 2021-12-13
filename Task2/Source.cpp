#include "mpi.h" 
#include <Windows.h>
#include <vector>
#include <random>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
  int vect_size = 0; //Размер векторов
  int recvcount = 0; //Количество принимаемых элементов

  int* sendcounts = nullptr;
  int* displs = nullptr;

  double local_res = 0, global_res = 0;

  double t_start, t_stop;

  int ProcRank, ProcNum;
  MPI_Status status;

  MPI_Init(&argc, &argv); //Инизиализация МПИ
  MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
  MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

  std::vector<double> vect1, vect2, pers_vect1, pers_vect2;

  if (ProcRank == 0)
  {
    cout << "Vect size: ";
    cin >> vect_size;

    std::uniform_real_distribution<double> dist(0, 10000);
    std::default_random_engine re;

    vect1.reserve(vect_size);
    vect2.reserve(vect_size);

    for (int i = 0; i < vect_size; i++)
    {
      vect1[i] = dist(re);
      vect2[i] = dist(re);
    }

    int len = vect_size / ProcNum;
    int rem = vect_size % ProcNum;

    sendcounts = new int[ProcNum];
    displs = new int[ProcNum];

    sendcounts[0] = len + rem;
    displs[0] = 0;
    recvcount = sendcounts[0];

    for (int i = 1; i < ProcNum; i++)
    {
      sendcounts[i] = len;
      displs[i] = displs[i-1] + sendcounts[i-1];

      if (ProcRank == 0)
      {
        MPI_Send(&sendcounts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      }
    }
  }

  if (ProcRank != 0)
  {
    MPI_Recv(&recvcount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
  }

  pers_vect1.resize(recvcount);
  pers_vect2.resize(recvcount);

  MPI_Scatterv(&vect1[0], sendcounts, displs, MPI_DOUBLE, &pers_vect1[0], recvcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(&vect2[0], sendcounts, displs, MPI_DOUBLE, &pers_vect2[0], recvcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if (ProcRank == 0)
  {
    t_start = MPI_Wtime();
  }
  
  for (int i = 0; i < recvcount; i++)
  {
    local_res += pers_vect1[i] * pers_vect2[i];
  }

  MPI_Reduce(&local_res, &global_res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (ProcRank == 0)
  {
    t_stop = MPI_Wtime();
    std::cout << "Parallel  res: " << global_res << endl;
    std::cout << "Parallel time: " << t_stop - t_start << endl;

    double serial_res = 0;

    t_start = MPI_Wtime();

    for (int i = 0; i < vect_size; i++)
    {
      serial_res += vect1[i] * vect2[i];
    }

    t_stop = MPI_Wtime();
    std::cout << "  Serial  res: " << serial_res << endl;
    std::cout << "  Serial time: " << t_stop - t_start << endl;
  }

  MPI_Finalize();
  return 0;
}