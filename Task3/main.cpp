#include "mpi.h" 
#include <iostream>
#include <random>

using namespace std;

void Rubbish()
{
  std::uniform_real_distribution<double> dist(800000, 900000);
  std::default_random_engine re;
  double a_random_double = dist(re);

  double** matrix = new double* [8000];

  for (int i = 0; i < 5000; i++)
  {
    matrix[i] = new double[5000];
    for (int j = 0; j < 5000; j++)
    {
      matrix[i][j] = dist(re);
    }
  }

  double* vect = new double[5000];
  for (int i = 0; i < 5000; i++)
  {
    vect[i] = dist(re);
  }

  double* res = new double[5000] {0};

  for (int i = 0; i < 5000; i++)
  {
    for (int j = 0; j < 5000; j++)
    {
      res[i] += matrix[i][j] * vect[j];
    }
  }

  delete[] matrix;
  delete[] vect;
  delete[] res;
}

void GetData(int* mass, int size_m, int size_n, int ProcNum, int* recvcounts_gatherv, int* displ_gatherv)
{

  int control = 0; //������� ����������� ���������

  if (size_m % ProcNum != 0) //���� ���������� ��������� �� ������� ��� �������
  {
    int temp = size_m - 1;

    for (int i = 1; i < size_m; i++) //���������� �����, ������� ����� �������� ��� �������
    {
      if (temp % ProcNum == 0)
      {
        break;
      }
      else
      {
        temp = temp - 1;
      }
    }

    for (int i = 0; i < ProcNum; i++) //��������� ������� ��� ������� ��������
    {
      mass[i] = temp / ProcNum;
      control += mass[i];
    }

    //���������� ������������ �������
    int iter = 0;
    int residue = size_m - temp;

    while (1)
    {
      if (residue != 0) //���� ������� ��� �������
      {
        if (iter < ProcNum)
        {
          mass[iter] += 1;
        }
        else //���� ������� ����� �� ������� - ��������
        {
          iter = 0;
        }

        control++;
        iter++;
        residue--;
      }
      else //���� ���� ������� ����������� - �������
      {
        break;
      }
    }
  }
  else //���� ���������� ��������� ������� ��� �������
  {
    for (int i = 0; i < ProcNum; i++)
    {
      mass[i] = size_m / ProcNum;
      control += mass[i];
    }
  }

  if (control != size_m) //���������, �������� �� ���������� ����������� ��������� � �������� �����������
  {
    cout << control << " = " << size_m << endl;
    throw - 1;
  }

  for (int i = 0; i < ProcNum; i++)
  {
    recvcounts_gatherv[i] = mass[i];
    mass[i] = mass[i] * size_n;
  }

  displ_gatherv[0] = 0;
  for (int i = 1; i < ProcNum; i++)
  {
    displ_gatherv[i] = recvcounts_gatherv[i - 1] + displ_gatherv[i - 1];
  }

}

void Displs(int* displ, int* counts, int ProcNum)
{
  displ[0] = 0;

  for (int i = 1; i < ProcNum; i++)
  {
    displ[i] = counts[i - 1] + displ[i - 1];
  }
}

void SerialAlg(double* matrix, double* vect, double* res, int size_m, int size_n)
{
  int q = 0;
  for (int i = 0; i < size_m; i++)
  {
    for (int j = 0; j < size_n; j++)
    {
      res[i] += matrix[q] * vect[j];
      q++;
    }
  }
}

int main(int argc, char** argv)
{
  double t_start, t_stop;

  int ProcRank, ProcNum;
  MPI_Status status;
  MPI_Init(&argc, &argv); //������������� ���
  MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
  MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

  double* global_matrix = nullptr;    //�������� �������
  double* local_matrix = nullptr;     //�������������� ������� ��������
  int recvcount = 0;                  //������ ����������� �������
  int* sendcounts = nullptr;          //���������� ��������� �������, ������������ ������� ��������
  int* displs = nullptr;              //�������� ��� ��������
  int size_m = 0, size_n;             //������ �������

  double* vect = nullptr;             //������

  double* local_res = nullptr;        //���������, ����������� ����� ���������
  double* global_res = nullptr;       //����
  int* recvcounts_gatherv = nullptr;  //������� ������
  int* displ_gatherv = nullptr;       //�������� ��� �����

  if (ProcRank == 0)
  {
    std::uniform_real_distribution<double> dist(800000, 900000);
    std::default_random_engine re;

    cout << "Row: ";
    cin >> size_m; //������
    cout << "Column: ";
    cin >> size_n; //�������

    if (size_m < ProcNum)
    {
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    //-----------------------------------------------------------------������������ ������� � �������
    vect = new double[size_n];
    for (int i = 0; i < size_n; i++)
    {
      vect[i] = dist(re);
    }

    int size_matrix = size_m * size_n;

    global_matrix = new double[size_matrix];
    int q = 0;
    for (int i = 0; i < size_m; i++)
    {
      for (int j = 0; j < size_n; j++)
      {
        global_matrix[q] = dist(re);
        q++;
      }
    }

    cout << "Start parallel alg" << endl;
    t_start = MPI_Wtime();

    sendcounts = new int[ProcNum];  //���������� ��������� ��� ������� ��������
    displs = new int[ProcNum];      //�������� (��������)

    recvcounts_gatherv = new int[ProcNum]; //���������� ��������� �� ������� �������� (����)
    displ_gatherv = new int[ProcNum];   //�������� (����)

    GetData(sendcounts, size_m, size_n, ProcNum, recvcounts_gatherv, displ_gatherv);  //���������� ���������� ��������� ��� �����, ������, ������ �����
    Displs(displs, sendcounts, ProcNum);  //���������� ������ ��������

    recvcount = sendcounts[ProcRank]; //���������� ����� ����������� �������

    char buff[1000];
    for (int i = 1; i < ProcNum; i++)
    {
      int pos = 0;
      MPI_Pack(&size_n, 1, MPI_INT, buff, 1000, &pos, MPI_COMM_WORLD);          //���������� ��������
      MPI_Pack(&sendcounts[i], 1, MPI_INT, buff, 1000, &pos, MPI_COMM_WORLD);   //����� ���������� ��������� �����

      MPI_Send(buff, pos, MPI_PACKED, i, 0, MPI_COMM_WORLD);
    }

    global_res = new double[size_m];  //������ ��� ���������� 
  }

  if (ProcRank != 0) //����� ���������� �������� � ���������� ��������� �������
  {
    int temp[2];
    MPI_Recv(&temp, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

    size_n = temp[0];
    recvcount = temp[1];

    vect = new double[size_n] {0};
  }

  local_matrix = new double[recvcount] { 0 };

  MPI_Scatterv(&global_matrix[0], sendcounts, displs, MPI_DOUBLE, &local_matrix[0], recvcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&vect[0], size_n, MPI_DOUBLE, 0, MPI_COMM_WORLD); //�������� �������

  //----------------------------------------------------------��������� ����� ������� �� ������
  int count_row = recvcount / size_n;             //���������� ���������� �����
  local_res = new double[count_row] {0};          //���������

  int q = 0;
  for (int i = 0; i < count_row; i++)
  {
    for (int j = 0; j < size_n; j++)
    {
      local_res[i] += local_matrix[q] * vect[j];
      q++;
    }
  }
  
  //���� ������
  MPI_Gatherv(&local_res[0], count_row, MPI_DOUBLE, &global_res[0],
    recvcounts_gatherv, displ_gatherv, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (ProcRank == 0)
  {
    t_stop = MPI_Wtime();
    cout << "Parallel time: " << t_stop - t_start << endl;

    //������������ ���������
    
    //for (int i = 0; i < size_m; i++)
    //{
    //  cout << global_res[i] << " ";
    //}
    //cout << endl;

    Rubbish(); //�����

    double* res = new double[size_m] {0};
    cout << "Start serial alg" << endl;
    t_start = MPI_Wtime();
    SerialAlg(global_matrix, vect, res, size_m, size_n);
    t_stop = MPI_Wtime();
    cout << "Serial time: " << t_stop - t_start << endl;

    //���������������� ���������
 
    //for (int i = 0; i < size_m; i++)
    //{
    //  cout << res[i] << " ";
    //}
    //cout << endl;

    delete[] global_res;
    delete[] global_matrix;
    delete[] sendcounts;
    delete[] displs;
    delete[] recvcounts_gatherv;
    delete[] displ_gatherv;
  }

  delete[] vect;
  delete[] local_res;
  delete[] local_matrix;

  MPI_Finalize();
  return 0;
}