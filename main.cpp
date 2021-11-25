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

  int control = 0; //Подсчет разделенных элементов

  if (size_m % ProcNum != 0) //Если количество элементов не делится без остатка
  {
    int temp = size_m - 1;

    for (int i = 1; i < size_m; i++) //Вычисление числа, которое будет делиться без остатка
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

    for (int i = 0; i < ProcNum; i++) //Вычисляем размеры для каждого процесса
    {
      mass[i] = temp / ProcNum;
      control += mass[i];
    }

    //Равномерно распределяем остаток
    int iter = 0;
    int residue = size_m - temp;

    while (1)
    {
      if (residue != 0) //Если остаток еще имеется
      {
        if (iter < ProcNum)
        {
          mass[iter] += 1;
        }
        else //Если счетчик вышел за пределы - обнулить
        {
          iter = 0;
        }

        control++;
        iter++;
        residue--;
      }
      else //Если весь остаток распределен - выходим
      {
        break;
      }
    }
  }
  else //Если количество элементов делится без остатка
  {
    for (int i = 0; i < ProcNum; i++)
    {
      mass[i] = size_m / ProcNum;
      control += mass[i];
    }
  }

  if (control != size_m) //Проверяем, сходится ли количество разделенных элементов с исходным количеством
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
  MPI_Init(&argc, &argv); //Инизиализация МПИ
  MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
  MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

  double* global_matrix = nullptr;    //Основная матрица
  double* local_matrix = nullptr;     //Индивидуальная матрица процесса
  int recvcount = 0;                  //Размер принимаемой матрицы
  int* sendcounts = nullptr;          //Количество элементов матрицы, отправляемых каждому процессу
  int* displs = nullptr;              //Смещение для рассылки
  int size_m = 0, size_n;             //Размер матрицы

  double* vect = nullptr;             //Вектор

  double* local_res = nullptr;        //Результат, посчитанный одним процессом
  double* global_res = nullptr;       //Итог
  int* recvcounts_gatherv = nullptr;  //Размеры приема
  int* displ_gatherv = nullptr;       //Смещение для сбора

  if (ProcRank == 0)
  {
    std::uniform_real_distribution<double> dist(800000, 900000);
    std::default_random_engine re;

    cout << "Row: ";
    cin >> size_m; //Строки
    cout << "Column: ";
    cin >> size_n; //Столбцы

    if (size_m < ProcNum)
    {
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    //-----------------------------------------------------------------Формирование матрицы и вектора
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

    sendcounts = new int[ProcNum];  //Количество элементов для каждого процесса
    displs = new int[ProcNum];      //Смещение (отправка)

    recvcounts_gatherv = new int[ProcNum]; //Количество элементов от каждого процесса (сбор)
    displ_gatherv = new int[ProcNum];   //Смещение (сбор)

    GetData(sendcounts, size_m, size_n, ProcNum, recvcounts_gatherv, displ_gatherv);  //Вычисление количества элементов для сбора, приема, сдвига сбора
    Displs(displs, sendcounts, ProcNum);  //Вычисление сдвига отправки

    recvcount = sendcounts[ProcRank]; //Количество строк принимаемой матрицы

    char buff[1000];
    for (int i = 1; i < ProcNum; i++)
    {
      int pos = 0;
      MPI_Pack(&size_n, 1, MPI_INT, buff, 1000, &pos, MPI_COMM_WORLD);          //Количество столбцов
      MPI_Pack(&sendcounts[i], 1, MPI_INT, buff, 1000, &pos, MPI_COMM_WORLD);   //Общее количество элементов строк

      MPI_Send(buff, pos, MPI_PACKED, i, 0, MPI_COMM_WORLD);
    }

    global_res = new double[size_m];  //Массив для результата 
  }

  if (ProcRank != 0) //Прием количества столбцов и количества элементов матрицы
  {
    int temp[2];
    MPI_Recv(&temp, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

    size_n = temp[0];
    recvcount = temp[1];

    vect = new double[size_n] {0};
  }

  local_matrix = new double[recvcount] { 0 };

  MPI_Scatterv(&global_matrix[0], sendcounts, displs, MPI_DOUBLE, &local_matrix[0], recvcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&vect[0], size_n, MPI_DOUBLE, 0, MPI_COMM_WORLD); //Рассылка вектора

  //----------------------------------------------------------Умножение части матрицы на вектор
  int count_row = recvcount / size_n;             //Количество полученных строк
  local_res = new double[count_row] {0};          //Результат

  int q = 0;
  for (int i = 0; i < count_row; i++)
  {
    for (int j = 0; j < size_n; j++)
    {
      local_res[i] += local_matrix[q] * vect[j];
      q++;
    }
  }
  
  //Сбор данных
  MPI_Gatherv(&local_res[0], count_row, MPI_DOUBLE, &global_res[0],
    recvcounts_gatherv, displ_gatherv, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (ProcRank == 0)
  {
    t_stop = MPI_Wtime();
    cout << "Parallel time: " << t_stop - t_start << endl;

    //Параллельный результат
    
    //for (int i = 0; i < size_m; i++)
    //{
    //  cout << global_res[i] << " ";
    //}
    //cout << endl;

    Rubbish(); //Мусор

    double* res = new double[size_m] {0};
    cout << "Start serial alg" << endl;
    t_start = MPI_Wtime();
    SerialAlg(global_matrix, vect, res, size_m, size_n);
    t_stop = MPI_Wtime();
    cout << "Serial time: " << t_stop - t_start << endl;

    //Последовательный результат
 
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