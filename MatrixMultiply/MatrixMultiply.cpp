#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include "mpi.h"
#include <malloc.h>
#include <fstream>
#include <iostream>

#define M 4
#define N 4
int main()
{
	int my_rank;/*My process rank*/
	int comm_sz;/*Number of processes*/
	int local_M;
	int i, j, k;
	double start, finish;/*timer*/
	int tem;
	//Инициализировать MPI
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	//Количество строк, выделенных для каждой матрицы
	local_M = M / comm_sz;
	//Матрица присваивается каждому процессу
	int* local_Matrix_one = (int*)malloc(local_M * N * sizeof(int));
	//Определить две матрицы
	int* Matrix_one = NULL;
	int* Matrix_two = (int*)malloc(M * N * sizeof(int));
	//Матрица результатов в каждом процессе
	int* local_result = (int*)malloc(local_M * N * sizeof(int));
	//Матрица результатов
	int* result_Matrix = NULL;
	if (my_rank == 0)
	{
		//printf("process %d of %d\n",my_rank,comm_sz);
		FILE* fp;
		//Прочитайте первую матрицу
		Matrix_one = (int*)malloc(M * N * sizeof(int));
		//Matrix_one[M][N]={0};
		fp = fopen("a.txt", "r");
		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
				fscanf(fp, "%d ", &Matrix_one[i * N + j]);
			fscanf(fp, "\n");
		}
		fclose(fp);
		/*for(i=0;i<M;i++)
		{
		for(j=0;j<N;j++)
		   printf("%d ",Matrix_one[i*N+j]);
		   printf("\n");
		}*/
		// Читаем вторую матрицу
		start = MPI_Wtime();
		fp = fopen("b.txt", "r");
		for (j = 0; j < N; j++)
		{
			for (i = 0; i < M; i++)
				fscanf(fp, "%d ", &Matrix_two[i * N + j]);
			fscanf(fp, "\n");
		}
		fclose(fp);
		/*for(j=0;j<N;j++)
		{
		for(i=0;i<M;i++)
		   printf("%d ",Matrix_two[i*M+j]);
		   printf("\n");
		}*/
		// Распределение данных
		MPI_Scatter(Matrix_one, local_M * N, MPI_INT, local_Matrix_one, local_M * N, MPI_INT, 0, MPI_COMM_WORLD);

		//Передача данных
		MPI_Bcast(Matrix_two, M * N, MPI_INT, 0, MPI_COMM_WORLD);
		//расчетlocalрезультат
		for (i = 0; i < local_M; i++)
			for (j = 0; j < M; j++) {
				tem = 0;
				for (k = 0; k < N; k++)
					tem += local_Matrix_one[i * M + k] * Matrix_two[j * M + k];
				local_result[i * M + j] = tem;
			}
		free(local_Matrix_one);
		result_Matrix = (int*)malloc(M * N * sizeof(int));
		//Сбор результатов
		MPI_Gather(local_result, local_M * N, MPI_INT, result_Matrix, local_M * N, MPI_INT, 0, MPI_COMM_WORLD);
		//Обработка оставшейся строки (обработка ситуации, которая не может быть делимой)
		int rest = M % comm_sz;
		if (rest != 0)
			for (i = M - rest - 1; i < M; i++)
				for (j = 0; j < M; j++) {
					tem = 0;
					for (k = 0; k < N; k++)
						tem += Matrix_one[i * M + k] * Matrix_two[j * M + k];
					result_Matrix[i * M + j] = tem;
				}
		finish = MPI_Wtime();
		free(Matrix_one);
		free(Matrix_two);
		free(local_result);

		printf("Proc %d > Elapsed time = %e seconds\n", my_rank, finish - start);
		//Записать результаты в файл
		std::ofstream fileout;
		fileout.open("c.txt", std::ios::out | std::ios::trunc);
		if (fileout.is_open())
		{
			for (i = 0; i < M; i++)
			{
				for (j = 0; j < N; j++)
					fileout << result_Matrix[i * M + j] << " ";
				fileout << std::endl;
			}
		}
		fileout.close();
	}
	else {
		//printf("process %d of %d\n",my_rank,comm_sz);
		//Распределение данных
		MPI_Scatter(Matrix_one, local_M * N, MPI_INT, local_Matrix_one, local_M * N, MPI_INT, 0, MPI_COMM_WORLD);
		//Передача данных
		MPI_Bcast(Matrix_two, M * N, MPI_INT, 0, MPI_COMM_WORLD);
		//расчетlocalрезультат
		for (i = 0; i < local_M; i++)
			for (j = 0; j < M; j++) {
				tem = 0;
				for (k = 0; k < N; k++)
					tem += local_Matrix_one[i * M + k] * Matrix_two[j * M + k];
				local_result[i * M + j] = tem;
			}
		free(local_Matrix_one);
		free(Matrix_two);
		//Сбор результатов
		MPI_Gather(local_result, local_M * N, MPI_INT, result_Matrix, local_M * N, MPI_INT, 0, MPI_COMM_WORLD);
		free(local_result);
		//printf("%d %d\n",local_M,my_rank);
	}
	MPI_Finalize();
	return 0;
}