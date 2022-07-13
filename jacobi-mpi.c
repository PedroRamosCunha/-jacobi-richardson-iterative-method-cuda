//////////////////////////////////////////////////////////////////////
//                                                                  //
//      Trabalho 2 - SSC0903 - Computação de alto desempenho        //
//                                                                  //
//      10716550 - Diego da Silva Parra                             //
//      10691971 - Mateus Fernandes Doimo                           //
//      10892248 - Pedro Ramos Cunha                                //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define ROOT 0

void input(int rank, double ***matrizA, double **vetorA, double **vetorB, int N);
void verifica(int N, int no_procs, int rank);
void interacoes(double *x_new, double *x_old, double *x_bloc, double *a_recv, double *b_recv, int no_linha_blocos, int N, int rank);
void visualizar(int N, double **matrizA, double *vetorB, double *x, int rank);
void solucao(int N, double **matrizA, double *vetorB, double *x, int rank);
int check(double *x_old, double *x_new, int N);
double numRand();
double abs1(double x);

int T = 0;
int N = 0;
int P = 0;

int main(int argc, char *argv[])
{
    double tempo;
    tempo = MPI_Wtime();

    N = atoi(argv[1]);
    P = atoi(argv[2]);
    T = atoi(argv[3]);

    int no_linha_blocos;
    int no_procs, rank;
    double **matrizA, *vetorA, *vetorB, *a_recv, *b_recv;
    double *x_new, *x_old, *x_bloc;
    int provided;

    // Inicializar MPI
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &no_procs);

    input(rank, &matrizA, &vetorA, &vetorB, N);
    verifica(N, no_procs, rank);

    no_linha_blocos = N / no_procs;
    x_new = (double *)malloc(N * sizeof(double));
    x_old = (double *)malloc(N * sizeof(double));
    x_bloc = (double *)malloc(no_linha_blocos * sizeof(double));
    a_recv = (double *)malloc(no_linha_blocos * N * sizeof(double));
    b_recv = (double *)malloc(no_linha_blocos * N * sizeof(double));
    MPI_Scatter(vetorA, no_linha_blocos * N, MPI_DOUBLE, a_recv, no_linha_blocos * N, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Scatter(vetorB, no_linha_blocos, MPI_DOUBLE, b_recv, no_linha_blocos, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    interacoes(x_new, x_old, x_bloc, a_recv, b_recv, no_linha_blocos, N, rank);

    if (rank == ROOT)
    {
        tempo = MPI_Wtime() - tempo;
        printf("\n>> Tempo de execução: %f \n\n", tempo);

        solucao(N, matrizA, vetorB, x_new, rank);
    }

    MPI_Finalize();
    return 0;
}

// Entrada dos valores na matriz
void input(int rank, double ***matrizA, double **vetorA, double **vetorB, int N)
{
    int i, j;

    if (rank == ROOT)
    {

        *matrizA = (double **)malloc(N * sizeof(double *));
        *vetorB = (double *)malloc(N * sizeof(double));

        
        for (i = 0; i < N; i++)
        {
            (*vetorB)[i] = numRand();
            (*matrizA)[i] = (double *)malloc(N * sizeof(double));

            for (j = 0; j < N; j++)
            {
                (*matrizA)[i][j] = numRand();

                if (i == j){
                    (*matrizA)[i][j] = (*matrizA)[i][j] + N * 8;
                }
            }
        }

        *vetorA = (double *)malloc(N * N * sizeof(double));

        for (i = 0; i < N; i++)
        {
            #pragma omp parallel for num_threads(T)
            for (j = 0; j < N; j++)
            {
                (*vetorA)[(i*N) + j] = (*matrizA)[i][j];
            }
        }
    }
}

// Verifica número de processos com o tamanho da matriz
void verifica(int N, int no_procs, int rank)
{

    if (N % no_procs != 0)
    {
        MPI_Finalize();
        if (rank == ROOT)
        {
            printf(">> Erro: Número de linhas tem que ser divisível pelo número de processos!\n");
        }
        exit(-1);
    }
}

// Interações do método 
void interacoes(double *x_new, double *x_old, double *x_bloc, double *a_recv, 
                double *b_recv, int no_linha_blocos, int N, int rank)
{
    int i, j, global_rows_idx, iter = 0;

    #pragma omp parallel for num_threads(T)
    for (i = 0; i < no_linha_blocos; i++)
    {
        x_bloc[i] = b_recv[i];
    }

    MPI_Allgather(x_bloc, no_linha_blocos, MPI_DOUBLE, x_new, no_linha_blocos, MPI_DOUBLE, MPI_COMM_WORLD);
    
    do
    {
        #pragma omp parallel for num_threads(T)
        for (i = 0; i < N; i++)
        {
            x_old[i] = x_new[i];
        }

        for (i = 0; i < no_linha_blocos; i++)
        {
            global_rows_idx = (rank * no_linha_blocos) + i;
            x_bloc[i] = b_recv[i];

            for (j = 0; j < N; j++)
            {
                if (j == global_rows_idx)
                    continue;
                x_bloc[i] -= x_old[j] * a_recv[i * N + j];
            }
            x_bloc[i] /= a_recv[i * N + global_rows_idx];
        }

        // Receber todos os x_bloc de todos os processos 
        MPI_Allgather(x_bloc, no_linha_blocos, MPI_DOUBLE, x_new, no_linha_blocos, MPI_DOUBLE, MPI_COMM_WORLD);
        iter++;
    } while (check(x_old, x_new, N));
}

// Verifica se chegou no critério de parada
int check(double *x_old, double *x_new, int N)
{
    int i = 0;
    double max1 = 0.0;
    double max2 = 0.0;
    double m = 0.0;

    for (i = 0; i < N; i++)
    {
        if (abs1(x_new[i] - x_old[i]) >= max1)
        {
            max1 = abs1((x_new[i] - x_old[i]));
        }

        if (abs1(x_new[i]) >= max2)
        {
            max2 = abs1((x_new[i]));
        }
    }

    m = max1 / max2;

    if (m > 0.001)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

// Imprimir equações
void visualizar(int N, double **matrizA, double *vetorB,
             double *x, int rank)
{
    int i, j;
    printf(">> Método Iterativo de Jacobi-Richardson\n\n");
    printf(">> Equações: \n");
    for (i = 0; i < N; i++)
    {
        for (j = 0; j <= N; j++)
        {
            if (j == N - 1)
                printf("%0.2f * X%d = ", matrizA[i][j], j);
            else if (matrizA[i][j] < 0)
                printf("%0.2f * X%d", matrizA[i][j], j);
            else if (j == N)
                printf("%0.2f", vetorB[i]);
            else
                printf("%0.2f * X%d + ", matrizA[i][j], j);
        }
        printf("\n");
    }
    printf("\n");

    printf("Solução: \n");
    for (i = 0; i < N; i++)
    {
        printf("x[%d] = %lf\n", i, x[i]);
    }
}

// Verificar solução interativa
void solucao(int N, double **matrizA, double *vetorB,
             double *x, int rank)
{
    int linha = 0;
    float resultado = 0.0;

    printf("\n>> Escolha uma linha para resolver a equação: \n");
    while (1)
    {
        scanf("%d", &linha);

        if (linha <= N && linha >= 1)
            break;

        printf("\n>> Escolha uma linha para resolver a equação: \n");
    }

    linha--;

    printf("\n>> Cálculado pela interação: \n");
    for (int j = 0; j < N; j++)
    {
        resultado += matrizA[linha][j] * x[j];
    }

    printf(">> Resultado encontrado: %0.4lf\n", resultado);

    printf("\n>> Resposta esperada: %0.4lf\n", vetorB[linha]);
}

// Valores pseudoaleatória
double numRand()
{
    return ((rand() % 8) * (1.0));
}

// Retorna valor absoluto
double abs1(double x)
{
    if (x > 0.0)
        return x;
    else
        return x * (-1.0);
}