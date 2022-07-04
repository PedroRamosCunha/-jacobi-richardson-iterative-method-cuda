#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

float *varprev;
float *vetorB;
float *varcurr;
float **matrizA;

int T = 0;

// Valores pseudoaleatória
float numRand()
{
    return ((rand() % 8) * (1.0));
}

// Retorna valor absoluto
float abs1(float x)
{
    if (x > 0.0)
        return x;
    else
        return x * (-1.0);
}

// Condição de convergência, se a matriz A for diagonalmente dominante
void diagonalmenteDominante(int N)
{

    int i, j = 0;
    float somaNaoDiagonalLinha, somaNaoDiagonalColuna = 0.0;
    int flag1 = 0;
    int flag2 = 0;

#pragma omp parallel for num_threads(T)
    for (i = 0; i < N; i++)
    {
        somaNaoDiagonalLinha = 0.0;
        somaNaoDiagonalColuna = 0.0;

#pragma omp parallel for num_threads(T) reduction(+                                   \
                                                  : somaNaoDiagonalLinha) reduction(+ \
                                                                                    : somaNaoDiagonalColuna)
        for (j = 0; j < N; j++)
        {
            if (i != j)
            {
                somaNaoDiagonalLinha += abs1(matrizA[i][j]);
                somaNaoDiagonalColuna += abs1(matrizA[j][i]);
            }
        }

        if (matrizA[i][i] < somaNaoDiagonalLinha)
        {
            flag1 = 1;
        }
        if (matrizA[i][i] < somaNaoDiagonalColuna)
        {
            flag2 = 1;
        }
    }

    if (flag1 && flag2)
    {
        printf("\n>> Não é posssível garantir a convergência da solução!\n");
        exit(0);
    }

    return;
}

void inputManual(int N)
{
    int i, j = 0;
    for (i = 0; i < N; i++)
    {
        printf("\nEnter coefficients for equation %d:\n", (i + 1));
        for (j = 0; j <= N; j++)
        {
            if (j == N)
            {
                printf("\nEnter result value for equation %d:\n", (i + 1));
                scanf("%f", &vetorB[i]);
            }
            else
            {
                scanf("%f", &matrizA[i][j]);
            }
        }
    }
}

void input(int N)
{
    int i, j = 0;

    for (i = 0; i < N; i++)
    {
        vetorB[i] = numRand();

        for (j = 0; j < N; j++)
        {

            matrizA[i][j] = numRand();

            if (i == j)
            {
                matrizA[i][j] = matrizA[i][j] + N * 8;
            }
        }
    }
}

void preview(int N)
{
    int i, j = 0;
    printf("\nEqução de entrada:\n");
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
}

int check(int N)
{
    int i = 0;
    float sum = 0.0;
    float max1 = 0.0;
    float max2 = 0.0;
    float m = 0.0;

#pragma omp parallel for num_threads(T) private(i) reduction(max                   \
                                                             : max1) reduction(max \
                                                                               : max2)
    for (i = 0; i < N; i++)
    {
        if (abs1(varcurr[i] - varprev[i]) >= max1)
        {
            max1 = abs1((varcurr[i] - varprev[i]));
        }

        if (abs1(varcurr[i]) >= max2)
        {
            max2 = abs1((varcurr[i]));
        }
    }

    m = max1 / max2;

    if (m <= 0.001)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void show_solution(int N)
{
    int i = 0;
    printf("\nSolução: \n");
    for (i = 0; i < N; i++)
    {
        printf("X%d = %0.3f\n", i, varcurr[i]);
    }
}

float reverse(float x)
{
    return x * (-1);
}

void solve(int N)
{
    int i, j, interacao = 0;

    // while (1)
    // {

    //     for (i = 0; i < N; i++)
    //     {
    //         varcurr[i] = vetorB[i];
    //         for (j = 0; j < N; j++)
    //         {
    //             if (i != j)
    //             {
    //                 varcurr[i] = varcurr[i] + reverse(matrizA[i][j] * varprev[j]);
    //             }
    //         }
    //         varcurr[i] = varcurr[i] / matrizA[i][i];
    //     }
    //     if (check(N) == 1)
    //     {
    //         show_solution(N);
    //         break;
    //     }
    //     else
    //     {
    //         for (i = 0; i < N; i++)
    //             varprev[i] = varcurr[i];
    //     }
    //     interacao++;
    // }

    //     while (1)
    //     {

    // #pragma omp parallel for reduction(+ \
//                                    : varcurr[:N]) num_threads(T) private(i, j)
    //         for (i = 0; i < N; i++)
    //         {
    //             varcurr[i] = vetorB[i];
    //             for (j = 0; j < N; j++)
    //             {
    //                 if (i != j)
    //                 {
    //                     varcurr[i] = varcurr[i] + reverse(matrizA[i][j] * varprev[j]);
    //                 }
    //             }
    //             varcurr[i] = varcurr[i] / matrizA[i][i];
    //         }
    //         if (check(N) == 1)
    //         {
    //             show_solution(N);
    //             break;
    //         }
    //         else
    //         {
    //             for (i = 0; i < N; i++)
    //                 varprev[i] = varcurr[i];
    //         }
    //         interacao++;
    //     }

    // #pragma omp parallel
    //     {

    float somatudo = 0.0;

    do
    {
        //#pragma omp parallel for private(i)
        for (i = 0; i < N; i++)
        {

            varprev[i] = varcurr[i];
            // varcurr[i] = vetorB[i];

#pragma omp parallel for private(j) num_threads(T) reduction(+ \
                                                             : varcurr[:N])
            for (j = 0; j < N; j++)
            {
                if (i != j)
                {
                    varcurr[i] += reverse(matrizA[i][j] * varprev[j]);
                }
            }
            /////

            varcurr[i] += vetorB[i];
            varcurr[i] = varcurr[i] / matrizA[i][i];
        }

        interacao++;
    } while (!check(N));
    //}

    // show_solution(N);
}

void resolve(int N)
{
    int linha = 0;
    float resultado = 0.0;

    while (1)
    {
        printf("\n>> Escolha uma linha para resolver a equação: ");
        scanf("%d", &linha);

        if (linha <= N && linha >= 1)
            break;
    }

    linha--;

    printf(">> Cálculado pela interação: \n");
    for (int j = 0; j < N; j++)
    {
        resultado += matrizA[linha][j] * varcurr[j];
        // if (j == N - 1)
        //     printf("%0.2f * %0.3f = ", matrizA[linha][j], varcurr[j]);
        // else
        //     printf("%0.2f * %0.3f + ", matrizA[linha][j], varcurr[j]);
    }

    printf(">> Resposta encontrada: %0.4f\n", resultado);

    printf("\n>> Resposta esperada: %0.2f\n", vetorB[linha]);
}

int main(int argc, char **argv)
{
    double t;
    int N = atoi(argv[1]);
    T = atoi(argv[2]);

    omp_set_num_threads(T);

    printf("T = %d", T);
    t = omp_get_wtime();
    matrizA = (float **)calloc(N, sizeof(float *));

    if (matrizA == NULL)
    {
        printf("** Erro: Memoria Insuficiente **");
        exit(0);
    }

    for (int i = 0; i < N; i++)
    {
        matrizA[i] = (float *)calloc(N, sizeof(float));
        if (matrizA[i] == NULL)
        {
            printf("** Erro: Memoria Insuficiente **");
            exit(0);
        }
    }

    vetorB = (float *)calloc(N, sizeof(float));
    if (vetorB == NULL)
    {
        printf("** Erro: Memoria Insuficiente **");
        exit(0);
    }
    varcurr = (float *)calloc(N, sizeof(float));
    if (varcurr == NULL)
    {
        printf("** Erro: Memoria Insuficiente **");
        exit(0);
    }
    varprev = (float *)calloc(N, sizeof(float));
    if (varprev == NULL)
    {
        printf("** Erro: Memoria Insuficiente **");
        exit(0);
    }

    // Semente para o número aleatório
    srand(0);

    printf("N = %d", N);

    input(N);
    // preview(N);
    diagonalmenteDominante(N);
    solve(N);
    t = omp_get_wtime() - t;
    resolve(N);

    printf("\n>> Tempo de execução: %f \n", t);

    return 0;
}
