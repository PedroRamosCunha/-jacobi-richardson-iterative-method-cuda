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
#include <string.h>
#include <omp.h>

float *varprev;
float *vetorB;
float *varcurr;
float **matrizA;

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

void initialise(int N)
{
    int i = 0;
    for (i = 0; i < N; i++)
    {
        varprev[i] = 0.0;
        varcurr[i] = 0.0;
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
    printf("\n>> Equção de entrada:\n");
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

void solucao(int N)
{
    int i = 0;
    printf("\n>> Solução: \n");
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

    while (1)
    {
        for (i = 0; i < N; i++)
        {
            varcurr[i] = vetorB[i]; // vetor B
            for (j = 0; j < N; j++)
            {
                if (i != j)
                {
                    varcurr[i] = varcurr[i] + reverse(matrizA[i][j] * varprev[j]);
                }
            }
            varcurr[i] = varcurr[i] / matrizA[i][i];
        }
        if (check(N) == 1)
        {
            break;
        }
        else
        {
            for (i = 0; i < N; i++)
                varprev[i] = varcurr[i];
        }
        interacao++;
    }
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
    }

    printf(">> Resultado encontrado: %0.4f\n", resultado);

    printf("\n>> Resposta esperada: %0.2f\n", vetorB[linha]);
}

int main(int argc, char **argv)
{
    double t;
    int N = atoi(argv[1]);
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

    initialise(N);
    input(N);
    solve(N);
    t = omp_get_wtime() - t;
    printf("\n>> Tempo de execução: %f \n\n", t);
    resolve(N);
    return 0;
}
