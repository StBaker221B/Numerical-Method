
#include<stdio.h>

#define M 3
#define N 3

int LU(double a[M][N+1], double x[M], double L[M][N], double U[M][N], double y[M]);

int main()
{
    // printf("hello");

    // const int M = 5, N = 5;
    double coefficients[M][N];
    double augmentation[M][N+1];
    double L[M][N];
    double U[M][N];
    double x[M];
    double y[M];
    
    // set the coefficients
    for(int i = 0;i<M;i++)
    {
        printf("put in the coeff of line %d ", i);
        for(int j = 0;j<N;j++)
        {
            scanf("%lf", &augmentation[i][j]);
            printf("%lf \n",augmentation[i][j]);
        }
        // y[i] = 0;
        printf("put in the y of line %d ", i);
        scanf("%lf", &augmentation[i][N]);
    }
    // printf("initial over \n");
    
    // check
    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < N+1; j++)
        {
            printf("%f  ", augmentation[i][j]);
        }
        printf("\n");
    }

    // solve
    LU(augmentation, x, L, U, y);
    for(int i = 0; i < M; i++)
    {
        printf(" %f  ", x[i]);
    }
}

int LU(double a[M][N+1], double x[M], double L[M][N], double U[M][N], double y[M])
{
    //decompose
    for(int i = 0; i < M; i++)
    {
        
        // solve U by row
        for(int j = 0; j < N; j++)
        {
            if( i > j )
            {
                U[i][j] = 0;
                continue;
            }
            else if( i == 0)
            {
                U[i][j] = a[i][j];
                continue;
            }
            else
            {
                double t = 0;                
                for(int k = 0; k < i; k++)
                {
                    t += L[i][k] * U[k][j];
                }
                U[i][j] = a[i][j] - t;
                continue;
            }
        }

        // solve L by column
        for(int j = 0; j < N; j++)
        {
            if( j == i)
            {
                L[j][i] = 1;
                continue;
            }
            else if( j < i )
            {
                L[j][i] = 0;
                continue;
            }
            else if(i == 0)
            {
                L[j][i] = a[j][i] / U[i][i];
                continue;
            }
            else
            {
                double t = 0;
                for(int k = 0; k < i; k++)
                {
                    t += L[j][k] * U[k][i];
                }
                L[j][i] = (a[j][i] - t) / U[i][i];
                continue;
            }
        }
    }


    // for(int i = 0; i < M-1; i++)
    // {
    //     for(int k = i+1; k < M; k++ )
    //     {
    //         double t = a[k][i] / a[i][i];
            
    //         for(int j = i; j < N+1; j++)
    //         {
    //             a[k][j] = a[k][j] - a[i][j] * t;
    //         }
    //     }
    // }

    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < N; j++)
        {
            printf("%f  ", L[i][j]);
        }
        printf("\n");
    }

    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < N; j++)
        {
            printf("%f  ", U[i][j]);
        }
        printf("\n");
    }

/*
    // substitution
    for(int i = M-1; i >= 0; i--)
    {
        double t = 0;
        for(int j = N-1; j >= 0; j--)
        {
            if( i == (M-1) )
            {
                x[i] = a[i][N] / a[i][j];
                break;
            }

            else if( (j > 0) && a[i][j-1] == 0)
            {
                x[j] = ( a[i][N] - t ) / a[i][j];
                break;
            }
            else if( j == 0 )
            {
                x[j] = ( a[i][N] - t ) / a[i][j];
                break;
                
            }
            else 
            {
                t += x[j] * a[i][j];
            }
        }
    }
*/

    // solution
    // for y
    for(int i = 0; i < M; i++)
    {
        if(i == 0)
        {
            y[i] = a[i][N];
        }
        else
        {
            double t = 0;
            for(int k = 0; k < i; k++)
            {
                t += y[k] * L[i][k];
            }
            y[i] = a[i][N] - t;
        }
    }

    // for x
    for(int i = M-1; i >= 0; i--)
    {
        if(i == (M-1))
        {
            x[i] = y[i] / U[i][i];
        }
        else
        {
            double t = 0;
            for(int k = M-1; k > i; k--)
            {
                t += x[k] * U[i][k];
            }
            x[i] = (y[i] - t) / U[i][i];
        }        
    }

    return 0;
}